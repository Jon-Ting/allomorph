"""Merge and reformat NCPac-extracted features with MD outputs."""

import csv
import logging
import math
import os
from multiprocessing import Pool
from os.path import exists

import pandas as pd

from np_gen.constants import ELE_DICT, parse_ele_comb

logger = logging.getLogger(__name__)

N_AVOGADRO = 6.02214076 * 10 ** 23  # (atoms/mol)
A3_PER_M3 = 10 ** 30  # Angstrom^3 per m^3 (dimensionless)


def _get_ele_prop_dict(ele_dict=None):
    """Return a property dictionary compatible with the feature-extraction pipeline.

    Falls back to the global ELE_DICT if none is provided.
    """
    if ele_dict is None:
        ele_dict = ELE_DICT
    prop_dict = {}
    for ele, props in ele_dict.items():
        # Ensure the required bulk-packing properties are present
        rho = props.get("rho")
        m = props.get("m")
        bulkE = props.get("bulkE")
        if rho is None or m is None or bulkE is None:
            raise ValueError(
                f"Element '{ele}' is missing one or more of 'rho', 'm', 'bulkE' "
                f"required by the feature-extraction pipeline. "
                f"Please extend ELE_DICT with these keys or supply a custom ele_dict."
            )
        prop_dict[ele] = {"rho": rho, "m": m, "bulkE": bulkE}
    return prop_dict


def generate_headers(elements=None):
    """Programmatically generate the full list of features."""
    if elements is None:
        elements = sorted(ELE_DICT.keys())
    headers = [
        "T", "P", "Potential_E", "Kinetic_E", "Total_E",
        "N_atom_total"
    ]
    for ele in elements:
        headers.append(f"N_{ele}")
    headers.extend([
        "N_atom_bulk", "N_atom_surface",
        "Vol_bulk_pack", "Vol_sphere",
        "R_min", "R_max", "R_diff", "R_avg", "R_std", "R_skew", "R_kurt",
        "S_100", "S_111", "S_110", "S_311"
    ])

    # Curvature features
    for i in range(1, 19):
        headers.append(f"Curve_{(i-1)*10+1}-{i*10}")

    # Elemental Coordination Numbers (MM and each element)
    combos = ["MM"] + [f"{ele}M" for ele in elements]
    cn_types = ["TCN", "BCN", "SCN", "SOCN", "TGCN", "BGCN", "SGCN", "SOGCN"]
    for combo in combos:
        for cn in cn_types:
            headers.append(f"{combo}_{cn}_avg")
            for i in range(21):
                headers.append(f"{combo}_{cn}_{i}")

    # Bond lengths
    headers.extend(["MM_BL_avg", "MM_BL_std", "MM_BL_max", "MM_BL_min", "MM_BL_num"])
    for i, e1 in enumerate(elements):
        for e2 in elements[i:]:
            headers.extend([f"{e1}{e2}_BL_avg", f"{e1}{e2}_BL_std", f"{e1}{e2}_BL_max", f"{e1}{e2}_BL_min", f"{e1}{e2}_BL_num"])

    # Bond fractions
    for i, e1 in enumerate(elements):
        for e2 in elements[i:]:
            headers.append(f"{e1}{e2}_frac")
    headers.append("N_bond")

    # Bond Angles (BA1 and BA2)
    ba_combos = ["MMM"]
    for e1 in elements:
        for e2 in elements:
            for e3 in elements:
                ba_combos.append(f"{e1}{e2}{e3}")

    for ba in ["BA1", "BA2"]:
        for combo in ba_combos:
            headers.extend([f"{combo}_{ba}_avg", f"{combo}_{ba}_std", f"{combo}_{ba}_max", f"{combo}_{ba}_min", f"{combo}_{ba}_num"])

    # Bond Torsion
    headers.extend(["MMMM_BTneg_avg", "MMMM_BTneg_std", "MMMM_BTneg_max", "MMMM_BTneg_min", "MMMM_BTneg_num",
                    "MMMM_BTpos_avg", "MMMM_BTpos_std", "MMMM_BTpos_max", "MMMM_BTpos_min", "MMMM_BTpos_num"])

    bt_combos = []
    for e1 in elements:
        for e2 in elements:
            for e3 in elements:
                for e4 in elements:
                    bt_combos.append(f"{e1}{e2}{e3}{e4}")

    for combo in bt_combos:
        headers.extend([f"{combo}_BTneg_avg", f"{combo}_BTneg_std", f"{combo}_BTneg_max", f"{combo}_BTneg_min", f"{combo}_BTneg_num",
                        f"{combo}_BTpos_avg", f"{combo}_BTpos_std", f"{combo}_BTpos_max", f"{combo}_BTpos_min", f"{combo}_BTpos_num"])

    # Steinhardt order parameters
    headers.extend(["q6q6_T_avg", "q6q6_B_avg", "q6q6_S_avg"])
    for space in ["T", "B", "S"]:
        for i in range(21):
            headers.append(f"q6q6_{space}_{i}")
        headers.append(f"q6q6_{space}_20+")

    # Structure classification
    headers.extend(["FCC", "HCP", "ICOS", "DECA"])

    # Surface characteristics
    for characteristic in ["defects", "micros", "facets"]:
        for ele in list(elements) + [""]:
            suffix = f"_{ele}" if ele else ""
            headers.extend([
                f"Surf_{characteristic}{suffix}",
                f"Surf_{characteristic}{suffix}_bulk_pack_conc",
                f"Surf_{characteristic}{suffix}_bulk_pack_ratio",
                f"Surf_{characteristic}{suffix}_sphere_conc",
                f"Surf_{characteristic}{suffix}_sphere_ratio"
            ])

    return headers



def _build_add_feat_list(elements=None):
    """Build the list of derived features that will be re-added after dropping."""
    if elements is None:
        elements = sorted(ELE_DICT.keys())
    add_feat_list = [
        "Vol_bulk_pack", "Vol_sphere",
        "Curve_1-10", "Curve_11-20", "Curve_21-30", "Curve_31-40", "Curve_41-50",
        "Curve_51-60", "Curve_61-70", "Curve_71-80", "Curve_81-90", "Curve_91-100",
        "Curve_101-110", "Curve_111-120", "Curve_121-130", "Curve_131-140",
        "Curve_141-150", "Curve_151-160", "Curve_161-170", "Curve_171-180",
    ]
    for characteristic in ["defects", "micros", "facets"]:
        for ele in list(elements) + [""]:
            suffix = f"_{ele}" if ele else ""
            add_feat_list.extend([
                f"Surf_{characteristic}{suffix}",
                f"Surf_{characteristic}{suffix}_bulk_pack_conc",
                f"Surf_{characteristic}{suffix}_bulk_pack_ratio",
                f"Surf_{characteristic}{suffix}_sphere_conc",
                f"Surf_{characteristic}{suffix}_sphere_ratio"
            ])
    return add_feat_list


def _calc_bulk_pack_vol(row, ele_comb, ele_dict=None):
    elements = parse_ele_comb(ele_comb, ele_dict)
    ele_prop_dict = _get_ele_prop_dict(ele_dict)
    total_vol = 0
    for element in elements:
        if f"N_{element}" not in row:
            continue
        vol_ele = (
            row[f"N_{element}"]
            * (ele_prop_dict[element]["m"] / N_AVOGADRO)
            / (ele_prop_dict[element]["rho"] / A3_PER_M3)
        )
        total_vol += vol_ele
    return total_vol


def _cnt_surf_site(row, site_type, element, ele_comb, ele_dict=None):
    if len(element) > 2 or element == "":
        elements = parse_ele_comb(ele_comb, ele_dict)
        return sum([row.get(f"Surf_{site_type}_{ele}", 0) for ele in elements])

    if site_type == "defects":
        return sum(
            (
                row.get(f"{element}M_SCN_1", 0),
                row.get(f"{element}M_SCN_2", 0),
                row.get(f"{element}M_SCN_3", 0),
            )
        )
    elif site_type == "micros":
        return sum(
            (
                row.get(f"{element}M_SCN_4", 0),
                row.get(f"{element}M_SCN_5", 0),
                row.get(f"{element}M_SCN_6", 0),
                row.get(f"{element}M_SCN_7", 0),
            )
        )
    elif site_type == "facets":
        return sum(
            (
                row.get(f"{element}M_SCN_9", 0),
                row.get(f"{element}M_SCN_10", 0),
                row.get(f"{element}M_SCN_11", 0),
            )
        )
    else:
        raise ValueError(f"Unknown site type: {site_type}")


def _calc_surf_site_conc(row, site_type, element, vol_type, ele_comb, ele_dict=None):
    if len(element) > 2 or element == "":
        elements = parse_ele_comb(ele_comb, ele_dict)
        return sum(
            [row.get(f"Surf_{site_type}_{ele}_{vol_type}_conc", 0) for ele in elements]
        )
    vol = row.get(f"Vol_{vol_type}", 0)
    if vol == 0:
        return 0
    return row.get(f"Surf_{site_type}_{element}", 0) / vol


def _calc_surf_site_ratio(row, site_type, element, vol_type, ele_comb, ele_dict=None):
    if len(element) > 2 or element == "":
        elements = parse_ele_comb(ele_comb, ele_dict)
        return sum(
            [row.get(f"Surf_{site_type}_{ele}_{vol_type}_ratio", 0) for ele in elements]
        )
    total_atoms = row.get("N_atom_total", 0)
    if total_atoms == 0:
        return 0
    return row.get(f"Surf_{site_type}_{element}", 0) / total_atoms


def _drop_feats(df, all_headers, add_feat_list):
    """Drop unused features using programmed offsets."""
    logger.debug("Dropping unused features...")

    # Offsets in the original NCPac output (relative to start of NCPac columns)
    CURV_START = 22
    RDF_START = 969
    SF_START = 1821
    BA_START = 2334
    BT_START = 12876
    CL_START = 14372

    # Remove re-added features from header list
    for col in add_feat_list:
        if col in all_headers:
            all_headers.remove(col)

    drop_indices = [5]  # Frame index

    # Add temporary curvature columns
    for i in range(180):
        all_headers.insert(CURV_START + i, f"Curve_Deg_{i+1}")

    drop_indices.extend(range(RDF_START, RDF_START + 852))
    drop_indices.extend(range(SF_START, SF_START + 513))

    ba_deg_cols = []
    for i in range(10528):
        col_idx = BA_START + i
        if (i % 188) >= 8:
            ba_deg_cols.append(col_idx)
    drop_indices.extend(ba_deg_cols)

    drop_indices.extend(range(BT_START, BT_START + 362))
    drop_indices.extend(range(CL_START, CL_START + 20))

    df.drop(df.columns[drop_indices], axis=1, inplace=True)
    df = df[df.columns.drop(list(df.filter(regex="Type")))]

    df_numeric = df.apply(pd.to_numeric, errors="coerce")
    nan_count = df_numeric.isna().sum().sum()
    if nan_count > 0:
        logger.warning(f"Found {nan_count} NaN values after numeric conversion.")

    df_numeric.columns = all_headers
    return df_numeric


def _add_feats(df, ele_comb, ele_dict=None):
    """Add derived features."""
    logger.debug("Adding derived features...")
    df["Vol_bulk_pack"] = df.apply(lambda row: _calc_bulk_pack_vol(row, ele_comb, ele_dict), axis=1)
    df["Vol_sphere"] = df.apply(lambda row: 3 / 4 * math.pi * row["R_avg"] ** 3, axis=1)

    # Bin curvature
    curv_col_idx = 22
    for i in range(1, 19):
        curv_col_name = f"Curve_{(i-1)*10+1}-{i*10}"
        df[curv_col_name] = df.iloc[:, curv_col_idx : curv_col_idx + 10].sum(axis=1)
        curv_col_idx += 10

    df.drop(df.columns[range(22, 22 + 180)], axis=1, inplace=True)

    elements = parse_ele_comb(ele_comb, ele_dict)
    for characteristic in ("defects", "micros", "facets"):
        for element in elements:
            df[f"Surf_{characteristic}_{element}"] = df.apply(
                lambda row: _cnt_surf_site(row, characteristic, element, ele_comb, ele_dict), axis=1
            )
            for vol_type in ("bulk_pack", "sphere"):
                df[f"Surf_{characteristic}_{element}_{vol_type}_conc"] = df.apply(
                    lambda row: _calc_surf_site_conc(row, characteristic, element, vol_type, ele_comb, ele_dict),
                    axis=1,
                )
                df[f"Surf_{characteristic}_{element}_{vol_type}_ratio"] = df.apply(
                    lambda row: _calc_surf_site_ratio(row, characteristic, element, vol_type, ele_comb, ele_dict),
                    axis=1,
                )
        df[f"Surf_{characteristic}"] = df.apply(
            lambda row: _cnt_surf_site(row, characteristic, "", ele_comb, ele_dict), axis=1
        )
        for vol_type in ("bulk_pack", "sphere"):
            df[f"Surf_{characteristic}_{vol_type}_conc"] = df.apply(
                lambda row: _calc_surf_site_conc(row, characteristic, "", vol_type, ele_comb, ele_dict),
                axis=1,
            )
            df[f"Surf_{characteristic}_{vol_type}_ratio"] = df.apply(
                lambda row: _calc_surf_site_ratio(row, characteristic, "", vol_type, ele_comb, ele_dict),
                axis=1,
            )
    return df


def merge_reformat_data(output_md, feat_source_path, feat_eng_path, ele_comb, all_headers_list, verbose=True, ele_dict=None):
    """Concatenate MD output with features extracted from NCPac."""
    conf_id = output_md[0]
    if verbose:
        logger.info(f"Concatenating CSV files for nanoparticle {conf_id}...")

    df1 = pd.DataFrame([output_md])
    df1.columns = ["confID", "T", "P", "Potential_E", "Kinetic_E", "Total_E"]

    feat_file = os.path.join(feat_source_path, f"{conf_id}.csv")
    if not os.path.exists(feat_file) or os.path.getsize(feat_file) == 0:
        logger.warning(f"{conf_id}.csv not found or empty! Skipping...")
        return

    df2 = pd.read_csv(feat_file, sep=",", header=1, index_col=None)
    df = pd.concat([df1, df2], axis="columns")
    df.set_index(keys="confID", inplace=True)

    elements = parse_ele_comb(ele_comb, ele_dict)
    add_feat_list = _build_add_feat_list(elements)
    df = _drop_feats(df, all_headers_list.copy(), add_feat_list)
    df = _add_feats(df, ele_comb, ele_dict)

    energy_cols = ["Potential_E", "Kinetic_E", "Total_E"]
    other_cols = [c for c in df.columns if c not in energy_cols]
    df = df[other_cols + energy_cols]

    output_file = os.path.join(feat_eng_path, f"{conf_id}.csv")
    os.makedirs(feat_eng_path, exist_ok=True)
    df.to_csv(output_file, sep=",", header=True)


def run_merge_reformat_parallel(md_out_fname, feat_source_path, feat_eng_path, ele_comb, all_headers_list, verbose=False, ele_dict=None):
    if verbose:
        logger.info("Merging information and reformating data in parallel...")

    output_mds = []
    with open(md_out_fname, "r") as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            conf_id = row[0]
            if not exists(os.path.join(feat_eng_path, f"{conf_id}.csv")):
                output_mds.append(row)

    if not output_mds:
        logger.info("No new nanoparticles to merge.")
        return

    with Pool() as p:
        p.starmap(
            merge_reformat_data,
            [
                (om, feat_source_path, feat_eng_path, ele_comb, all_headers_list, verbose, ele_dict)
                for om in output_mds
            ],
        )


def concat_np_feats(feat_eng_path, final_data_fname, verbose=False):
    """Concatenate all processed feature CSV files into one final dataset."""
    if verbose:
        logger.info("Concatenating processed feature CSV files...")

    feat_csvs = sorted([f for f in os.listdir(feat_eng_path) if f.endswith(".csv")])
    if not feat_csvs:
        logger.warning("No CSV files found to concatenate.")
        return

    with open(final_data_fname, "wb") as fout:
        for i, feat_csv in enumerate(feat_csvs):
            feat_file = os.path.join(feat_eng_path, feat_csv)
            if os.path.getsize(feat_file) == 0:
                continue
            with open(feat_file, "rb") as f:
                if i != 0:
                    next(f)
                fout.write(f.read())
    logger.info(f"Final dataset saved to {final_data_fname}")
