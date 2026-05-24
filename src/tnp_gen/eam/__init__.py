"""EAM (Embedded Atom Method) potential file generation."""

from tnp_gen.eam.create_eam import create_eam
from tnp_gen.eam.eam_database import AtType, Database

__all__ = ["AtType", "Database", "create_eam"]
