"""EAM (Embedded Atom Method) potential file generation."""

from np_gen.eam.create_eam import create_eam
from np_gen.eam.eam_database import AtType, Database

__all__ = ["AtType", "Database", "create_eam"]
