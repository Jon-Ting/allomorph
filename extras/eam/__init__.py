"""EAM (Embedded Atom Method) potential file generation."""

from eam.create_eam import create_eam
from eam.eam_database import AtType, Database

__all__ = ["AtType", "Database", "create_eam"]
