"""EAM (Embedded Atom Method) potential file generation."""

from allomorph.eam.create_eam import create_eam
from allomorph.eam.eam_database import AtType, Database

__all__ = ["AtType", "Database", "create_eam"]
