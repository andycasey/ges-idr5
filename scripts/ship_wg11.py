
"""
Ship a WG11 recommended SP product.
"""

import yaml
import logging

from code import GESDatabase, ship

# Initialize logging.
logger = logging.getLogger("ges")

# Create a database object.
db_filename = "db.yaml"
with open(db_filename, "r") as fp:
    credentials = yaml.load(fp)
database = GESDatabase(**credentials)


# Produce a PER-SETUP file even though they will be duplicates.
ship.wg_recommended_sp_template(database,
    "fits-templates/recommended-templates/GES_iDR5_WG11_RecommendedTemplate_16072016_PERSETUP_FixedLengthFLAGColumns.fits",
    "outputs/GES_iDR5_WG11_Recommended_PERSETUP.fits",
    11, overwrite=True)

raise a

# Produce a per CNAME file.
ship.wg_recommended_sp_template(database,
    "fits-templates/recommended-templates/GES_iDR5_WG11_RecommendedTemplate_16072016_FixedLengthFLAGColumns.fits",
    "outputs/GES_iDR5_WG11_Recommended.fits",
    11, overwrite=True)


