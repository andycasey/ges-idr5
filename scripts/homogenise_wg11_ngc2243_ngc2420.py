
"""
Do WG11 homogenisation for NGC2243 and NGC2420.
"""

import yaml
import logging
import numpy as np

from code import GESDatabase
from code.model.ensemble import EnsembleModel

# Initialize logging.
logger = logging.getLogger("ges")

# Create a database object.
db_filename = "db.yaml"
with open(db_filename, "r") as fp:
    credentials = yaml.load(fp)
database = GESDatabase(**credentials)

wg = 11
parameters = ("teff", "logg", "feh")

for parameter in parameters:

    model_path = "homogenisation-wg{wg}-{parameter}-poly-sys.model".format(
        wg=wg, parameter=parameter)

    model = EnsembleModel.read(model_path, database)
        
    # Homogenise this parameter for all stars in NGC2243
    model.homogenise_stars_matching_query(
        """ WITH n as (SELECT id FROM nodes WHERE wg = 11)
            SELECT DISTINCT ON (r.cname) r.cname, s.ges_fld
            FROM n, spectra as s, results as r
            WHERE r.node_id = n.id 
              AND s.cname = r.cname
              AND s.ges_fld LIKE 'NGC2243%'
            ORDER BY cname DESC""",
        update_database=True)

    # Homogenise this parameter for all stars in NGC2420
    model.homogenise_stars_matching_query(
        """ WITH n as (SELECT id FROM nodes WHERE wg = 11)
            SELECT DISTINCT ON (r.cname) r.cname, s.ges_fld
            FROM n, spectra as s, results as r
            WHERE r.node_id = n.id 
              AND s.cname = r.cname
              AND s.ges_fld LIKE 'NGC2420%'
            ORDER BY cname DESC""",
        update_database=True)


