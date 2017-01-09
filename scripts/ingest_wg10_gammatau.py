
"""
Ingest \gamma and \tau values from the second extension of the WG10 recommended
file (per setup).
"""

import logging
import numpy as np
import psycopg2 as pg
import yaml
from glob import glob
from astropy.io import fits


from code import GESDatabase

logger = logging.getLogger("ges")

# Get credentials.
db_filename = "db.yaml"
with open(db_filename, "r") as fp:
    credentials = yaml.load(fp)

# Create a database object.
database = GESDatabase(**credentials)

image = fits.open("recommended-results/GES_iDR5_WG10_Recommended_PERSETUP.fits")
data = image[2].data
N = len(data)


# Update records in the database for the given CNAME/FILENAME
N = len(data)
for i, row in enumerate(data):

    logger.info("Updating row {}/{}".format(i, N))
    result = database.execute(
        """ UPDATE  wg_recommended_results
            SET     gamma = %(gamma)s,
                    e_gamma = %(e_gamma)s,
                    tau = %(tau)s,
                    e_tau = %(e_tau)s
            WHERE   wg = 10
              AND   cname = %(cname)s
              AND   TRIM(filename) = %(filename)s
            LIMIT   1
        """, 
        dict(
            gamma=row["GAM_1"],
            e_gamma=row["GAM_1_ERR"],
            tau=row["TAU"],
            e_tau=row["TAU_ERR"],
            cname=row["CNAME"],
            filename=row["FILENAME"]))
    logger.info(result)

database.connection.commit()
   