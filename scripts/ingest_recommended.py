#!/usr/bin/python

"""
Ingest the WG-level recommended results.
"""

import logging
import numpy as np
import psycopg2 as pg
import yaml
from glob import glob

from code import GESDatabase

logger = logging.getLogger("ges")

# Get credentials.
db_filename = "db.yaml"
with open(db_filename, "r") as fp:
    credentials = yaml.load(fp)

# Create a database object.
database = GESDatabase(**credentials)


# Ingest results 
filenames = glob("recommended-results/GES_iDR5_WG??_Recommended.fits")
for filename in filenames:
    logger.info("Ingesting from {}".format(filename))
    database.ingest_recommended_results(filename, extension=1)
