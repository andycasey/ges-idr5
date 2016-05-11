#!/usr/bin/python

"""
Create database structure for iDR5 quality control (QC) phase.
"""

import logging
import psycopg2 as pg
import yaml
#from astropy.io import fits


db_filename = "db.yaml"
schema_filename = "schema.sql"


# Connect to database.
with open(db_filename, "r") as fp:
    connection = pg.connect(**yaml.load(fp))
logging.info("Connected to database.")


# Create the tables.
cursor = connection.cursor()
logging.info("Creating tables from {}...".format(schema_filename))
with open(schema_filename, "r") as fp:
    cursor.execute(fp.read())
cursor.close()
logging.info("Tables created.")

# Ingest from the master-list.
# TODO:



connection.commit()


