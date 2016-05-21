#!/usr/bin/python

"""
Create database structure for iDR5 quality control (QC) phase.
"""

import logging
import numpy as np
import psycopg2 as pg
import yaml

from code import GESDatabase

db_filename = "code/db.yaml"
nodes_filename = "nodes.yaml"
schema_filename = "code/schema.sql"
masterlist_filename \
    = "fits-templates/masterlist/GES_iDR5_spectra_masterlist_15052016.fits"

# Connect to database.
with open(db_filename, "r") as fp:
    credentials = yaml.load(fp)
    connection = pg.connect(**credentials)

logger = logging.getLogger("ges")
logger.info("Connected to database.")


# Create the tables.
cursor = connection.cursor()
logger.info("Creating tables from {}...".format(schema_filename))
with open(schema_filename, "r") as fp:
    cursor.execute(fp.read())
cursor.close()
logger.info("Tables created.")

connection.commit()
connection.close()


# Create a database object.
database = GESDatabase(**credentials)

# Create nodes.
with open(nodes_filename, "r") as fp:
    all_nodes = yaml.load(fp)

for wg, node_names in all_nodes.items():
    for node_name in node_names:
        node_id = database.create_or_retrieve_node_id(wg, node_name)

# Ingest the masterlist of spectra.
N_ingested = database.ingest_spectra_masterlist(masterlist_filename)

# Ingest faux results from two different nodes.


# Make plots.

