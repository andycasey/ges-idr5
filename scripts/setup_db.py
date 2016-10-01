#!/usr/bin/python

"""
Create database for the fifth internal data release of the Gaia-ESO Survey.
"""

import logging
import numpy as np
import psycopg2 as pg
import yaml
from glob import glob


# For fake data generation
import os
from astropy.io import fits

from code import GESDatabase


db_filename = "db.yaml"
nodes_filename = "nodes.yaml"
schema_filename = "code/schema.sql"
masterlist_filename \
    = "fits-templates/masterlist/GES_iDR5_spectra_masterlist_15052016.fits"

# Connect to database.
with open(db_filename, "r") as fp:
    credentials = yaml.load(fp)
    credentials["database"] = "ges_idr5_with_spectra_snr"
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

# Ingest results from the nodes.
for filename in glob("node-results/stellar-parameters/WG??/GES_iDR5_WG??_*.fits"):
    N = database.ingest_node_results(filename, extension=1)
    logger.info("Ingested {} results from {}".format(N, filename))


# Ingest additional photometric temperatures from Laura Magrini.
#for filename in glob("fits-templates/additional-tphots-magrini/*.fits"):
#    database.ingest_magrini_photometric_temperatures(filename)


# Ingest previous data release for comparison
database.ingest_recommended_results_from_previous_dr(
    "node-results/GES_iDR4_WG15_Recommended_Abundances_20042016.fits")


database.connection.commit()

logger.info("Ingestion complete.")


# Note that there is an issue with the CNAMEs of UVES benchmark spectra. There
# are two CNAME entries for alf_Cen_A, two for alf_Cet, and two for GJ880.

database.execute(
    """ UPDATE spectra
           SET cname = '14392972-6049560'
         WHERE ges_fld like 'alf_Cen_A%'""")

database.execute(
    """ UPDATE spectra
           SET cname = '03021676+0405219'
         WHERE ges_fld like 'alf_Cet%'""")

database.execute(
    """ UPDATE spectra
           SET cname = '22563384+1633085'
         WHERE ges_fld like 'GJ880%'""")

database.connection.commit()