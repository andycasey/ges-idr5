#!/usr/bin/python

"""
Create database structure for iDR5 quality control (QC) phase.
"""

import logging
import numpy as np
import psycopg2 as pg
import yaml

# For fake data generation
import os
from astropy.io import fits

from code import GESDatabase


db_filename = "code/db.yaml"
nodes_filename = "nodes.yaml"
schema_filename = "code/schema.sql"
masterlist_filename \
    = "fits-templates/masterlist/GES_iDR5_spectra_masterlist_15052016.fits"

raise AreYouSureError

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

# Create some fake data.
extension = -1
filename = "fits-templates/node-templates/GES_iDR5_WG11_NodeTemplate_QC_16052016.fits"
node_names = ("Bologna", "Lumba", "ULB", "MyGIsFOS", "CAUP", "Vilnius", "EPINARBO",)

fake_modes = {
    "Bologna": {
        "teff": (5000, 400),
        "e_teff": (100, 20),
        "logg": (3.5, 0.5),
        "e_logg": (0.3, 0.05),
        "mh": (0, 0.5),
        "e_mh": (0.1, 0.02),
        "xi": (1, 0.2),
    },
    "Lumba": {
        "teff": (5000, 400),
        "e_teff": (100, 10),
        "logg": (3.5, 0.5),
        "mh": (0, 0.5),
        "xi": (1, 0.2)
    },
    "ULB": {
        "teff": (5000, 400),
        "logg": (3.5, 0.5),
        "mh": (0, 0.5),
        "xi": (1, 0.2)
    },
    "MyGIsFOS": {
        "teff": (5000, 400),
        "logg": (3.5, 0.5),
        "mh": (0, 0.5),
        "xi": (1, 0.2)
    },
    "CAUP": {
        "teff": (5000, 400),
        "logg": (3.5, 0.5),
        "mh": (0, 0.5),
        "xi": (1, 0.2)
    },
    "Vilnius": {
        "teff": (5000, 400),
        "logg": (3.5, 0.5),
        "mh": (0, 0.5),
        "xi": (1, 0.2)
    },
    "EPINARBO": {
        "teff": (5000, 400),
        "logg": (3.5, 0.5),
        "mh": (0, 0.5),
        "xi": (1, 0.2)
    },
}

for node_name in node_names:
    fake_filename = "GES_iDR5_WG11_{}_FAKE.fits".format(node_name)

    os.system("cp {} {}".format(filename, fake_filename))    
    image = fits.open(fake_filename)

    N = len(image[extension].data)
    for key, (mu, sigma) in fake_modes[node_name].items():
        image[extension].data[key] = np.random.normal(mu, sigma, size=N)

    image.writeto(fake_filename, clobber=True)
    logger.info("Created fake node result file {}".format(fake_filename))

# Ingest faux results from the nodes.
from glob import glob
for filename in glob("GES_iDR5_WG??_*_FAKE.fits"):
    N = database.ingest_node_results(filename)
    logger.info("Ingested {} fake results from {}".format(N, filename))

# Make plots.


