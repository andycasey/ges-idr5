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


# Ingest WG14 BIN results
data = Table.read("node-results/WG14/iDR5_SB234.txt", names=("cname", "peculi"), format="ascii")
node_id = database.retrieve_node_id(14, "BIN")

N = len(data)
for i, row in enumerate(data):
    row_data = dict(node_id=node_id, filename="ALL", setup="ALL")
    row_data.update(cname=row["cname"], peculi=row["peculi"])

    logger.info("Inserting row {}/{}: {}".format(i, N, row_data))
   
    database.execute(
        "INSERT INTO results ({}) VALUES ({})".format(
            ", ".join(row_data.keys()),
            ", ".join(["%({})s".format(column) for column in row_data.keys()])),
        row_data)


# Ingest WG14 Halpha results
node_id = database.retrieve_node_id(14, "Halpha")
for filename in glob("node-results/WG14/Halpha_data_for_template_*_by_CNAME.csv"):

    data = Table.read(filename, format="csv")
    N = len(data)

    _DEFAULTS = {
        "setup": "ALL",
        "filename": "ALL",
        "vel": np.nan,
        "e_vel": np.nan,
        "vrot": np.nan,
        "e_vrot": np.nan
    }

    for i, row in enumerate(data):
        
        tech = "{}".format(row["TECH"]).strip()
        peculi = "{}".format(row["PECULI"]).strip()
        tech = tech if tech != "--" else ""
        peculi = peculi if peculi != "--" else ""

        tech = "|".join(list(set(tech.split("|"))))
        peculi = "|".join(list(set(peculi.split("|"))))

        row_data = dict(node_id=node_id)
        row_data.update(
            cname=row["CNAME"],
            setup=row["SETUP"],
            filename=row["FILENAME"],
            vel=row["VEL"],
            e_vel=row["E_VEL"],
            vrot=row["VROT"],
            e_vrot=row["E_VROT"],
            tech=tech, peculi=peculi)

        # Check for masked data
        for key in row_data.keys():
            if isinstance(row_data[key], np.ma.core.MaskedArray):
                row_data[key] = _DEFAULTS[key]

        logger.info("Inserting row {}/{}: {}".format(i, N, row_data))

        database.execute(
            "INSERT INTO results ({}) VALUES ({})".format(
                ", ".join(row_data.keys()),
                ", ".join(["%({})s".format(column) for column in row_data.keys()])),
            row_data)


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


# Fix the MH/FEH issue:
database.execute(
    """ UPDATE results
           SET feh = mh
         WHERE feh = 'NaN'
           AND mh <> 'NaN';""")

database.connection.commit()