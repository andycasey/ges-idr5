
import logging
import numpy as np
import psycopg2 as pg
import yaml
from glob import glob
from astropy.table import Table

from code import GESDatabase

logger = logging.getLogger("ges")

# Get credentials.
db_filename = "db.yaml"
with open(db_filename, "r") as fp:
    credentials = yaml.load(fp)

# Create a database object.
database = GESDatabase(**credentials)


"""
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

"""

# Ingest WG14 Halpha results
node_id = database.retrieve_node_id(14, "Halpha")
for filename in glob("node-results/WG14/Halpha_data_for_template_*_by_CNAME.csv"):

    data = Table.read(filename, format="csv")
    N = len(data)

    for i, row in enumerate(data):

        tech = "|".join(row["TECH"].strip().split("|"))
        peculi = "|".join(row["PECULI"].strip().split("|"))

        row_data = dict(node_id=node_id)
        row_data.update(
            cname=row["CNAME"],
            filename=row["FILENAME"],
            vel=row["VEL"],
            e_vel=row["E_VEL"],
            vrot=row["VROT"],
            e_vrot=row["E_VROT"],
            tech=tech, peculi=peculi)

        logger.info("Inserting row {}/{}: {}".format(i, N, row_data))

        database.execute(
            "INSERT INTO results ({}) VALUES ({})".format(
                ", ".join(row_data.keys()),
                ", ".join(["%({})s".format(column) for column in row_data.keys()])),
            row_data)


database.connection.commit()
