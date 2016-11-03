
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


database.connection.commit()
