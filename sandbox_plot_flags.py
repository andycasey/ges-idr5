#!/usr/bin/python


import yaml
import logging
import os
import matplotlib.pyplot as plt
from glob import glob

from code import (GESDatabase, plot, summary)

# Initialize logging.
logger = logging.getLogger("ges.idr5.qc")
logger.setLevel(logging.DEBUG)

handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter(
    "%(asctime)s [%(levelname)-8s] %(message)s"))
logger.addHandler(handler)


# Create a database object.
db_filename = "db.yaml"
with open(db_filename, "r") as fp:
    credentials = yaml.load(fp)
database = GESDatabase(**credentials)


for wg in (10, 11, 12, 13):

    fig = plot.flags.heatmap(database, wg,
        show_multiple_flags_per_node=True, group_by="node",
        use_cuthill_mckee=True)
    fig.savefig(
        "figures/wg{}/flags-heatmap-node.pdf".format(wg), dpi=300)

    fig = plot.flags.heatmap(database, wg,
        show_multiple_flags_per_node=True, group_by="node",
        use_cuthill_mckee=False)
    fig.savefig(
        "figures/wg{}/flags-heatmap-node-ordered.pdf".format(wg), dpi=300)


    print("Done {}".format(wg))