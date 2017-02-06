#!/usr/bin/python


import yaml
from code import (GESDatabase, plot)

# Create a database object.
db_filename = "db.yaml"
with open(db_filename, "r") as fp:
    credentials = yaml.load(fp)
database = GESDatabase(**credentials)


for wg in (10, 11, 12, 13):

    fig, meta = plot.flags.heatmap(database, wg,
        show_multiple_flags_per_node=True, group_by="node",
        use_cuthill_mckee=True)
    fig.savefig(
        "figures/wg{}/flags-heatmap-node.pdf".format(wg), dpi=300)

    fig, meta = plot.flags.heatmap(database, wg,
        show_multiple_flags_per_node=True, group_by="node",
        use_cuthill_mckee=False)
    fig.savefig(
        "figures/wg{}/flags-heatmap-node-ordered.pdf".format(wg), dpi=300)


    print("Done {}".format(wg))