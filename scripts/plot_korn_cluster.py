
import yaml
import logging
import os
import matplotlib.pyplot as plt
from glob import glob

from code import (GESDatabase, plot, summary)
from astropy.table import Table

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



debug = False
wgs = (10, 11, 12, 13)
parameters = ("teff", "logg", "mh", "xi")
isochrones = glob("isochrones/*.dat")

isochrones = ["isochrones/M67_Parsec_4.0Gyr_Z0.017.dat"]

# Plot against clusters.
cluster_velocities = Table.read("fits-templates/oc_gc_radial_vel.dat", format="ascii")

nodes = database.retrieve_table("SELECT id, wg, TRIM(name) as name FROM nodes WHERE wg = 11")

fig, axes = plt.subplots(2, 3)
axes = np.array(axes).flatten()
index = 0

node_names = []
for node in nodes:

    try:
        ax = axes[index]
    except IndexError:
        print("Whoa we ran outa axes dude!")
        break


    for isochrone in isochrones:
        cluster = isochrone.split("/")[-1].split("_")[0]
        if cluster == "gamma2-Vel":
            cluster = "gamma2_Vel"

        print("Plotting cluster {} for node {}".format(cluster, node))

        match = cluster_velocities["id"] == cluster
        if not any(match):
            vel_range = None

        else:
            vrad = cluster_velocities["rv"][match]
            e_vrad = cluster_velocities["erv"][match]
            vel_range = (vrad - 2*e_vrad, vrad + 2*e_vrad)
            vel_range = (vel_range[0][0], vel_range[1][0])

        try:

            r = plot.cluster(database, node["wg"], node["name"], cluster,
                isochrone_filename=isochrone, vel_range=vel_range, ax=ax,
                no_tech_flags=True, limit_to_isochrone_range=True,
                show_legend=False, vmin=-0.3, vmax=0.15)

        except:
            if debug: raise
            logger.exception(
                "Could not create cluster {} figure for {}/{}:".format(
                    cluster, node["wg"], node["name"]))

        else:
            if r is not None:
                node_names.append(node["name"])
                index += 1
            continue

xticks = [3000, 4000, 5000, 6000, 7000]
for ax, name in zip(axes, node_names):

    ax.text(0.05, 0.90, name, color="k", transform=ax.transAxes)
    if not ax.is_first_col():
        ax.set_ylabel("")
        ax.set_yticklabels([])

    ax.set_xticks(xticks)

    if ax.is_last_row():
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks)
        ax.set_xlabel(r"$T_{\rm eff}$ $[{\rm K}]$")

fig.tight_layout()
cbar = plt.colorbar(r, ax=list(axes))
cbar.set_label(r"[Fe/H]")

raise a