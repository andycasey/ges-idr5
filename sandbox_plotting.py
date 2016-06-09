#!/usr/bin/python


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
db_filename = "code/db.yaml"
with open(db_filename, "r") as fp:
    credentials = yaml.load(fp)
database = GESDatabase(**credentials)




wgs = (10, 11, 12, 13)
parameters = ("teff", "logg", "mh", "xi")
isochrones = glob("isochrones/*.dat")


# Create node-level folders.
nodes = database.retrieve_table("SELECT * FROM nodes")
for node in nodes:
    folder = "figures/node-level/{}/{}".format(node["wg"], node["name"])
    if not os.path.exists(folder):
        os.mkdir(folder)

# Plot against clusters.
cluster_velocities = Table.read("fits-templates/oc_gc_radial_vel.dat", format="ascii")
for node in nodes:

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

            fig = plot.cluster(database, node["wg"], node["name"], cluster,
                isochrone_filename=isochrone, vel_range=vel_range)

        except:
            logger.exception(
                "Could not create cluster {} figure for {}/{}:".format(
                    cluster, node["wg"], node["name"]))

        else:
            if fig is not None:
                basename = "figures/node-level/{}/{}/cluster-{}".format(
                    node["wg"], node["name"], cluster)

                fig.savefig("{}.pdf".format(basename))
                fig.savefig("{}.png".format(basename))




        try:
            fig = plot.cluster(database, node["wg"], node["name"], cluster,
                isochrone_filename=isochrone, vel_range=vel_range,
                limit_to_isochrone_range=True)

        except:
            logger.exception(
                "Could not create cluster {} figure for {}/{}:".format(
                    cluster, node["wg"], node["name"]))

        else:
            if fig is not None:
                basename = "figures/node-level/{}/{}/cluster-{}-limited".format(
                    node["wg"], node["name"], cluster)

                fig.savefig("{}.pdf".format(basename))
                fig.savefig("{}.png".format(basename))


        # Plot against parameters for stars in this cluster.
        try:
            fig = plot.param_vs_param(database, node["wg"], node["name"],
                cluster, "teff", vel_range=vel_range)

        except:
            logger.exception(
                "Could not create param_vs_param plot for cluster {} for {}/{}:"\
                .format(cluster, node["wg"], node["name"]))

        else:
            basename = "figures/node-level/{}/{}/cluster-{}-vs-teff".format(
                node["wg"], node["name"], cluster)

            fig.savefig("{}.pdf".format(basename))
            fig.savefig("{}.png".format(basename))


        try:
            fig = plot.param_vs_param(database, node["wg"], node["name"],
                cluster, "logg", vel_range=vel_range)

        except:
            logger.exception(
                "Could not create param_vs_param plot for cluster {} for {}/{}:"\
                .format(cluster, node["wg"], node["name"]))

        else:
            basename = "figures/node-level/{}/{}/cluster-{}-vs-logg".format(
                node["wg"], node["name"], cluster)

            fig.savefig("{}.pdf".format(basename))
            fig.savefig("{}.png".format(basename))


        plt.close("all")


raise a


# H-R diagrams by setup.
for node in nodes:

    print("Plotting HR diagrams for node {}".format(node))

    try:
        fig = plot.hrd_by_setup(database, node["wg"], node["name"])

    except:
        raise 
        logger.exception(
            "Could not create hrd_by_setup figure for {}/{}:"\
            .format(node["wg"], node["name"]))
        continue

    else:
        if fig is None: continue
        basename = "figures/node-level/{wg}/{name}/hrd-by-setup".format(
            wg=node["wg"], name=node["name"])

        fig.savefig("{}.pdf".format(basename))
        fig.savefig("{}.png".format(basename))

plt.close("all")





# Plot the performance on the benchmarks.
for node in nodes:

    print("Plotting benchmark performance for {}".format(node))
    try:
        fig = plot.node_benchmark_performance(database, node["wg"], node["name"])

    except:
        logger.exception(
            "Could not create node_benchmark_performance figure for {}/{}:"\
            .format(node["wg"], node["name"]))
        continue

    else:
        basename = "figures/node-level/{wg}/{name}/benchmarks".format(
            wg=node["wg"], name=node["name"])

        fig.savefig("{}.pdf".format(basename))
        fig.savefig("{}.png".format(basename))

plt.close("all")




# Node-to-node comparisons within a WG
for wg in wgs:
    for parameter in parameters:

        print("Plotting node-to-node comparison of {} for wg {}".format(
            parameter, wg))
        try:
            fig = plot.compare_nodes_within_wg(database, wg, parameter)

        except:
            logger.exception(
                "Could not create compare_nodes_within_wg figure for {}/{}:"\
                .format(wg, parameter))
            continue

        else:
            fig.savefig("figures/wg-level/{}/{}.pdf".format(wg, parameter))
            fig.savefig("figures/wg-level/{}/{}.png".format(wg, parameter))

plt.close("all")



# Compare temperatures to photometric temperatures.
for node in nodes:

    print("Plotting photometric temperature comparisons for {}".format(node))

    try:
        fig = plot.compare_to_photometric_teff(database, node["wg"], node["name"])

    except:
        logger.exception(
            "Could not create compare_to_photometric_teff figure for {}/{}"\
            .format(node["wg"], node["name"]))
        continue

    else:
        basename = "figures/node-level/{}/{}/photometric-teff".format(
            node["wg"], node["name"])

        fig.savefig("{}.pdf".format(basename))
        fig.savefig("{}.png".format(basename))


plt.close("all")



# Compare to previous data release.
for node in nodes:
    for parameter in parameters:

        print("Comparing node {} parameters of {} to previous DR".format(node,
            parameter))
        try:
            fig = plot.compare_to_previous_dr(
                database, node["wg"], node["name"], parameter)

        except:
            logger.exception(
                "Could not create compare_to_previous_dr figure for {}/{}/{}:"\
                .format(node["wg"], node["name"], parameter))
            continue

        else:
            basename = "figures/node-level/{}/{}/idr4-compare-{}".format(
                node["wg"], node["name"], parameter)

            fig.savefig("{}.pdf".format(basename))
            fig.savefig("{}.png".format(basename))


plt.close("all")




# Total HRD.
for node in nodes:

    print("Plotting total HRD for {}".format(node))
    try:
        fig = plot.hrd(database, node["wg"], node["name"])

    except:
        logger.exception(
            "Could not create hrd figure for {}/{}:".format(
                node["wg"], node["name"]))
        continue

    else:
        basename = "figures/node-level/{wg}/{name}/hrd".format(
            wg=node["wg"], name=node["name"])

        fig.savefig("{}.pdf".format(basename))
        fig.savefig("{}.png".format(basename))

plt.close("all")


# Plot results for the Sun.
for node in nodes:

    print("PLotting solar results for {}".format(node))

    try:
        fig = plot.hrd(database, node["wg"], node["name"], mark=(5777, 4.4),
            where="CNAME = 'ssssssss-sssssss'")

    except:
        logger.exception(
            "Could not create hrd (SUN) figure for {}/{}:".format(
                node["wg"], node["name"]))
        continue

    else:
        basename = "figures/node-level/{wg}/{name}/solar".format(
            wg=node["wg"], name=node["name"])

        fig.savefig("{}.pdf".format(basename))
        fig.savefig("{}.png".format(basename))


plt.close("all")



# Plot distributions of stellar parameters.
for node in nodes:

    print("Plotting stellar parameter distributions for {}".format(node))

    try:
        fig = plot.stellar_parameter_histograms(database, node["wg"], node["name"])

    except:
        logger.exception(
            "Could not create stellar_parameter_histograms for {}/{}:".format(
                node["wg"], node["name"]))
        continue

    else:
        basename = "figures/node-level/{}/{}/parameter-histogram".format(
            node["wg"], node["name"])

        fig.savefig("{}.pdf".format(basename))
        fig.savefig("{}.png".format(basename))


plt.close("all")

# Plot distributions of errors in stellar parameters.
for node in nodes:

    print("Plotting stellar parameter error distributions for {}".format(node))
    
    try:
        fig = plot.stellar_parameter_error_histograms(database, node["wg"], node["name"])

    except:
        logger.exception(
            "Could not create stellar_parameter_error_histograms for {}/{}:".format(
                node["wg"], node["name"]))
        continue

    else:
        basename = "figures/node-level/{}/{}/parameter-error-histogram".format(
            node["wg"], node["name"])

        fig.savefig("{}.pdf".format(basename))
        fig.savefig("{}.png".format(basename))

plt.close("all")





# Create summary tables.
parameter_summary = summary.stellar_parameter_summary(database)
parameter_summary.write("figures/wg-level/parameter-summary.txt",
    format="ascii")

parameter_ranges = summary.stellar_parameter_range(database)
parameter_ranges.write("figures/wg-level/parameter-range-summary.txt",
    format="ascii")

for node in nodes:

    tech = summary.tech_flags(database, node["wg"], node["name"])
    if tech is not None:
        tech.write(
            "figures/tech-summary-{}-{}.txt".format(node["wg"], node["name"]),
            format="ascii")



