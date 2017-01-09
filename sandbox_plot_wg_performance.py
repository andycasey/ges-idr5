
import yaml
from glob import glob

from code import (GESDatabase, plot, summary)
from astropy.table import Table

# Create a database object.
db_filename = "db.yaml"
with open(db_filename, "r") as fp:
    credentials = yaml.load(fp)
database = GESDatabase(**credentials)

# Clean up bits and pieces...


savefig_kwds = dict(dpi=300, bbox_inches="tight")

benchmarks = Table.read("fits-templates/benchmarks/GES_iDR5_FGKM_Benchmarks_ARC_29092016.fits")
benchmarks = benchmarks[benchmarks["TEFF"] < 8000]

isochrones = glob("isochrones/*.dat")

debug = True


cluster_velocities = Table.read("fits-templates/oc_gc_radial_vel.dat", format="ascii")

for wg in (10, 11, 20, ):

    for isochrone in isochrones:

        cluster = isochrone.split("/")[-1].split("_")[0]
        if cluster == "gamma2-Vel":
            cluster = "gamma2_Vel"

        match = cluster_velocities["id"] == cluster
        if not any(match):
            vel_range = None

        else:
            vrad = cluster_velocities["rv"][match]
            e_vrad = cluster_velocities["erv"][match]
            vel_range = (vrad - 2*e_vrad, vrad + 2*e_vrad)
            vel_range = (vel_range[0][0], vel_range[1][0])

        try:

            fig = plot.cluster(database, cluster, wg,
                isochrone_filename=isochrone, vel_range=vel_range)

        except:
            if debug: raise
            logger.exception(
                "Could not create cluster {} figure for {}:".format(cluster, wg))

        else:
            if fig is not None:
                basename = "figures/wg{wg}/wg{wg}-cluster-{cluster}-poly".format(
                    wg=wg, cluster=cluster)

                fig.savefig("{}.pdf".format(basename), **savefig_kwds)
                fig.savefig("{}.png".format(basename), **savefig_kwds)



        try:
            fig = plot.cluster(database, cluster, wg,
                isochrone_filename=isochrone, vel_range=vel_range,
                limit_to_isochrone_range=True)

        except:
            if debug: raise
            logger.exception(
                "Could not create cluster {} figure for {}:".format(cluster, wg))

        else:
            if fig is not None:
                basename = "figures/wg{wg}/wg{wg}-cluster-{cluster}-limited-poly"\
                            .format(wg=wg, cluster=cluster)

                fig.savefig("{}.pdf".format(basename), **savefig_kwds)
                fig.savefig("{}.png".format(basename), **savefig_kwds)


    # Plot benchmarks first.
    fig = plot.wg_benchmark_performance(database, wg, benchmarks, 
        show_recommended=True, ylims=dict(TEFF=1000, LOGG=1, FEH=1),

        recommended_table="wg_recommended_results")
    fig.savefig("figures/wg{wg}/wg{wg}-benchmarks-zoom-poly.pdf".format(wg=wg), **savefig_kwds)
    fig.savefig("figures/wg{wg}/wg{wg}-benchmarks-zoom-poly.png".format(wg=wg), **savefig_kwds)


    fig = plot.wg_benchmark_performance(
        database, wg, benchmarks, show_recommended=True,
         recommended_table="wg_recommended_results")
    fig.savefig("figures/wg{wg}/wg{wg}-benchmarks-poly.pdf".format(wg=wg), **savefig_kwds)
    fig.savefig("figures/wg{wg}/wg{wg}-benchmarks-poly.png".format(wg=wg), **savefig_kwds)


raise a

