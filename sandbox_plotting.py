#!/usr/bin/python


import yaml
from code import (GESDatabase, plot, summary)

# Create a database object.
db_filename = "code/db.yaml"
with open(db_filename, "r") as fp:
    credentials = yaml.load(fp)
database = GESDatabase(**credentials)



# Node-to-node comparisons within a WG.

#fig = plot.compare_nodes_within_wg(database, 11, "logg")



#fig = plot.compare_to_photometric_teff(database, 11, "EPINARBO")


#fig = plot.compare_to_previous_dr(database, 11, "Vilnius", "teff")

#fig = plot.hrd_by_setup(database, 11, "Vilnius")


#fig = plot.node_benchmark_performance(database, 11, "Vilnius")


#fig = plot.hrd(database, 11, "EPINARBO")


#fig_sun = plot.hrd(database, 11, "CAUP", mark=(5777, 4.4),
#    where="CNAME = 'ssssssss-sssssss'")


#bar = summary.stellar_parameter_range(database)
#foo = plot.param_vs_param(database, 11, "EPINARBO", "NGC6705", "xi")

#foo = plot.cluster(database, 11, "EPINARBO", "NGC6705",
#    isochrone_filename="isochrones/NGC6705_Parsec_0.3Gyr_Z0.018.dat")

fig = plot.stellar_parameter_histograms(database, 11, "Lumba")


fig = plot.stellar_parameter_error_histograms(database, 11, "EPINARBO")



foo = summary.stellar_parameter_summary(database)

raise a
