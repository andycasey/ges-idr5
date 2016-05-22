#!/usr/bin/python


import yaml
from code import GESDatabase

# Create a database object.
db_filename = "code/db.yaml"
with open(db_filename, "r") as fp:
    credentials = yaml.load(fp)
database = GESDatabase(**credentials)


from code import plot

#fig = plot.compare_nodes_within_wg(database, 11, "teff")
#fig = plot.compare_to_photometric_teff(database, 11, "Bologna")

#fig = plot.compare_to_previous_dr(database, 11, "Bologna", "teff")

fig = plot.hrd_by_setup(database, 11, "Bologna")


raise a

fig = plot.hrd(database, 11, "Bologna")

fig_sun = plot.hrd(database, 11, "Bologna", mark=(5777, 4.4),
    where="CNAME = 'ssssssss-sssssss'")

from code import summary

bar = summary.stellar_parameter_range(database)

foo = plot.param_vs_param(database, 11, "Bologna", "NGC6705", "xi")
raise a

foo = plot.cluster(database, 11, "Bologna", "NGC6705",
    isochrone_filename="isochrones/NGC6705_Parsec_0.3Gyr_Z0.018.dat")

raise a
fig = plot.stellar_parameter_histograms(database, 11, "Lumba")


fig = plot.stellar_parameter_error_histograms(database, 11, "EPINARBO")



foo = summary.stellar_parameter_summary(database)