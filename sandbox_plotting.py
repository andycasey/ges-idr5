#!/usr/bin/python


import yaml
from code import GESDatabase

# Create a database object.
db_filename = "code/db.yaml"
with open(db_filename, "r") as fp:
    credentials = yaml.load(fp)
database = GESDatabase(**credentials)


from code import plot

fig = plot.hrd(database, 11, "Bologna")

fig_sun = plot.hrd(database, 11, "Bologna", mark=(5777, 4.4),
    where="CNAME = 'ssssssss-sssssss'")


raise a
fig = plot.stellar_parameter_histograms(database, 11, "Lumba")


fig = plot.stellar_parameter_error_histograms(database, 11, "EPINARBO")


from code import summary

foo = summary.stellar_parameter_summary(database)