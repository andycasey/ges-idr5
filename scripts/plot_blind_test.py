
"""
Make plots from the blind test.
"""

import yaml
import logging
import numpy as np
import os
from astropy.table import Table
from collections import OrderedDict


from code import GESDatabase, plot

logger = logging.getLogger("ges")

db_filename = "db.yaml"
with open(db_filename, "r") as fp:
    credentials = yaml.load(fp)
database = GESDatabase(**credentials)

# Load benchmark truths
benchmarks = Table.read("fits-templates/benchmarks/GES_iDR5_FGKM_Benchmarks_ARC_29092016.fits")

# Load giraffe blind spectra.
blind_giraffe = Table.read("fits-templates/blind-test/blind-test-giraffe.txt",
    format="ascii", names=("original", "blind"))
blind_giraffe = [(original, blind) for original, blind in \
    zip(blind_giraffe["original"], blind_giraffe["blind"])]


blind_uves = Table.read("fits-templates/blind-test/blind-test-uves.txt",
    format="ascii", names=("original", "blind"))
blind_uves = [(original, blind) for original, blind in \
    zip(blind_uves["original"], blind_uves["blind"])]


fig_uves = plot.blind.explain_differences(database, blind_uves,
    __comparison_filename="ges-dr5-uves-blind-test-comparison.fits",
    group_by=("wg", "name"))
fig_uves.savefig("figures/blind-test/ges-dr5-uves-blind-test-explaination.pdf", dpi=150)
fig_uves.savefig("figures/blind-test/ges-dr5-uves-blind-test-explaination.png", dpi=150)


fig_giraffe = plot.blind.explain_differences(database, blind_giraffe,
    __comparison_filename="ges-dr5-giraffe-blind-test-comparison.fits",
    group_by=("wg", "name"))
fig_giraffe.savefig("figures/blind-test/ges-dr5-giraffe-blind-test-explaination.pdf", dpi=150)
fig_giraffe.savefig("figures/blind-test/ges-dr5-giraffe-blind-test-explaination.png", dpi=150)


fig_giraffe, comparison_giraffe = plot.blind.differences_wrt_snr(
    database, benchmarks, blind_giraffe,
    group_by=("wg", "name", ), ylims=dict(teff=1000, logg=2.0, feh=2.0),
    __comparison_filename="ges-dr5-giraffe-blind-test-comparison.fits")
fig_giraffe.savefig("figures/blind-test/ges-dr5-giraffe-blind-test.pdf", dpi=150)
fig_giraffe.savefig("figures/blind-test/ges-dr5-giraffe-blind-test.png", dpi=150)


fig_uves, comparison_uves = plot.blind.differences_wrt_snr(
    database, benchmarks, blind_uves,
    group_by=("wg", "name", ), ylims=dict(teff=1000, logg=2.0, feh=2.0),
    __comparison_filename="ges-dr5-uves-blind-test-comparison.fits")
fig_uves.savefig("figures/blind-test/ges-dr5-uves-blind-test.pdf", dpi=150)
fig_uves.savefig("figures/blind-test/ges-dr5-uves-blind-test.png", dpi=150)


comparison_giraffe.write("ges-dr5-giraffe-blind-test-comparison.fits", overwrite=True)
comparison_uves.write("ges-dr5-uves-blind-test-comparison.fits", overwrite=True)




N = 20
xlims = dict(teff=500, logg=1, feh=1)
fig_uves_delta = plot.blind.differences_histogram(
    database, blind_uves, group_by=("wg", "name"),
    xlims=xlims, bins={ p: np.linspace(-l, +l, N) for p, l in xlims.items() },
    __comparison_filename="ges-dr5-uves-blind-test-comparison.fits")
fig_uves_delta.savefig("figures/blind-test/ges-dr5-uves-blind-test-delta.pdf", dpi=150)
fig_uves_delta.savefig("figures/blind-test/ges-dr5-uves-blind-test-delta.png", dpi=150)


fig_giraffe_delta = plot.blind.differences_histogram(
    database, blind_giraffe, group_by=("wg", "name"),
    xlims=xlims, bins={ p: np.linspace(-l, +l, N) for p, l in xlims.items() },
    __comparison_filename="ges-dr5-giraffe-blind-test-comparison.fits")
fig_giraffe_delta.savefig("figures/blind-test/ges-dr5-giraffe-blind-test-delta.pdf", dpi=150)
fig_giraffe_delta.savefig("figures/blind-test/ges-dr5-giraffe-blind-test-delta.png", dpi=150)
