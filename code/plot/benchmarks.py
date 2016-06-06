


import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.table import Table
from matplotlib.ticker import MaxNLocator
from matplotlib import gridspec

import utils

__all__ = ["node_benchmark_performance"]

benchmark_filename = "fits-templates/benchmarks/GES_iDR5_FGKMCoolWarm_Benchmarks_AcceptedParams_11052016.fits"


def node_benchmark_performance(database, wg, node_name, sort_by="TEFF"):
    """
    Show a box-and-whisker plot for the benchmark parameters reported by a given
    node.

    :param database:
        A database for transactions.

    :param wg:
        The working group.

    :param node_name:
        The node name.

    """

    width = 0.45
    colors=("#000000", "g", "#0874D4")


    node_id = database.retrieve_node_id(wg, node_name)

    benchmarks = Table.read(benchmark_filename)
    ok = np.isfinite(benchmarks["TEFF"] * benchmarks["LOGG"] * benchmarks["FEH"])
    benchmarks = benchmarks[ok]


    fig, axes = plt.subplots(3, figsize=(16.5, 7.5))
    parameters = ("teff", "logg", "mh")
    fits_parameters = ("TEFF", "LOGG", "FEH")
    ylabels = {
        "teff": r"$\Delta{}T_{\rm eff}$ $({\rm K})$",
        "logg": r"$\Delta\log{g}$",
        "mh": r"$\Delta[{\rm Fe}/{\rm H}]$",
        "feh": r"$\Delta[{\rm Fe}/{\rm H}]$"
    }
    benchmarks.sort(sort_by)
    for i, (ax, parameter) in enumerate(zip(axes, parameters)):

        ax.axhline(0, c="k", zorder=-1)

        fits_parameter = fits_parameters[i]
        
        x = benchmarks[fits_parameter]
        diff = []
        for benchmark in benchmarks:

            results = database.retrieve_table(
                """ SELECT  DISTINCT ON (r.filename) 
                            r.cname, ges_fld, feh, e_feh, {0}, e_{0}
                    FROM    results r, spectra s
                    WHERE   r.cname = s.cname
                    AND     r.node_id = %s
                    AND     TRIM(s.ges_fld) = %s""".format(parameter),
                (node_id, benchmark["GES_FLD"].strip()))

            if parameter == "feh":
                raise a


            if results is None:
                diff.append([])

            else:
                if parameter == "mh":
                    data = results[utils.mh_or_feh(results)]
                else:
                    data = results[parameter]

                data -= benchmark[fits_parameter]

                if np.any(np.isfinite(data)):
                    diff.append(data[np.isfinite(data)])
        
        # Show box-and-whisker plot.
        bp = ax.boxplot(diff, widths=width, patch_artist=True)
        ax.set_ylabel(ylabels.get(parameter))

        # Put numbers.
        if ax.is_last_row():
            y_loc = ax.get_ylim()[0] + 0.05 * np.ptp(ax.get_ylim())
            for j, d in enumerate(diff):
                ax.text(j + 1, y_loc, r"${}$".format(len(d)), color="k",
                    horizontalalignment="center")

            ax.set_xticklabels([each.strip() for each in benchmarks["GES_FLD"]])
            [l.set_rotation(90) for l in ax.get_xticklabels()]
        else:
            ax.set_xticklabels([])

        ax.yaxis.set_major_locator(MaxNLocator(5))
        ax.spines["left"]._linewidth = 0.5
        ax.spines["bottom"]._linewidth = 0.0
        ax.spines["top"]._linewidth = 0.0
        ax.spines["right"]._linewidth = 0.0

        opposite_ax = ax.twinx()
        opposite_ax.set_yticks([])

        plt.setp(bp["medians"], color="k", linewidth=2)
        plt.setp(bp["fliers"], color="k")
        plt.setp(bp["caps"], visible=False)
        plt.setp(bp["whiskers"], color="k", linestyle="solid", linewidth=0.5)
        plt.setp(bp["boxes"], color="k", alpha=0.5, linewidth=1)

    fig.tight_layout()
    return fig
