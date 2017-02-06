

import logging
import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.table import Table
from matplotlib.ticker import MaxNLocator
from matplotlib import gridspec

import utils

__all__ = ["node_benchmark_performance", "wg_benchmark_performance"]

benchmark_filename = "fits-templates/benchmarks/GES_iDR5_FGKMCoolWarm_Benchmarks_AcceptedParams_01082016.fits"

logger = logging.getLogger("ges")

def node_benchmark_performance(database, wg, node_name, sort_by="TEFF",
    ylims=None):
    """
    Show a box-and-whisker plot for the benchmark parameters reported by a given
    node.

    :param database:
        A database for transactions.

    :param wg:
        The working group.

    :param node_name:
        The node name.

    :param ylims: [optional]
        A dictionary containing the absolute y-limit for each label (teff, logg,
        mh).
    """

    width = 0.45
    colors=("#000000", "g", "#0874D4")


    node_id = database.retrieve_node_id(wg, node_name)

    benchmarks = Table.read(benchmark_filename)
    ok = np.isfinite(benchmarks["TEFF"] * benchmarks["LOGG"] * benchmarks["FEH"])
    benchmarks = benchmarks[ok]


    parameters = ("teff", "logg", "mh")
    fits_parameters = ("TEFF", "LOGG", "FEH")
    
    fig, axes = plt.subplots(3, figsize=(16.5, 7.5))
    
    ylabels = {
        "teff": r"$\Delta{}T_{\rm eff}$ $({\rm K})$",
        "logg": r"$\Delta\log{g}$",
        "mh": r"$\Delta[{\rm Fe}/{\rm H}]$",
        "feh": r"$\Delta[{\rm Fe}/{\rm H}]$"
    }
    benchmarks.sort(sort_by)
    N = 0
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
                else:
                    diff.append([])
        
        # Show box-and-whisker plot.
        N += np.hstack(diff).size
        bp = ax.boxplot(diff, widths=width, patch_artist=True)
        ax.set_xlim(-0.5, len(benchmarks) + 0.5)
        
        ax.set_ylabel(ylabels.get(parameter))
        ylim = ylims.get(parameter, None)
        if ylim is not None:
            ax.set_ylim(-ylim, ylim)

        # Put numbers.
        if ax.is_last_row():
            y_loc = ax.get_ylim()[0] + 0.05 * np.ptp(ax.get_ylim())
            for j, d in enumerate(diff):
                ax.text(j + 1, y_loc, r"${}$".format(len(d)), color="k",
                    horizontalalignment="center")

            # Show how many are outside of each frame.
            if ylim:
                y_loc = ax.get_ylim()[1] - 0.05 * np.ptp(ax.get_ylim())
                for j, d in enumerate(diff):
                    N_bad = np.sum(np.abs(d) > ylim)
                    if N_bad == 0: continue
                    ax.text(j + 1, y_loc, r"${}$".format(N_bad), color="r",
                        horizontalalignment="center",
                        verticalalignment="top")

            ax.set_xticklabels(
                [each.strip().replace("_", "-") for each in benchmarks["GES_FLD"]])
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
    
    return fig if N > 0 else None



def wg_benchmark_performance(database, wg, truths, show_recommended=False,
    sort_by=None, ylims=None, node_sql_constraint=None, skip_missing=True,
    recommended_sql_constraint=None,
    show_num_estimates=None, xlabels=None, ylabels=None, **kwargs):
    """
    Show a box-and-whisker plot for the benchmark parameters reported by a given
    node.

    :param database:
        A database for transactions.

    :param wg:
        The working group.

    :param truths:
        An :astropy.table.Table: with the accepted "truth" values for the
        benchmark stars. A `GES_FLD` column is required for matching.

    :param show_recommended: [optional]
        Show the recommended values from the working group.

    :param sort_by: [optional]
        A column in the `truths` table to sort the benchmarks by.

    :param ylims: [optional]
        A dictionary containing the absolute y-limit for each label (teff, logg,
        mh).

    :param node_sql_constraint: [optional]
        An additional SQL constraint to apply to the node query.

    :param skip_missing: [optional]
        Skip benchmarks with zero node estimates.

    :param show_num_estimates: [optional]
        Show the number of node estimates for each benchmark star. If `None` is
        specified, then the number of estimates will be shown if no `ylims` are
        given. If `ylims` are given (and there are likely data points outside
        the range), then the number of estimates will only be shown if
        `show_num_estimates = True`.

    :param xlabels: [optional]
        A dictionary containing labels for each benchmark value, where the 
        `GES_FLD` entries are keys and the labels to use are values.

    :param ylabels: [optional]
        A dictionary containing labels for each parameter, where the parameters
        are keys and the labels to use are values.
    """

    width = kwargs.get("width", 0.45)
    colors = kwargs.get("colors", ("#000000", "#4daf4a", "#e41a1c"))

    
    show_num_estimates = True if ylims is None or show_num_estimates else False
    xlabels = xlabels or {}
    ylims = ylims or {}

    parameters = kwargs.get("parameters", ("TEFF", "LOGG", "FEH"))
    fig, axes = plt.subplots(
        len(parameters), figsize=kwargs.get("figsize", (12, 2.5 * len(parameters))))
    
    default_ylabels = {
        "TEFF": r"$\Delta{}T_{\rm eff}$ $({\rm K})$",
        "LOGG": r"$\Delta\log{g}$",
        "FEH": r"$\Delta[{\rm Fe}/{\rm H}]$"
    }
    default_ylabels.update(ylabels or {})
    ylabels = default_ylabels

    truths = truths.copy()
    truths.sort(sort_by or parameters[0])
    N = 0
    diffs = {}
    recommended_diffs = {}

    #for i, (ax, parameter) in enumerate(zip(axes, parameters)):
    for i, parameter in enumerate(parameters):

        x = truths[parameter]

        diffs[parameter] = []
        recommended_diffs[parameter] = []

        for benchmark in truths:

            kwds = dict(
                parameter=parameter, wg=wg, ges_fld=benchmark["GES_FLD"].strip(),
                node_sql_constraint_str=" AND {}".format(node_sql_constraint) \
                    if node_sql_constraint is not None else "",
                recommended_sql_constraint_str = "" if recommended_sql_constraint is None \
                                    else " AND {} ".format(recommended_sql_constraint),
                recommended_table=kwargs.get(
                    "recommended_table", "wg_recommended_results"))

            node_estimates = database.retrieve_table(
                """ SELECT DISTINCT ON (results.id)
                            {parameter}, e_{parameter}
                      FROM  results, nodes, spectra
                     WHERE  nodes.wg = '{wg}'
                       AND  results.node_id = nodes.id
                       AND  results.cname = spectra.cname
                       AND  TRIM(spectra.ges_fld) = '{ges_fld}'
                       AND  {parameter} <> 'NaN'
                       {node_sql_constraint_str};""".format(**kwds))

            if not np.isfinite(benchmark[parameter]):
                logger.warn("No finite {} `truth` value for {}!".format(
                    parameter, benchmark["GES_FLD"].strip()))

            if node_estimates is None:
                # No node estimates for this object
                diffs[parameter].append([])

            else:
                diffs[parameter].append(
                    node_estimates[parameter.lower()] - benchmark[parameter])
                
            if show_recommended:
                wg_recommended = database.retrieve_table(
                    """ SELECT  {parameter}, e_{parameter}, 
                                e_pos_{parameter}, e_neg_{parameter}
                          FROM  {recommended_table} as wgr, spectra
                         WHERE  wgr.wg = {wg}
                           AND  wgr.cname = spectra.cname
                           AND  TRIM(spectra.ges_fld) = '{ges_fld}'
                           {recommended_sql_constraint_str};
                    """.format(**kwds))

                if wg_recommended is None:
                    recommended_diffs[parameter].append([])

                else:
                    recommended_diffs[parameter].append([
                        wg_recommended[parameter.lower()][0] - benchmark[parameter],
                        wg_recommended["e_{}".format(parameter.lower())][0]
                    ])

    # Show benchmarks with missing entries?
    if skip_missing:
        keep = np.sum(np.array([map(len, v) for v in diffs.values()]), axis=0) > 0

        truths = truths[keep]
        for parameter in parameters:
            diffs[parameter] = np.array(diffs[parameter])[keep]
            recommended_diffs[parameter] = np.array(recommended_diffs[parameter])[keep]

    axes = np.array([axes]).flatten()

    for i, (ax, parameter) in enumerate(zip(axes, parameters)):

        ax.axhline(0, c="k", zorder=-1)

        diff = diffs[parameter]
        recommended_diff = recommended_diffs[parameter]

        # Show box-and-whisker plot.
        N += np.hstack(diff).size
        bp = ax.boxplot(diff, widths=width, patch_artist=True)
        
        # Show the recommended values.
        if show_recommended:
            for j, result in enumerate(recommended_diff):
                if not result: continue
                mu_diff, sigma = result
                ax.fill_between(
                    [j + 1 - width/2, j + 1 + width/2],
                    [mu_diff - sigma, mu_diff - sigma],
                    [mu_diff + sigma, mu_diff + sigma],
                    linewidth=0, alpha=0.5, facecolor=colors[1])
                ax.plot(
                    [j + 1 - width/2, j + 1 + width/2],
                    [mu_diff, mu_diff],
                    lw=2, color=colors[1])

        ax.set_ylabel(ylabels.get(parameter))
        ax.set_xlim(0.5, len(truths) + 0.5)
        ylim = ylims.get(parameter, None)
        if ylim is not None:
            try:
                ylim[0]
            except:
                ax.set_ylim(-ylim, ylim)
            else:
                ax.set_ylim(*ylim)

            # Show how many are outside of each frame.
            y_loc = ax.get_ylim()[1] - 0.075 * np.ptp(ax.get_ylim())
            for j, d in enumerate(diff):
                N_bad = np.sum(d > ax.get_ylim()[1])
                if N_bad == 0: continue
                ax.plot([j + 1], [y_loc], marker="^", 
                    c=colors[2], linewidth=0)
                ax.text(j + 1.4, y_loc, r"${}$".format(N_bad), 
                    color=colors[2], fontsize=9, horizontalalignment="center",
                    verticalalignment="center")
            
            y_loc = ax.get_ylim()[0] + 0.075 * np.ptp(ax.get_ylim())
            for j, d in enumerate(diff):
                N_bad = np.sum(d < ax.get_ylim()[0])
                if N_bad == 0: continue
                ax.plot([j + 1], [y_loc], marker="v", 
                    c=colors[2], linewidth=0)
                ax.text(j + 1.4, y_loc, r"${}$".format(N_bad), 
                    color=colors[2], fontsize=9, horizontalalignment="center",
                    verticalalignment="center")

        if show_num_estimates and ax.is_last_row():
            y_loc = ax.get_ylim()[0] + 0.05 * np.ptp(ax.get_ylim())
            for j, d in enumerate(diff):
                ax.text(j + 1, y_loc, r"${}$".format(len(d)), color="k",
                    horizontalalignment="center")

        if ax.is_last_row():
            #ax.set_xticks(range(len(truths)))
            ax.set_xticklabels([
                xlabels.get(each, each.strip().replace("_", "-")) \
                    for each in truths["GES_FLD"]])
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

        plt.setp(bp["medians"], color=colors[0], linewidth=2)
        plt.setp(bp["fliers"], color=colors[0])
        plt.setp(bp["caps"], visible=False)
        plt.setp(bp["boxes"], color=colors[0], alpha=0.5, linewidth=1)
        plt.setp(
            bp["whiskers"], color=colors[0], linestyle="solid", linewidth=0.5)
        
        ax.xaxis.set_tick_params(width=0)
        
    fig.tight_layout()
    
    return fig if N > 0 else None
