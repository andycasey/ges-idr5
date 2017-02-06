
"""
Produce figures related to the blind tests.
"""

import logging
import numpy as np
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import yaml
from matplotlib.ticker import MaxNLocator
from matplotlib import gridspec
from astropy.table import Table
from collections import Counter, OrderedDict

import utils

__all__ = []

logger = logging.getLogger("ges.idr5")

def _retrieve_blind_test_result(database, original_filename, blind_filename,
    additional_columns=None, sql_constraint=None):
    """
    Retrieve a blind test result.

    :param original_filename:
        The original filename.

    :param blind_filename:
        The blind filename that was used in the test.

    :param sql_constraint: [optional]
        Any additional SQL constraint to place on the query.

    :returns:
        A two-length tuple containing a table of results based on the original
        filename, and a table of results based on the blind filename. If no
        results were available for either query, then `None` will be returned
        instead of a table.
    """

    sql_constraint \
        = "" if sql_constraint is None else "AND {}".format(sql_constraint)

    # Get the original and blind CNAMEs first.
    #original_cname, blind_cname = [database.retrieve(
    #    """ SELECT cname
    #        FROM spectra 
    #        WHERE filename LIKE '%{}%';""".format(os.path.basename(path)))[0][0] \
    #    for path in (original_filename, blind_filename)]

    # Get results for each.
    # Add columns that we need for reasons.
    columns = [
        "node_id", "nodes.name", "nodes.wg", 
        "teff", "e_teff", "logg", "e_logg", "feh", "e_feh",
        "trim(results.tech) as tech"]
    if additional_columns is not None:
        columns = list(set(columns + list(additional_columns)))

    original_results = database.retrieve_table(
        """ SELECT DISTINCT ON (results.id) {columns}
            FROM results, nodes, spectra
            WHERE results.cname = spectra.cname
              AND results.filename LIKE '%{basename}%'
              AND results.node_id = nodes.id
            {sql_constraint};""".format(
                basename=os.path.basename(original_filename),
                columns=", ".join(columns),
                sql_constraint=sql_constraint))

    blind_results = database.retrieve_table(
        """ SELECT DISTINCT ON (results.id) {columns}
            FROM results, nodes, spectra
            WHERE results.cname = spectra.cname
              AND results.filename LIKE '%{basename}%'
              AND results.node_id = nodes.id
            {sql_constraint};""".format(
                basename=os.path.basename(blind_filename),
                columns=", ".join(columns),
                sql_constraint=sql_constraint))

    return (original_results, blind_results)


def explain_differences(database, blind_test_filenames, group_by=None, 
    sql_constraint=None, **kwargs):
    """
    Show a figure demonstrating when results were reported in the blind test,
    or why they were not reported.

    :param database:
        A database to connect to.

    :param blind_test_filenames:
        A list of two-length tuples containing: (1) the original filename, and
        (2) the blind test filename.

    :param group_by: [optional]
        Column(s) to group the results by in each panel.

    :param sql_constraint: [optional]
        Any additional SQL constraints to place on the original or blind test
        queries.
    """

    if group_by is not None and isinstance(group_by, (str, unicode)):
        group_by = [group_by]

    # For each blind test object, get the tables data.
    comparison_filename = kwargs.pop("__comparison_filename", None)
    if comparison_filename is not None and os.path.exists(comparison_filename):
        logger.info("Comparisons read from {}".format(comparison_filename))
        comparisons = Table.read(comparison_filename)

    else:
        additional_columns = kwargs.pop(
            "additional_columns", (
                "results.snr", 
                "trim(results.filename) as filename",
                "trim(results.setup) as setup", 
                "trim(results.tech) as tech",
                "ges_fld")
            )
        comparison_rows = []

        for i, (original_filename, blind_filename) in enumerate(blind_test_filenames):

            original_results, blind_results = _retrieve_blind_test_result(
                database, original_filename, blind_filename,
                additional_columns=additional_columns,
                sql_constraint=sql_constraint)

            if original_results is None:
                logger.warn("No original results found for {}".format(
                    original_filename))
                ges_fld = original_filename.split("/")[0]

            else:
                # Match the GES_FLD of the original results to the benchmark one.
                ges_fld = original_results["ges_fld"][0].strip()
            
            match = (truths["GES_FLD"] == ges_fld)
            if not any(match):
                logger.warning("No benchmark record for {}!".format(ges_fld))
                _truths = { "truth_{}".format(p): np.nan for p in parameters}

            else:   
                _truths = \
                    { "truth_{}".format(p): truths[p.upper()][match][0] for p in parameters }


            for row in blind_results:

                comparison_row = _truths.copy()
                comparison_row["snr"] = row["snr"]
                comparison_row.update(
                    { "blind_{}".format(p): row[p] for p in parameters })
                
                if group_by is not None:
                    comparison_row.update({ k: row[k] for k in group_by })

                # Include the original results, if possible.           
                if original_results is None:
                    comparison_row.update({ "original_{}".format(p): np.nan \
                        for p in parameters })

                else:
                    # Match each blind result to the equivalent original result by
                    match_row = (original_results["setup"] == row["setup"]) \
                              * (original_results["node_id"] == row["node_id"])
                    assert match_row.sum() == 1

                    comparison_row.update(
                        { "original_{}".format(p): original_results[p][match_row][0] \
                            for p in parameters })

                comparison_rows.append(comparison_row)

        comparisons = Table(rows=comparison_rows)

    if group_by is not None:
        comparisons = comparisons.group_by(group_by)

    # The number of rows.
    parameter = "teff" # The parameter to determine whether a value was given.
    R = len(comparisons.groups)


    fig, ax = plt.subplots()

    # Assign colors for all segments.
    subset = np.isfinite(comparisons["original_{}".format(parameter)]) \
           * ~np.isfinite(comparisons["blind_{}".format(parameter)])
    unique_reasons = ["Value reported"] \
        + list(set([tech.split("-")[0].strip() for tech in \
            comparisons["blind_tech"][subset]]).difference(["NaN"]))

    cmap = kwargs.get("cmap", mpl.cm.Accent)
    C = float(len(unique_reasons))
    segment_colors = { k: cmap(i/C) for i, k in enumerate(unique_reasons) }

    segment_sums = {}

    G = 0
    ylabels = []
    for i, group in enumerate(comparisons.groups):

        subset = np.isfinite(group["original_{}".format(parameter)])

        blind_reported = np.isfinite(group["blind_{}".format(parameter)][subset])

        segments = OrderedDict([
            ("Value reported", sum(blind_reported))
        ])

        # When it wasn't reported, why was that?
        techs = Counter([tech.split("-")[0].strip() \
            for tech in group["blind_tech"][subset][~blind_reported]])
        sorted(techs, key=techs.get, reverse=True)
        segments.update(techs)

        if "NaN" in segments:
            del segments["NaN"]

        # Convert to fractions.
        total = float(sum(segments.values()))
        if total == 0: continue

        for k in segments.keys():
            segments[k] /= total

            segment_sums.setdefault(k, 0)
            segment_sums[k] += segments[k]

        # Show as a horizontal bar.
        total = 0
        for key, fraction in segments.items():
            ax.add_patch(mpl.patches.Rectangle(
                (total, G), total + fraction, 0.5, 
                facecolor=segment_colors[key], edgecolor="none"))
            total += fraction
        G += 1
        if group_by is not None:
            ylabels.append(" ".join([str(group[column][0]) for column in group_by]))



    for unique_reason in unique_reasons:
        # Don't show things that are less than 1%.
        if segment_sums.get(unique_reason, 0.0) < 0.01: continue
        ax.scatter([-1], [-1], facecolor=segment_colors[unique_reason],
            edgecolor="none", marker="s", s=50, 
            label=unique_reason if len(unique_reason) > 0 else "[None supplied]")
    
    if group_by is not None:
        ax.set_yticks(0.25 + np.arange(G))
        ax.set_yticklabels(ylabels)
        ax.set_ylim(0, G - 0.5)

    ax.set_xlim(0, 1)
    ax.yaxis.set_tick_params(width=0)

    legend = plt.legend(loc="upper center", ncol=3, fontsize=12, frameon=False,
        scatterpoints=1)
    ax.set_ylim(0, G - 0.5 + 1) 

    ax.set_xlabel("Fraction of blind test spectra with results")

    fig.tight_layout()

    return fig

    raise a






def differences_histogram(database, blind_test_filenames,
    parameters=("teff", "logg", "feh"), group_by=None, sql_constraint=None,
    xlims=None, bins=None, **kwargs):
    """
    Show histograms of the differences in results from the exact same spectrum.

    :param database:
        A database to connect to.

    :param blind_test_filenames:
        A list of two-length tuples containing: (1) the original filename, and
        (2) the blind test filename.

    :param parameters:
        A list of parameters to show differences with respect to.

    :param group_by: [optional]
        Column(s) to group the results by in each panel.

    :param sql_constraint: [optional]
        Any additional SQL constraints to place on the original or blind test
        queries.
    """

    if group_by is not None and isinstance(group_by, (str, unicode)):
        group_by = [group_by]

    # For each blind test object, get the tables data.
    comparison_filename = kwargs.pop("__comparison_filename", None)
    if comparison_filename is not None and os.path.exists(comparison_filename):
        logger.info("Comparisons read from {}".format(comparison_filename))
        comparisons = Table.read(comparison_filename)

    else:
        additional_columns = kwargs.pop(
            "additional_columns", (
                "results.snr", 
                "trim(results.filename) as filename",
                "trim(results.setup) as setup", 
                "trim(results.tech) as tech",
                "ges_fld")
            )
        comparison_rows = []

        for i, (original_filename, blind_filename) in enumerate(blind_test_filenames):

            original_results, blind_results = _retrieve_blind_test_result(
                database, original_filename, blind_filename,
                additional_columns=additional_columns,
                sql_constraint=sql_constraint)

            if original_results is None:
                logger.warn("No original results found for {}".format(
                    original_filename))
                ges_fld = original_filename.split("/")[0]

            else:
                # Match the GES_FLD of the original results to the benchmark one.
                ges_fld = original_results["ges_fld"][0].strip()
            
            match = (truths["GES_FLD"] == ges_fld)
            if not any(match):
                logger.warning("No benchmark record for {}!".format(ges_fld))
                _truths = { "truth_{}".format(p): np.nan for p in parameters}

            else:   
                _truths = \
                    { "truth_{}".format(p): truths[p.upper()][match][0] for p in parameters }


            for row in blind_results:

                comparison_row = _truths.copy()
                comparison_row["snr"] = row["snr"]
                comparison_row.update(
                    { "blind_{}".format(p): row[p] for p in parameters })
                
                if group_by is not None:
                    comparison_row.update({ k: row[k] for k in group_by })

                # Include the original results, if possible.           
                if original_results is None:
                    comparison_row.update({ "original_{}".format(p): np.nan \
                        for p in parameters })

                else:
                    # Match each blind result to the equivalent original result by
                    match_row = (original_results["setup"] == row["setup"]) \
                              * (original_results["node_id"] == row["node_id"])
                    assert match_row.sum() == 1

                    comparison_row.update(
                        { "original_{}".format(p): original_results[p][match_row][0] \
                            for p in parameters })

                comparison_rows.append(comparison_row)

        comparisons = Table(rows=comparison_rows)

    if group_by is not None:
        comparisons = comparisons.group_by(group_by)

    # Show histograms of each delta.
    N = len(comparisons.groups) 
    K = len(parameters)

    fig, axes = plt.subplots(K, N, figsize=kwargs.get("figsize", (24, 6)))
    axes = np.atleast_2d(axes)

    _xlims = { p: -np.inf for p in parameters }
    for n, group in enumerate(comparisons.groups):

        for k, parameter in enumerate(parameters):
            ax = axes[k, n]

            # Show the real results.
            delta = group["original_{}".format(parameter)] \
                  - group["blind_{}".format(parameter)]

            finite = np.isfinite(delta)
            if any(finite):
                ax.hist(delta[finite], bins=bins.get(parameter, None), facecolor="k", alpha=0.5)

            ax.text(0.05, 0.80, r"$N = {:.0f}$".format(sum(finite)), transform=ax.transAxes)
            ax.text(0.05, 0.65, r"$\mu = {:.2f}$".format(np.nanmean(delta)), transform=ax.transAxes)
            ax.text(0.05, 0.50, r"$\sigma = {:.2f}$".format(np.nanstd(delta)), transform=ax.transAxes)

            if xlims is None:
                _xlims[parameter] = np.max([np.abs(ax.get_xlim()).max(), _xlims[parameter]])

            else:
                ax.set_xlim(-xlims.get(parameter, None), +xlims.get(parameter, None))
            
            if group_by is not None:
                ax.set_title(" ".join([str(group[k][0]) for k in group_by]))


    for k, parameter in enumerate(parameters):
        for n in range(N):
            ax = axes[k, n]
            ax.xaxis.set_major_locator(MaxNLocator(4))
            
            if xlims is None:
                ax.set_xlim(-_xlims[parameter], +_xlims[parameter])

            ax.set_xlabel("delta {}".format(parameter))
            
            if ax.is_first_col():
                ax.set_ylabel("N")
            else:
                ax.set_yticklabels([])

    fig.tight_layout()

    return fig


def differences_wrt_snr(database, truths, blind_test_filenames, 
    parameters=("teff", "logg", "feh"), group_by=None, sql_constraint=None, 
    xlims=None, ylims=None,
    **kwargs):
    """
    Plot the differences between the accepted value and the reported value for
    the blind test results.

    :param database:
        A database to connect to.

    :param truths:
        A table of accepted benchmark values. These will be matched against the
        database (by GES_FLD) to assign a known value.

    :param blind_test_filenames:
        A list of two-length tuples containing: (1) the original filename, and
        (2) the blind test filename.

    :param parameters:
        A list of parameters to show differences with respect to.

    :param group_by: [optional]
        Column(s) to group the results by in each panel.

    :param sql_constraint: [optional]
        Any additional SQL constraints to place on the original or blind test
        queries.
    """

    if group_by is not None and isinstance(group_by, (str, unicode)):
        group_by = [group_by]


    additional_columns = kwargs.pop(
        "additional_columns", (
            "results.id",
            "results.snr", 
            "trim(results.filename) as filename",
            "trim(results.setup) as setup", 
            "trim(results.tech) as tech",
            "ges_fld")
        )

    # For each blind test object, get the tables data.
    comparison_filename = kwargs.pop("__comparison_filename", None)
    if comparison_filename is not None and os.path.exists(comparison_filename):
        logger.info("Comparisons read from {}".format(comparison_filename))
        comparisons = Table.read(comparison_filename)

    else:
        comparison_rows = []

        for i, (original_filename, blind_filename) in enumerate(blind_test_filenames):

            original_results, blind_results = _retrieve_blind_test_result(
                database, original_filename, blind_filename,
                additional_columns=additional_columns,
                sql_constraint=sql_constraint)

            if original_results is None:
                logger.warn("No original results found for {}".format(
                    original_filename))
                ges_fld = original_filename.split("/")[0]

            else:
                # Match the GES_FLD of the original results to the benchmark one.
                ges_fld = original_results["ges_fld"][0].strip()
            
            match = (truths["GES_FLD"] == ges_fld)
            if not any(match):
                logger.warning("No benchmark record for {}!".format(ges_fld))
                _truths = { "truth_{}".format(p): np.nan for p in parameters}

            else:   
                _truths = \
                    { "truth_{}".format(p): truths[p.upper()][match][0] for p in parameters }


            for row in blind_results:

                comparison_row = _truths.copy()
                comparison_row["snr"] = row["snr"]
                comparison_row["blind_tech"] = row["tech"]
                comparison_row["blind_id"] = row["id"]
                comparison_row.update(
                    { "blind_{}".format(p): row[p] for p in parameters })
                
                if group_by is not None:
                    comparison_row.update({ k: row[k] for k in group_by })

                # Include the original results, if possible.           
                if original_results is None:
                    comparison_row.update({ "original_{}".format(p): np.nan \
                        for p in parameters })

                else:
                    # Match each blind result to the equivalent original result by
                    match_row = (original_results["setup"] == row["setup"]) \
                              * (original_results["node_id"] == row["node_id"])
                    assert match_row.sum() == 1

                    comparison_row.update(
                        { "original_{}".format(p): original_results[p][match_row][0] \
                            for p in parameters })

                comparison_rows.append(comparison_row)

        comparisons = Table(rows=comparison_rows)


    if group_by is not None:
        comparisons = comparisons.group_by(group_by)

    # Only consider groups with finite values in *something*
    skip_groups = []
    for i, group in enumerate(comparisons.groups):
        for parameter in parameters:
            values = group["original_{}".format(parameter)] \
                   * group["truth_{}".format(parameter)] \
                   * group["blind_{}".format(parameter)] \
                   * group["snr"]
            if np.any(np.isfinite(values)): break
        else:
            # Didn't break: none are finite.
            #if group["wg"][0] == 10:
            #    raise a

            #skip_groups.append(i)
            print("Adding group {} {}".format(i, group["name"][0]))

    N = len(comparisons.groups) - len(skip_groups)
    K = len(parameters)

    fig, axes = plt.subplots(K, N, figsize=kwargs.get("figsize", (24, 6)))
    axes = np.atleast_2d(axes)

    _xlims = (np.inf, -np.inf)
    _ylims = { p: -np.inf for p in parameters }
    n = 0
    for g, group in enumerate(comparisons.groups):

        if g in skip_groups:
            continue

        for k, parameter in enumerate(parameters):
            ax = axes[k, n]

            # Show the real results.
            x = group["snr"]
            y_original = group["original_{}".format(parameter)] \
                       - group["truth_{}".format(parameter)]

            y_blind = group["blind_{}".format(parameter)] \
                    - group["truth_{}".format(parameter)]


            ax.scatter(x, y_blind, facecolor="r", edgecolor="r", alpha=0.5)
            ax.scatter(x, y_original, facecolor="k", edgecolor="k", alpha=0.5)

            ax.axhline(0, c="#666666", zorder=-1, linestyle=":")

            if xlims is None:
                _ = np.vstack([ax.get_xlim(), _xlims])
                _xlims = (np.min(_[:, 0]), np.max(_[:, 1]))

            else:
                ax.set_xlim(xlims)

            if ylims is None:
                _ylims[parameter] = np.max([np.abs(ax.get_ylim()).max(), _ylims[parameter]])

            else:
                ax.set_ylim(-ylims.get(parameter, None), +ylims.get(parameter, None))
            
            if group_by is not None:
                ax.set_title(" ".join([str(group[k][0]) for k in group_by]))

        n += 1

    for k, parameter in enumerate(parameters):
        for n in range(N):
            ax = axes[k, n]
            ax.xaxis.set_major_locator(MaxNLocator(4))
            ax.yaxis.set_major_locator(MaxNLocator(3))

            if xlims is None:
                ax.set_xlim(_xlims)
            if ylims is None:
                ax.set_ylim(-_ylims[parameter], +_ylims[parameter])

            if ax.is_last_row():
                ax.set_xlabel("snr")
            else:
                ax.set_xticklabels([])

            if ax.is_first_col():
                ax.set_ylabel(parameter)
            else:
                ax.set_yticklabels([])

    fig.tight_layout()

    return (fig, comparisons)



if __name__ == "__main__":

    from astropy.table import Table

    paths = Table.read("blind-test-giraffe.txt", format="ascii")
    
differences_wrt_snr
