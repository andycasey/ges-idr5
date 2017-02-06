


import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib import gridspec

import corner

import utils

__all__ = ["compare_nodes_within_wg", "compare_to_photometric_teff",
    "compare_to_previous_dr"]


def compare_to_previous_dr(database, wg, node_name, parameter):
    """
    Compare the reported node parameter against the recommended results from
    the previous data release.

    :param database:
        A database for the connection.

    :param wg:
        The working group.

    :param node_name:
        The name of the node.

    :param parameter:
        The parameter to compare against.
    """

    parameter = str(parameter).lower().strip()
    if parameter not in ("teff", "logg", "mh", "xi"):
        raise ValueError("parameter '{}' not recognised".format(parameter))

    node_id = database.retrieve_node_id(wg, node_name)
    results = database.retrieve_table(
        """ SELECT  DISTINCT ON (r.cname) r.cname,
                    r.{0}, r.e_{0},
                    p.{0} as recommended_{0}, p.e_{0} as e_recommended_{0}
            FROM    results r, recommended_idr4 p
            WHERE   r.cname = p.cname
            AND     r.node_id = %s""".format(parameter), (node_id, ))

    labels = {
        "teff": r"$T_{\rm eff}$",
        "logg": r"$\log{g}$",
        "mh":   r"$[{\rm Fe}/{\rm H}]$",
        "xi":   r"$\xi$"
    }

    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 4])

    fig = plt.figure()
    ax_diff = plt.subplot(gs[0])
    ax_main = plt.subplot(gs[1])

    y, x = (results[parameter], results["recommended_{}".format(parameter)])
    yerr, xerr = (
        results["e_{}".format(parameter)],
        results["e_recommended_{}".format(parameter)]
        )
    ax_main.scatter(x, y, facecolor="#666666")
    ax_main.errorbar(x, y, xerr=xerr, yerr=yerr,
        fmt=None, ecolor="#666666", alpha=0.5, zorder=-1)

    # Propagate nans.
    _ = x * y
    limits = [
        np.nanmin(np.hstack([_/x, _/y]).flatten()),
        np.nanmax(np.hstack([_/x, _/y]).flatten())
    ]
    _ = np.ptp(limits)
    limits = [limits[0] - 0.05 * _, limits[1] + 0.05 * _]

    ax_main.plot(limits, limits, linestyle=":", c="#666666", zorder=-100)
    ax_main.set_xlim(limits)
    ax_main.set_ylim(limits)

    ax_main.set_xlabel("{} (iDR4)".format(labels.get(parameter, parameter)))
    ax_main.set_ylabel("{} ({})".format(labels.get(parameter, parameter),
        node_name))
    
    y = y - x
    yerr = np.sqrt(xerr**2 + yerr**2)
    ax_diff.scatter(x, y, facecolor="#666666")
    ax_diff.errorbar(x, y, xerr=xerr, yerr=yerr,
        fmt=None, ecolor="#666666", alpha=0.5, zorder=-1)
    ax_diff.set_xlim(limits)
    ax_diff.set_xticklabels([])
    ax_diff.set_ylabel(r"$\Delta{}$" + labels.get(parameter, parameter))
    limit = np.abs(ax_diff.get_ylim()).max()
    ax_diff.set_ylim(-limit, +limit)
    ax_diff.axhline(0, linestyle=":", c="#666666", zorder=-100)
    ax_diff.yaxis.set_major_locator(MaxNLocator(3))

    return fig


def compare_to_photometric_teff(database, wg, node_name):
    """
    Compare the reported temperatures against the photometric temperatures.

    :param database:
        A database for connections.

    :param wg:
        The working group.

    :param node_name:
        The name of the node.
    """

    node_id = database.retrieve_node_id(wg, node_name)
    results = database.retrieve_table(
        """ SELECT  DISTINCT ON (r.cname) r.cname, r.filename, 
                    r.teff, r.e_teff, s.teff_irfm, s.e_teff_irfm
            FROM    results r, spectra s
            WHERE   r.cname = s.cname
              AND   r.node_id = %s""", (node_id, ))

    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 4])

    fig = plt.figure()
    ax_diff = plt.subplot(gs[0])
    ax_main = plt.subplot(gs[1])

    ax_main.scatter(results["teff_irfm"], results["teff"],
        facecolor="#666666")
    ax_main.errorbar(results["teff_irfm"], results["teff"],
        xerr=results["e_teff_irfm"], yerr=results["e_teff"],
        fmt=None, ecolor="#666666", alpha=0.5, zorder=-1)
    limits = ax_main.get_xlim()
    ax_main.plot(limits, limits, linestyle=":", c="#666666", zorder=-100)
    ax_main.set_xlim(limits)
    ax_main.set_ylim(limits)

    ax_main.set_xlabel(r"$T_{\rm eff,photometry}$ $({\rm K})$")
    ax_main.set_ylabel(r"$T_{\rm eff}$ $({\rm K})$")

    x, y = (results["teff_irfm"], results["teff"] - results["teff_irfm"])
    ax_diff.scatter(x, y, facecolor="#666666")
    ax_diff.errorbar(x, y,
        xerr=results["e_teff_irfm"],
        yerr=np.sqrt(results["e_teff"]**2 + results["e_teff_irfm"]**2),
        fmt=None, ecolor="#666666", alpha=0.5, zorder=-1)
    ax_diff.set_xlim(limits)
    ax_diff.set_xticklabels([])
    ax_diff.set_ylabel(r"$\Delta{}T_{\rm eff}$ $({\rm K})$")
    limit = np.abs(ax_diff.get_ylim()).max()
    ax_diff.set_ylim(-limit, +limit)
    ax_diff.set_yticks([-limit, 0, +limit])
    ax_diff.axhline(0, linestyle=":", c="#666666", zorder=-100)

    fig.tight_layout()

    return fig






def compare_nodes_within_wg(database, wg, parameter, extent=None,
    show_one_to_one=True):
    """
    Show a corner plot comparing all of the nodes within a single working group.

    :param database:
        A database for connections.

    :param wg:
        The working group.

    :param parameter:
        The parameter to compare.

    :param extent: [optional]
        The (lower, upper) limits to show in each axis.

    :param show_one_to_one: [optional]
        Show a dashed line marking the `y=x` relation.
    """

    wg = utils.wg_as_int(wg)
    parameter = parameter.lower()
    if parameter not in ("teff", "logg", "mh", "xi"):
        raise ValueError("parameter '{}' not recognised".format(parameter))

    # Get the nodes.
    nodes = database.retrieve_table("""SELECT id, name FROM nodes
        WHERE wg = %s""", (wg, ))
    N_nodes = len(nodes)

    # Get the data.
    results = database.retrieve_table(
        """SELECT r.node_id, r.cname, r.filename, r.{0}, r.e_{0}, r.feh, r.e_feh
        FROM results r, nodes n 
        WHERE n.wg = %s and n.id = r.node_id
        """.format(parameter),
        (wg, ))

    results = results.group_by("cname")
    N_groups = len(results.groups)

    data = np.nan * np.ones((N_nodes, N_groups))
    error = np.nan * np.ones_like(data)
    node_ids = np.sort(np.unique(results["node_id"]))

    indices = results.groups.indices
    for i, start_index in enumerate(indices[:-1]):
        end_index = indices[i + 1]

        for index in range(start_index, end_index):
            j = np.where(results["node_id"][index] == node_ids)[0][0]

            if parameter == "mh" and \
            not np.any(np.isfinite(results["mh"][index])) \
            and np.any(results["feh"][index]):
                data[j, i] = results["feh"][index]
                error[j, i] = results["e_feh"][index]

            else:
                data[j, i] = results[parameter][index]
                error[j, i] = results["e_{}".format(parameter)][index]

    # Remove axes without any data.
    use = np.any(np.isfinite(data), axis=1)
    data = data[use, :]
    error = error[use, :]
    node_ids = node_ids[use]
    labels = [nodes["name"][nodes["id"] == _][0].strip() for _ in node_ids]

    # How many nodes to plot?
    K = data.shape[0]
    assert K > 0, "Need more than one node to compare against."

    factor = 2.0           # size of one side of one panel
    lbdim = 0.5 * factor   # size of left/bottom margin
    trdim = 0.5 * factor   # size of top/right margin
    whspace = 0.15         # w/hspace size
    plotdim = factor * (K - 1.) + factor * (K - 2.) * whspace
    dim = lbdim + plotdim + trdim

    fig, axes = plt.subplots(K - 1, K - 1, figsize=(dim, dim))
    if 3 > K:
        axes = np.atleast_2d([axes])

    lb = lbdim / dim
    tr = (lbdim + plotdim) / dim
    fig.subplots_adjust(left=lb, bottom=lb, right=tr, top=tr,
        wspace=whspace, hspace=whspace)
    
    # Match all of the nodes
    lim = [np.inf, -np.inf]
    for i in range(1, K):
        for j in range(K):
            if j >= i:
                try:
                    ax = axes[i-1, j]
                except IndexError:
                    continue
                ax.set_visible(False)
                ax.set_frame_on(False)
                continue
            if K > 1:
                ax = axes[i-1, j]
            else:
                ax = axes

            ax.scatter(data[j, :], data[i, :], facecolor="#666666")
            ax.errorbar(data[j, :], data[i, :],
                xerr=error[j, :], yerr=error[i, :],
                fmt=None, ecolor="#666666", alpha=0.5, zorder=-1)

            if ax.is_last_row():
                ax.set_xlabel(labels[j])
            else:
                ax.set_xticklabels([])

            if ax.is_first_col():
                ax.set_ylabel(labels[i])
            else:
                ax.set_yticklabels([])

            if extent is not None:
                if show_one_to_one:
                    ax.plot(extent, extent, c="#666666", zorder=-100)
                ax.set_xlim(extent)
                ax.set_ylim(extent)

                ax.xaxis.set_major_locator(MaxNLocator(5))
                ax.yaxis.set_major_locator(MaxNLocator(5))

            else:
                lim[0] = np.nanmin([ax.get_xlim()[0], ax.get_ylim()[0], lim[0]])
                lim[1] = np.nanmax([ax.get_xlim()[1], ax.get_ylim()[1], lim[1]])

    # Ensure all have the same limits and ticks.
    if extent is None:
        for ax in np.array(axes).flatten():
            if ax.get_visible():
                if show_one_to_one:
                    ax.plot(lim, lim, c="#666666", zorder=-100)

                ax.set_xlim(lim)
                ax.set_ylim(lim)

                ax.xaxis.set_major_locator(MaxNLocator(5))
                ax.yaxis.set_major_locator(MaxNLocator(5))

    return fig
