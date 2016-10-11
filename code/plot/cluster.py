

import logging
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib import gridspec

import utils

__all__ = ["cluster", "param_vs_param"]

logger = logging.getLogger("ges.idr5.qc")

def param_vs_param(database, wg, node_name, ges_fld, reference_parameter,
    vel_range=None):
    """
    Show parameters vs parameters for stars in a cluster.
    
    :param database:
        A database for transactions.

    :param wg:
        The working group.

    :param node_name:
        The name of the node to show results for.

    :param ges_fld:
        The name of the cluster (as listed in GES_FLD).

    :param reference_parameter:
        The name of the reference parameter to show on the x-axis.

    :param vel_range: [optional]
        The (lower, upper) range of velocities (VEL) to select as bonafide
        cluster members.
    """

    vel_range = vel_range or (-1000, 1000)

    parameters = ["teff", "logg", "mh", "xi"]
    reference_parameter = reference_parameter.lower()
    if reference_parameter not in parameters:
        raise ValueError("unrecognized reference parameter")
    parameters.remove(reference_parameter)

    labels = {
        "teff": r"$T_{\rm eff}$ $({\rm K})$",
        "logg": r"$\log{g}$",
        "feh":   r"$[{\rm Fe}/{\rm H}]$",
        "mh":   r"$[{\rm M}/{\rm H}]$",
        "xi":   r"$\xi$ $({\rm km}$ ${\rm s}^{-1})$"
    }

    # Get the data.
    node_id = database.retrieve_node_id(wg, node_name)

    # Collect the results for this node.
    results = database.retrieve_table(
        """ SELECT  DISTINCT ON (r.cname)
                    r.cname, r.node_id, r.setup, s.filename, s.vel, s.e_vel, 
                    teff, e_teff, logg, e_logg, feh, e_feh, mh, e_mh, xi, e_xi
            FROM    results r, spectra s 
            WHERE   r.cname = s.cname
                AND TRIM(s.ges_fld) = %s
                AND s.vel > %s
                AND s.vel < %s
                AND node_id = %s""",
                (ges_fld, min(vel_range), max(vel_range), node_id))

    if results is None:
        return None

    if reference_parameter.lower() == "mh":
        reference_parameter = utils.mh_or_feh(results)

    fig, axes = plt.subplots(3, 1)
    for i, (ax, parameter) in enumerate(zip(axes, parameters)):

        if ax.is_first_row():
            ax.set_title("({0:.1f}, {1:.1f})".format(*vel_range))

        if parameter == "mh":
            parameter = utils.mh_or_feh(results)

        ax.scatter(results[reference_parameter], results[parameter],
            facecolor="#666666", zorder=1)
        ax.errorbar(results[reference_parameter], results[parameter],
            xerr=results["e_{}".format(reference_parameter)],
            yerr=results["e_{}".format(parameter)],
            fmt=None, ecolor="k", alpha=0.5, zorder=-1)

        if ax.is_last_row():
            ax.set_xlabel(labels.get(reference_parameter, reference_parameter))
        else:
            ax.set_xticklabels([])

        ax.set_ylabel(labels.get(parameter, parameter))

    return fig



def cluster(database, ges_fld, wg, node_name=None, vel_range=None,
    isochrone_filename=None, limit_to_isochrone_range=False, ax=None,
    no_tech_flags=False, show_legend=True, **kwargs):
    """
    Show a Hertzsprung-Russell diagram for cluster stars, with isochrones
    optionally shown.

    :param database:
        A database for transactions.

    :param ges_fld:
        The `GES_FLD` entry for the cluster.

    :param wg:
        The working group.

    :param node_name: [optional]
        The name of the node to show results for. If `None` is provided, then
        recommended results will be shown for the specified working group `wg`.

    :param ges_fld:
        The name of the cluster (as listed in GES_FLD).

    :param vel_range: [optional]
        The (lower, upper) range of velocities (VEL) to select as bonafide
        cluster members.

    :param isochrone_filename: [optional]
        The path of an isochrone_filename file to show for this cluster.

    :param limit_to_isochrone_range: [optional]

    :param ax: [optional]
        Provide an axes to plot the HRD in.

    :param no_tech_flags: [optional]
        Require that all results have no TECH flags.
    """

    vel_range = vel_range or (-1000, 1000)

    vmin = kwargs.pop("vmin", None)
    vmax = kwargs.pop("vmax", None)

    if node_name is not None:
        # Get the data.
        table = "results"
        node_id = database.retrieve_node_id(wg, node_name)
        sql_constraint = "AND r.node_id = '{}'".format(node_id)

    else:
        table = "wg_recommended_results"
        sql_constraint = "AND r.wg = '{}'".format(wg)

    
    if no_tech_flags:
        tech_flag_constraint = " AND TRIM(r.TECH) = ''"

    else:
        tech_flag_constraint = ""

    # Collect the results for this node.
    results = database.retrieve_table(
        """ SELECT  DISTINCT ON (r.cname)
                    r.cname, s.vel, s.e_vel, 
                    r.teff, r.e_teff, r.logg, r.e_logg, r.feh, r.e_feh, r.mh, r.e_mh
            FROM    {table} r, spectra s 
            WHERE   r.cname = s.cname
                AND TRIM(s.ges_fld) = '{ges_fld}'
                AND s.vel > '{lower_vel}'
                AND s.vel < '{upper_vel}'
                {sql_constraint} {tech_flag_constraint}""".format(
                    table=table, ges_fld=ges_fld, 
                    sql_constraint=sql_constraint,
                    tech_flag_constraint=tech_flag_constraint,
                    lower_vel=min(vel_range), upper_vel=max(vel_range)))

    if results is None:
        logger.warn("No cluster data on {} from {}/{}".format(
            ges_fld, wg, node_name))
        return None

    # Draw velocity/metallicity. Highlight members.
    gs = gridspec.GridSpec(3, 1, height_ratios=[1, 12, 4])

    if ax is None:        
        fig = plt.figure()
        ax_hrd = plt.subplot(gs[1])
        ax_diff = plt.subplot(gs[2])

    else:
        ax_diff = None
        ax_hrd = ax
        fig = ax.figure

    mh_col = utils.mh_or_feh(results)

    # Draw HRD, distinguishing markers by setup.
    scat = ax_hrd.scatter(results["teff"], results["logg"], c=results[mh_col],
        label=None, cmap="viridis", vmin=vmin, vmax=vmax)
    ax_hrd.errorbar(results["teff"], results["logg"],
        xerr=results["e_teff"], yerr=results["e_logg"],
        fmt=None, ecolor="k", alpha=0.5, zorder=-1, cmap="viridis",
        label=None)

    if isochrone_filename is not None:
        
        isochrone = utils.parse_isochrone(isochrone_filename)

        label, _ = os.path.splitext(os.path.basename(isochrone_filename))
        label = label.replace("_", "-")
        label = "{0} ({1:.1f}, {2:.1f})".format(label, vel_range[0], vel_range[1])
        ax_hrd.plot(isochrone["teff"], isochrone["logg"],
            c="k", lw=2, zorder=-1, label=label)

        if limit_to_isochrone_range:
            xlimits = (7000, 3000)
            ylimits = (5.5, 0)
            #xlimits = ax_hrd.get_xlim()[::-1]
            #ylimits = ax_hrd.get_ylim()[::-1]

        # Draw different wrt to isochrone.
        x = results["teff"]
        y = []
        for i in range(len(x)):
            distance = np.sqrt((results["teff"][i] - isochrone["teff"])**2 \
                + (1000*(results["logg"][i] - isochrone["logg"]))**2)

            index = np.argmin(distance)
            y.append(results["logg"][i] - isochrone["logg"][index])


        if not limit_to_isochrone_range:
            xlimits = ax_hrd.get_xlim()[::-1]
            ylimits = ax_hrd.get_ylim()[::-1]

        if ax_diff is not None:

            ax_diff.scatter(x, y, c=results[mh_col], cmap="viridis", s=50)
            ax_diff.errorbar(x, y, xerr=results["e_teff"], yerr=results["e_logg"],
                fmt=None, ecolor="k", alpha=0.5, zorder=-1)

            ax_diff.axhline(0, linestyle=":", c="#666666", zorder=-2)

            ax_diff.set_xlabel(r"$T_{\rm eff}$ $({\rm K})$")
            ax_diff.set_ylabel(r"$\Delta\log{g}$")
            #ax_diff.set_xlim(ax_hrd.get_xlim()[::-1])
            ax_diff.set_xlim(xlimits)
            ax_diff.set_ylim(ax_diff.get_ylim()[::-1])
            ax_diff.xaxis.set_major_locator(MaxNLocator(5))
            ax_diff.yaxis.set_major_locator(MaxNLocator(5))


        ax_hrd.set_xticklabels([])

        if show_legend:
            ax_hrd.legend(frameon=False, loc="upper left")

    else:
        ax_hrd.set_xlabel(r"$T_{\rm eff}$ $({\rm K})$")

    # Labels, etc.
    ax_hrd.xaxis.set_major_locator(MaxNLocator(5))
    ax_hrd.yaxis.set_major_locator(MaxNLocator(5))

    ax_hrd.set_ylabel(r"$\log{g}$")

    ax_hrd.set_xlim(xlimits)
    ax_hrd.set_ylim(ylimits)

    if ax is None:
        cb = plt.colorbar(
            cax=plt.subplot(gs[0]), mappable=scat, orientation='horizontal')
        cb.ax.xaxis.set_ticks_position('top')
        cb.ax.xaxis.set_label_position('top')
        cb.set_label(r"$[{\rm Fe}/{\rm H}]$")

    return fig