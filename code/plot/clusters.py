
"""
Functions to plot an overview of stellar parameter trends across all clusters.
"""

import logging
import numpy as np
import os
import matplotlib.pyplot as plt
import yaml
from matplotlib.ticker import MaxNLocator
from matplotlib import gridspec

import utils

__all__ = ["wg_metallicity_overview", "hertzsprung_russell_diagrams"]

logger = logging.getLogger("ges.idr5.qc")

# HARRIS CATALOG of metallicities
with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), "harris.yaml"), "r") as fp:
    harris_metallicities = yaml.load(fp.read())


def get_valid_clusters(database):

    t = database.retrieve_table(
        """ SELECT  DISTINCT ON (s.ges_fld) TRIM(s.ges_fld) AS ges_fld
              FROM  wg_recommended_results AS r,
                    spectra as s
             WHERE  s.cname = r.cname
               AND  teff <> 'NaN'
               AND  logg <> 'NaN'
               AND  feh <> 'NaN'
               AND  s.vel <> 'NaN'
               AND  s.ges_type IN ('GE_SD_GC', 'GE_SD_OC', 'AR_SD_GC', 'AR_SD_OC', 'GE_CL')
               AND  (s.ges_fld LIKE 'N%' OR s.ges_fld LIKE 'M%')""")

    return list(t["ges_fld"])


def hertzsprung_russell_diagrams(database, cluster_names=None,
    isochrone_filenames=None, velocity_constraints=None, metallicity_constraints=None,
    wgs=None, sql_constraint=None, styles=None, figsize=None, 
    xlim=(7000, 3000), ylim=(5, 0), **kwargs):
    """
    Produce a multi-panel plot showing the Hertszprung-Russell diagram for
    multiple clusters, showing the homogenised results from different
    working groups.

    :param database:
        A database for transactions.

    :param cluster_names: [optional]
        A list of clusters to show.

    :param isochrone_filenames: [optional]
        A dictionary containing the cluster names as keys and an isochrone
        filename as values.

    :param velocity_constraints: [optional]
        A dictionary containing cluster names as keys, and a two-length tuple
        of values specifying the lower and upper radial velocities required for
        cluster membership. If provided, then only stars matching the cluster
        name (by GES_FLD) and within the radial velocity limits will be shown.

    :param metallicity_constraints: [optional]
        A dictionary containing cluster names as keys, and a two-length tuple
        of values specifying the lower and upper metallicities required for
        cluster membership. Cluster candidates with metallicities outside of
        this range will not be shown.

    :param wgs: [optional]
        A tuple of integers specifying the WGs to show. If None is given then
        all WGs will be shown.
    """


    velocity_constraints = velocity_constraints or {}
    
    if cluster_names is None:
        cluster_names = get_valid_clusters(database)

        # Order the cluster names by metallicity.
        cluster_names = sorted(
            cluster_names, key=lambda k: literature.get(k, -np.inf))

    sql_constraint = "" if sql_constraint is None else " AND {}".format(sql_constraint)

    #"#e74c3c"

    default_styles = {
        10: dict(facecolor="#e67e22"), #1abc9c"),
        11: dict(facecolor="#2ecc71",),
        12: dict(facecolor="#3498db",),
        13: dict(facecolor="#9b59b6",),
        20: dict(facecolor="#f1c40f"),#34495e",)
    }
    styles = styles or {}
    for k, v in default_styles.items():
        styles.setdefault(k, v)

    N = len(cluster_names)
    K = int(np.ceil(N**0.5))
    L = K if K*(K - 1) < N else K - 1

    fig, axes = plt.subplots(K, L, figsize=figsize)
    axes = np.array(axes).flatten()

    x_axis, y_axis = ("teff", "logg")

    for i, (ax, cluster_name) in enumerate(zip(axes, cluster_names)):

        # Get the data for this cluster.
        candidates = database.retrieve_table(
            """ SELECT  s.vel, r.wg, 
                        r.teff, r.logg, r.feh, 
                        r.e_teff, r.e_logg, r.e_feh
                  FROM  wg_recommended_results AS r, spectra AS s
                 WHERE  s.cname = r.cname
                   AND  r.teff <> 'NaN' AND r.logg <> 'NaN' AND r.feh <> 'NaN'
                   AND  s.ges_fld = '{}' {}""".format(cluster_name, sql_constraint))


        v_min, v_max = velocity_constraints.get(cluster_name, (-np.inf, +np.inf))
        v_min, v_max = (v_max, v_min) if v_min > v_max else (v_min, v_max)

        f_min, f_max = metallicity_constraints.get(cluster_name, (-np.inf, +np.inf))
        f_min, f_max = (f_max, f_min) if f_min > f_max else (f_min, f_max)

        # Require members to be within some velocity range.
        is_member = (v_max >= candidates["vel"]) * (candidates["vel"] >= v_min) \
                  * (f_max >= candidates["feh"]) * (candidates["feh"] >= f_min)

        members = candidates[is_member].group_by("wg")

        for group in members.groups:
            wg = group["wg"][0]
            if wgs is not None and wg not in wgs: continue

            style_kwds = styles[wg]
            print(cluster_name, wg, len(group))

            y_offset = np.zeros(len(group))

            if wg == 10:
                # WG-15 style corection # MAGIC
                m = 0.5/-2

                y_offset = \
                    np.ones(len(group)) * np.clip(m * (np.nanmedian(group["feh"]) - 0.5), 0, np.inf)
                y_offset[group[y_axis] > 3.5] = 0.0

            ax.scatter(group[x_axis], group[y_axis] + y_offset, s=30, alpha=0.75, 
                **style_kwds)
            ax.errorbar(group[x_axis], group[y_axis] + y_offset,
                xerr=group["e_" + x_axis], yerr=group["e_" + y_axis],
                fmt=None, ecolor="#666666", alpha=0.25, zorder=-1)

        ax.set_title(cluster_name)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        ax.xaxis.set_major_locator(MaxNLocator(5))
        ax.yaxis.set_major_locator(MaxNLocator(5))


    for ax in axes[N:-1]:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_frame_on(False)

    # If final axes is free, use it for legend.
    # If not, show legend outside of all axes.
    if len(axes) > N:
        ax = axes[-1]
        for wg, style_kwds in styles.items():
            if wgs is not None and wg not in wgs: continue
            ax.scatter([0], [0], label="WG{}".format(wg), s=30, **style_kwds)

        ax.set_xlim(1, 2)
        ax.set_ylim(1, 2)
        ax.legend(frameon=False, fontsize=10, loc=4)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_frame_on(False)

    else:
        raise WhatToDo


    for i, ax in enumerate(axes[:N]):
        
        if K > N - i:
            ax.set_xlabel(r"$T_{\rm eff}$ $({\rm K})$")
        else:
            ax.set_xticklabels([])

        if ax.is_first_col():
            ax.set_ylabel(r"$\log{g}$")
        else:
            ax.set_yticklabels([])


    return fig


def wg_metallicity_overview(database, x_axis, literature_metallicities=True,
    cluster_names=None, velocity_constraints=None, metallicity_constraints=None,
    wgs=None, sql_constraint=None, styles=None, xlims=(0, 6), ptp_ylim=1.0, 
    figsize=None):
    """
    Produce a multi-panel plot showing the metallicities (y-axis) for cluster
    stars with respect to a specified stellar parameter (x-axis).

    :param database:
        A database for transactions.

    :param x_axis:
        The name of the stellar parameter to show on the x-axis.

    :param literature_metallicities: [optional]
        Whether to show a horizontal line indicating the literature metallicities.
        This can either be a boolean value, or a dictionary specifying the
        cluster names as keys and the metallicities as values. If a boolean flag
        is given, the Harris (1996) catalog of metallicities will be used.

    :param cluster_names: [optional]
        A list of clusters to show.

    :param velocity_constraints: [optional]
        A dictionary containing cluster names as keys, and a two-length tuple
        of values specifying the lower and upper radial velocities required for
        cluster membership. If provided, then only stars matching the cluster
        name (by GES_FLD) and within the radial velocity limits will be shown.

    :param metallicity_constraints: [optional]
        A dictionary containing cluster names as keys, and a two-length tuple
        of values specifying the lower and upper metallicities required for
        cluster membership. Cluster candidates with metallicities outside of
        this range will not be shown.

    :param wgs: [optional]
        A tuple of integers specifying the WGs to show. If None is given then
        all WGs will be shown.

    """

    y_axis = "feh"
    x_axis, available_x_axis = (x_axis.lower(), ("teff", "logg"))
    if x_axis not in available_x_axis:
        raise ValueError(
            "x_axis must be one of: {}".format(", ".join(available_x_axis)))

    velocity_constraints = velocity_constraints or {}
    literature = harris_metallicities.copy()
    if isinstance(literature_metallicities, dict):
        literature.update(literature_metallicities)
    elif literature_metallicities is False:
        literature = {}

    if cluster_names is None:
        cluster_names = get_valid_clusters(database)

        # Order the cluster names by metallicity.
        cluster_names = sorted(
            cluster_names, key=lambda k: literature.get(k, -np.inf))


    if sql_constraint is None:
        sql_constraint = ""

    else:
        sql_constraint = " AND {}".format(sql_constraint)
    #"",
    #"",
    #"#e74c3c"

    default_styles = {
        10: dict(facecolor="#e67e22"), #1abc9c"),
        11: dict(facecolor="#2ecc71",),
        12: dict(facecolor="#3498db",),
        13: dict(facecolor="#9b59b6",),
        20: dict(facecolor="#f1c40f"),#34495e",)
    }
    styles = styles or {}
    for k, v in default_styles.items():
        styles.setdefault(k, v)

    N = len(cluster_names)
    K = int(np.ceil(N**0.5))
    L = K if K*(K - 1) < N else K - 1

    fig, axes = plt.subplots(K, L, figsize=figsize)
    axes = np.array(axes).flatten()


    for i, (ax, cluster_name) in enumerate(zip(axes, cluster_names)):

        # Get the data for this cluster.
        candidates = database.retrieve_table(
            """ SELECT  s.vel, r.wg, 
                        r.teff, r.logg, r.feh, 
                        r.e_teff, r.e_logg, r.e_feh
                  FROM  wg_recommended_results AS r, spectra AS s
                 WHERE  s.cname = r.cname
                   AND  r.teff <> 'NaN' AND r.logg <> 'NaN' AND r.feh <> 'NaN'
                   AND  s.ges_fld = '{}' {}""".format(cluster_name, sql_constraint))


        v_min, v_max = velocity_constraints.get(cluster_name, (-np.inf, +np.inf))
        v_min, v_max = (v_max, v_min) if v_min > v_max else (v_min, v_max)

        f_min, f_max = metallicity_constraints.get(cluster_name, (-np.inf, +np.inf))
        f_min, f_max = (f_max, f_min) if f_min > f_max else (f_min, f_max)

        # Require members to be within some velocity range.
        is_member = (v_max >= candidates["vel"]) * (candidates["vel"] >= v_min) \
                  * (f_max >= candidates["feh"]) * (candidates["feh"] >= f_min)

        members = candidates[is_member].group_by("wg")

        for group in members.groups:
            wg = group["wg"][0]
            if wgs is not None and wg not in wgs: continue

            style_kwds = styles[wg]
            print(cluster_name, wg, len(group))

            ax.scatter(group[x_axis], group[y_axis], s=30, alpha=0.75, 
                **style_kwds)
            ax.errorbar(group[x_axis], group[y_axis],
                xerr=group["e_" + x_axis], yerr=group["e_" + y_axis],
                fmt=None, ecolor="#666666", alpha=0.5, zorder=-1)

            median_y = np.nanmedian(group[y_axis])
            std_y = np.nanstd(group[y_axis])
            ax.axhline(median_y, 
                c=style_kwds["facecolor"], lw=2, zorder=-1, alpha=0.75)

            ax.axhspan(median_y - std_y, median_y + std_y,
                facecolor=style_kwds["facecolor"], alpha=0.75, zorder=-10,
                edgecolor="none")

        literature_feh = literature.get(cluster_name, np.nan)
        ax.axhline(literature_feh, c="#666666", linestyle=":", zorder=-1)

        #ax.text(0.5, 0.90, cluster_name, transform=ax.transAxes,
        #    horizontalalignment="center", verticalalignment="top")
        ax.set_title(cluster_name)

    for ax in axes[N:-1]:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_frame_on(False)

    # If final axes is free, use it for legend.
    # If not, show legend outside of all axes.
    if len(axes) > N:
        ax = axes[-1]
        for wg, style_kwds in styles.items():
            if wgs is not None and wg not in wgs: continue
            ax.scatter([0], [0], label="WG{}".format(wg), s=30, **style_kwds)

        ax.set_xlim(1, 2)
        ax.set_ylim(1, 2)
        ax.legend(frameon=False, fontsize=10, loc=4)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_frame_on(False)

    else:
        raise WhatToDo


    # Common axes, etc etc.
    if xlims is None:
        limits = np.array([ax.get_xlim() for ax in axes[:N]])
        xlims = (limits.min(), limits.max())

    for i, ax in enumerate(axes[:N]):
        ax.set_xlim(xlims)

        if K > N - i:
            ax.set_xlabel(x_axis)

        else:
            ax.set_xticklabels([])

        if ptp_ylim is not None:
            center = np.mean(ax.get_ylim())
            ax.set_ylim(center - ptp_ylim/2.0, center + ptp_ylim/2.0)

        ax.yaxis.set_major_locator(MaxNLocator(5))
        if ax.is_first_col():
            ax.set_ylabel(r"$[{\rm Fe}/{\rm H}]$")



    for i, (ax, cluster_name) in enumerate(zip(axes, cluster_names)):

        ylim = ax.get_ylim()

        f_min, f_max = metallicity_constraints.get(cluster_name, (-np.inf, +np.inf))
        f_min, f_max = (f_max, f_min) if f_min > f_max else (f_min, f_max)

        if np.all(np.isfinite([f_min, f_max])):
            ax.axhspan(f_min, f_max, facecolor="#CCCCCC", edgecolor="none", 
                alpha=0.5, zorder=-1)

            ax.set_ylim(ylim)

    fig.tight_layout()
    fig.subplots_adjust(wspace=0.10, hspace=0.10)

    raise a
