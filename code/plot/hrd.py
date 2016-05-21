
__all__ = ["hrd"]


import numpy as np
import matplotlib.pyplot as plt



def hrd(database, wg, node_name, where=None, isochrones=None,
    show_other_wg_nodes=False):
    """
    Return a Hertszprung-Russell Diagram (effective temperature versus surface
    gravity) for all stars reported by a single node.

    :param database:
        A database for transactions.

    :param wg:
        The working group.

    :param node_name:
        The name of the node to show results for.

    :param where: [optional]
        Additional SQL constraints, which would start from '... and `where`'.

    :param isochrones: [optional]
        A single path or list of paths of isochrones to display on the figure.

    :param show_other_wg_nodes: [optional]
        If set to `True`, the results from other nodes in the same working
        group will be shown in the background to provide context.
    """

    node_id = database.retrieve_node_id(wg, node_name)

    where_str = "" if where is None else " AND {}".format(where)
    if isochrones is not None and not isinstance(isochrones, (list, tuple)):
        isochrones = [isochrones]

    # Collect the results for this node.
    results = database.retrieve_table(
        """ SELECT node_id, cname, teff, e_teff, logg, e_logg, mh, e_mh, xi, e_xi
            FROM results WHERE node_id = %s {}""".format(where_str), (node_id, ))
    ok = np.isfinite(
        results["teff"].astype(float) * results["logg"].astype(float))

    # TODO Get full limit range on teff/logg/etc?
    #param_ranges = database.retrieve_table(
    #    """SELECT   max(teff) AS max_teff, min(teff) AS min_teff,
    #                max(logg) AS max_logg, min(logg) AS min_logg,
    #                max(mh) AS max_mh, min(mh) AS min_mh,
    #                max(xi) AS max_xi, min(xi) as min_xi
    #       FROM results
    #       WHERE node_id 

    fig, ax = plt.subplots()

    if not np.any(np.isfinite(results["xi"].astype(float))):
        # Don't show size-varying points.
        s = None
    else:
        s = 100 * results["xi"]
        ok *= np.isfinite(results["xi"].astype(float))

    # Error bars.
    ax.errorbar(results["teff"][ok], results["logg"][ok],
        xerr=results["e_teff"][ok], yerr=results["e_logg"][ok],
        fmt=None, c="#666666", zorder=-1)

    scat = ax.scatter(results["teff"][ok], results["logg"][ok],
        c=results["mh"][ok], s=s, cmap="plasma", zorder=2)

    cbar = plt.colorbar(scat)
    cbar.set_label(r"$[{\rm Fe}/{\rm H}]$")

    ax.set_xlim(8000, 3000)
    ax.set_ylim(5.5, 0)
    ax.set_xlabel(r"$T_{\rm eff}$ $({\rm K})$")
    ax.set_ylabel(r"$\log{g}$")

    fig.tight_layout()

    return fig

