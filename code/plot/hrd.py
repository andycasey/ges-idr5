
__all__ = ["hrd", "stellar_parameter_histograms", "stellar_parameter_error_histograms"]


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator



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



def histograms(database, wg, node_name, parameters, labels=None, where=None,
    **kwargs):
    """
    Show histograms of parameters.

    :param database:
        A database for transactions.

    :param wg:
        The working group.

    :param node_name:
        The name of the node to show results for.

    :param parameters:
        The names of the columns (parameters) to show.

    :param labels:
        The labels to show. If `None` are given, the `parameters` will be used.

    :param where: [optional]
        Additional SQL constraints, which would start from '... and `where`'.
    """

    node_id = database.retrieve_node_id(wg, node_name)
    where_str = "" if where is None else " AND {}".format(where)

    # Collect the results for this node.
    labels = labels or parameters
    results = database.retrieve_table(
        """SELECT {} FROM results WHERE node_id = %s {}""".format(
            ", ".join(parameters), where_str), (node_id, ))

    N = int(np.ceil(np.sqrt(len(parameters))))
    fig, axes = plt.subplots(N, N)
    axes = np.array(axes).flatten()
    for i, (ax, parameter, label) in enumerate(zip(axes, parameters, labels)):
        ok = np.isfinite(results[parameter])
        ax.hist(results[parameter][ok], bins=50, facecolor="#666666")
        ax.text(0.05, 0.95, r"${}$".format(ok.sum()),
            verticalalignment="top",
            color="k", transform=ax.transAxes)
        ax.set_xlabel(label)
        ax.set_ylabel(r"Count")
        ax.xaxis.set_major_locator(MaxNLocator(5))
        ax.yaxis.set_major_locator(MaxNLocator(5))

    fig.tight_layout()

    for ax in axes[len(parameters):]:
        ax.set_visible(False)

    return fig



def stellar_parameter_histograms(database, wg, node_name, where=None, **kwargs):
    """
    Show four panel historgrams of the TEFF, LOGG, MH, and XI for a given node.

    :param database:
        A database for transactions.

    :param wg:
        The working group.

    :param node_name:
        The name of the node to show results for.

    :param where: [optional]
        Additional SQL constraints, which would start from '... and `where`'.
    """

    parameters = ("teff", "logg", "mh", "xi")
    labels = (r"$T_{\rm eff}$ $({\rm K})$", r"$\log{g}$", 
        r"$[{\rm Fe}/{\rm H}]$", r"$\xi$ $({\rm km}$ ${\rm s}^{-1})$")
    
    return histograms(database, wg, node_name, parameters, labels,
        where=where, **kwargs)



def stellar_parameter_error_histograms(database, wg, node_name, where=None, **kwargs):
    """
    Show four panel historgrams of the E_TEFF, E_LOGG, E_MH, and E_XI for a 
    given node.

    :param database:
        A database for transactions.

    :param wg:
        The working group.

    :param node_name:
        The name of the node to show results for.

    :param where: [optional]
        Additional SQL constraints, which would start from '... and `where`'.
    """

    parameters = ("e_teff", "e_logg", "e_mh", "e_xi")
    labels = (r"$\sigma_{T_{\rm eff}}$ $({\rm K})$", r"$\sigma_{\log{g}}$", 
        r"$\sigma_{[{\rm Fe}/{\rm H}]}$", r"$\sigma_\xi$ $({\rm km}$ ${\rm s}^{-1})$")
    
    return histograms(database, wg, node_name, parameters, labels,
        where=where, **kwargs)

