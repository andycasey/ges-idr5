
"""
Plot things relevant to the ensemble model.
"""

import cPickle as pickle
import itertools
import logging
import numpy as np
import scipy.sparse

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

import brewer2mpl # For pretty colors.

logger = logging.getLogger("ges")

_DEFAULT_SAVEFIG_KWDS = dict(dpi=300, bbox_inches="tight")

def node_uncertainty_with_snr(model, quartiles=[2.5, 50, 97.5], N=1000,
    show_wg_limit=True, xlims=(1, 500), ylims=None, hide_node_ids=None, 
    show_data_points=True, **kwargs):
    """
    Plot the total node uncertainty as a function of S/N ratio.

    :param model:
        A trained ensemble model.

    :param xlims: [optional]
        A two-length tuple giving the lower and upper bounds to show in SNR.

    :param N: [optional]
        The number of draws to use when calculating the projected node uncertainty.
    """

    quartiles = sum([quartiles], [])
    Q = len(quartiles)
    if Q not in (1, 3):
        raise ValueError("quartiles must be length 1 or 3")

    chains = model._chains

    assert chains is not None, "Has the model been sampled?"


    x = np.linspace(xlims[0], xlims[1], 1000)

    # Get the sample indices.
    K = chains["truths"].shape[0]
    N = K if not K >= N > 0 else N
    indices = np.random.choice(range(K), N, replace=False)

    node_ids = model._data["node_ids"]
    hide_node_ids = hide_node_ids or []


    L = len(node_ids) + 1 if show_wg_limit else len(node_ids)
    colors = brewer2mpl.get_map("Set1", "qualitative", L).mpl_colors[::-1]

    fig, ax = plt.subplots()

    for i, node_id in enumerate(set(node_ids).difference(hide_node_ids)):

        k = np.where(node_id == node_ids)[0][0]

        y = np.sqrt(
            chains["var_sys_estimator"][indices, k].reshape(-1, 1)/np.ones_like(x) + \
            chains["alpha"][indices, k].reshape(-1, 1)**2/x)
        
        #y = np.sqrt(chains["var_sys_estimator"][indices, k].reshape(-1, 1)/np.ones_like(x))
        q = np.percentile(y, quartiles, axis=0)

        # Get node name.
        record = model._database.retrieve(
            "SELECT wg, name FROM nodes WHERE id = %s", (node_id, ))
        assert record is not None, "Node id {} doesn't exist?".format(node_id)
        wg, name = record[0]

        ax.plot(x, q[Q / 2], lw=2, c=colors[i], zorder=10,
            label=r"${{\rm {0}}}$".format(name.strip()))

        # Show filled region.
        if Q > 1:
            ax.fill_between(x, q[0], q[-1], facecolor=colors[i], 
                alpha=0.5, edgecolor="none")


    if show_data_points:

        # The i - 1 is to account for Stan's 1-indexing policy (ew!)
        truths = np.array([model._data["mu_calibrator"][i - 1] \
            for i in model._data["calibrator_index"]])

        data_x = model._data["snr_spectrum"]

        for i, node_id in enumerate(set(node_ids).difference(hide_node_ids)):

            k = np.where(node_id == node_ids)[0][0]
            estimates = model._data["estimates"][:, k] + np.median(chains["c0_estimators"][:, k])

            diffs = estimates - truths
            # Account for filled-in values.
            
            diffs[model._data["var_additive"][:, k] > 0] = np.nan
           
            ax.scatter(data_x, np.abs(diffs), facecolor=colors[i], s=50)

            # Show binned std.
            bins = np.linspace(xlims[0], xlims[1], 25)
            stds = []
            for j in range(bins.size -1):
                in_bin = (bins[j + 1] >= data_x) * (data_x > bins[j])

                value = np.nanstd(diffs[in_bin])
                if 2 > sum(in_bin): value = np.nan
                stds.append(value)

            ax.scatter(bins[:-1] + 0.5 * np.diff(bins)[0], stds, 
                alpha=0.5, lw=2, facecolor="k")
            #ax.plot(bins[:-1] + 0.5 * np.diff(bins)[0], stds, c='k')

    if show_wg_limit:
        # Calculate the minimum variance as a function of SNR using all the
        # information from every node.

        sigma_wg = np.zeros((N, x.size))
        for i, index in enumerate(indices):

            Sigma = np.zeros((x.size, L - 1, L - 1))

            # Fill up the covariance matrix.
            for j in range(L - 1):

                var_total = chains["var_sys_estimator"][index, j]/np.ones_like(x) + \
                            chains["alpha"][index, j].reshape(-1, 1)**2/x
                Sigma[:, j, j] = var_total

            # Fill up the correlation coefficients.
            a = 0
            for j in range(L - 1):
                for k in range(j + 1, L - 1):
                    Sigma[:, j, k] = Sigma[:, k, j] = \
                        chains["rho_estimators"][index, a] * np.sqrt(Sigma[:, j, j] * Sigma[:, k, k])
                    a += 1

            # Calculate the minimum variance.

            for n in range(x.size):

                W = np.ones((L - 1, 1))
                Cinv = np.linalg.inv(Sigma[n])
                sigma_wg[i, n] = np.dot(np.dot(W.T, Cinv), W)**-0.5


        q = np.nanpercentile(sigma_wg, quartiles, axis=0)

        ax.plot(x, q[Q / 2], lw=2, c=colors[-1], linestyle="--",
            label=\
                r"${\rm Homogenised}$"
                "\n"
                r"$({\rm Cram\'er$-${\rm Rao}$ ${\rm bound)}$")
    
        if Q > 1:
            ax.fill_between(x, q[0], q[-1], facecolor=colors[-1],
                alpha=0.5, edgecolor="none", zorder=10)

    ax.set_xlim(*xlims)
    ax.set_xlabel(r"${\rm Signal}$-${\rm to}$-${\rm noise}$ ${\rm ratio},$ $S/N$ $({\rm pixel}^{-1})$")

    default_ylabels = {
        "teff": r"${\rm Uncertainty}$ ${\rm in}$ ${\rm effective}$ ${\rm temperature},$ $\sigma_{T_{\rm eff}}$ $({\rm K})$",
        "logg": r"${\rm Uncertainty}$ ${\rm in}$ ${\rm surface}$ ${\rm gravity},$ $\sigma_\log{g}$ $({\rm dex})$",
        "feh": r"${\rm Uncertainty}$ ${\rm in}$ ${\rm metallicity},$ $\sigma_{\rm [Fe/H]}$ $({\rm dex})$",
    }
    ylabel = kwargs.get("ylabel", default_ylabels.get(model._parameter,""))
    ax.set_ylabel(ylabel)

    if ylims is None:
        default_ylims = dict(teff=(0, 500), logg=(0, 0.5), feh=(0, 0.5))
        ylims = default_ylims.get(model._parameter, None)

    ax.set_ylim(ylims)
    
    legend_kwds = dict(loc="upper right", frameon=False)
    legend_kwds.update(kwargs.get("legend_kwds", {}))
    plt.legend(**legend_kwds)
    
    ax.set(adjustable="box-forced")
    fig.tight_layout()

    # Monkey patch the savefig to incorporate the default keywords to ensure
    # the resulting figure is ready for publication.

    old_savefig = fig.savefig
    def new_savefig(self, *args, **kwargs):
        kwds = _DEFAULT_SAVEFIG_KWDS.copy()
        kwds.update(kwargs)
        return old_savefig(self, *args, **kwds)

    fig.savefig = new_savefig

    return fig



def node_correlations(model, reorder=True, plot_edges=True,
    **kwargs):
    """
    Show a lower-diagonal matrix coloured by the median correlation coefficient
    between each node.

    :param reorder: [optional]
        Re-order the correlation matrix to highlight structure.

    :param plot_edges: [optional]
        Plot edges surrounding the lower-diagonal matrix.
    """

    # Get the node names.
    node_names = []
    for node_id in model._metadata["node_ids"]:
        record = model._database.retrieve(
            "SELECT name FROM nodes WHERE id = {}".format(node_id))
        assert record is not None, "Node id {} is unknown".format(node_id)
        node_names.append(record[0][0].strip())
    node_names = np.array(node_names)
    
    k, N = (0, len(node_names))
    rho = np.nan * np.ones((N, N))

    for i in range(N):
        for j in range(i + 1, N):
            rho[j, i] = np.nanmedian(model._chains["rho_estimators"][:, k])
            rho[i, j] = rho[j, i]
            k += 1

    if reorder:
    
        permutations = list(itertools.permutations(range(N), N))
        score = np.nan * np.ones(len(permutations))
        for i, permutation in enumerate(permutations):

            _ = np.array(permutation)
            m = rho[:, _][_, :]

            score[i] = np.nansum(np.abs(np.diff(m, axis=0))) \
                     + np.nansum(np.abs(np.diff(m, axis=1)))

        matrix_indices = np.array(permutations[np.argmin(score)])

        rho = rho[:, matrix_indices][matrix_indices, :]
        node_names = node_names[matrix_indices]

    rho = rho[1:, :-1]        
    for i in range(N - 1):
        for j in range(i + 1, N - 1):
            rho[i, j] = np.nan        
    

    vrange = np.round(np.nanmax(np.abs(rho)), 1)

    vmin = kwargs.get("vmin", None) or -vrange
    vmax = kwargs.get("vmax", None) or +vrange
    cmap = kwargs.get("cmap", "coolwarm")

    fig, ax = plt.subplots()
    im = ax.imshow(rho, vmin=vmin, vmax=vmax, cmap=cmap, interpolation="nearest")
    
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.xaxis.set_ticks_position("none") # Seriously matplotlib what the fuck
    ax.yaxis.set_ticks_position("none")

    ax.set_xticks(range(N - 1))
    ax.set_yticks(range(N - 1))
    ax.set_xticklabels(node_names[:-1], rotation=90)
    ax.set_yticklabels(node_names[1:])
    
    # Draw a line around everything.
    if plot_edges:
        kwds = dict(lw=2, c="k")

        edge = [-0.5, N - 1 - 0.5]
        ax.plot(edge, [N - 1 - 0.5, N - 1 -0.5], **kwds)
        ax.plot([-0.5, -0.5], edge, **kwds)

        x = []
        y = []
        for i in range(N - 1):
            x.extend([-0.5 + i, i + 0.5, i + 0.5])
            y.extend([i - 0.5, i - 0.5, i + 0.5])

        ax.plot(x, y, **kwds)
        tolerance = kwargs.get("__tolerance", 0.10)
        ax.set_xlim(-0.5 - tolerance, N - 1.5 + tolerance)
        ax.set_ylim(N - 1.5 + tolerance, -0.5 - tolerance)

    ax.set(adjustable="box-forced", aspect="equal")

    p = ax.get_position()
    cbar = plt.colorbar(ax=[ax], mappable=im, orientation="horizontal")
    cbar.set_ticks(np.linspace(-vrange, vrange, 5))
    cbar.ax.xaxis.set_ticks_position("none")

    default_labels = {
        "teff": r"${\rm Correlation}$ ${\rm coefficient},$ $\rho_{T_{\rm eff}}$",
        "logg": r"${\rm Correlation}$ ${\rm coefficient},$ $\rho_\log{g}$",
        "feh": r"${\rm Correlation}$ ${\rm coefficient},$ $\rho_{\rm [Fe/H]}$"
    }
    cbar.set_label(
        kwargs.get("label", default_labels.get(model._parameter, "rho")))

    fig.tight_layout()
    # Monkey patch the savefig to incorporate the default keywords to ensure
    # the resulting figure is ready for publication.

    old_savefig = fig.savefig
    def new_savefig(self, *args, **kwargs):
        kwds = _DEFAULT_SAVEFIG_KWDS.copy()
        kwds.update(kwargs)
        return old_savefig(self, *args, **kwds)

    fig.savefig = new_savefig

    fig.subplots_adjust(bottom=0.35)
    plt.show()

    bottom, height = 0.05, 0.05
    co = kwargs.get("__cbar_offset", 0.01)
    cbar.ax.set_position([
        ax.get_position().x0 + co, 
        bottom, 
        ax.get_position().width - co*2, 
        bottom + height
    ])
    plt.show()

    raise a
    return fig