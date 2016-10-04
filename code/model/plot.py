
"""
Plot things relevant to the ensemble model.
"""

import cPickle as pickle
import logging
import numpy as np
import matplotlib.pyplot as plt

import brewer2mpl # For pretty colors.

logger = logging.getLogger("ges")

_DEFAULT_SAVEFIG_KWDS = dict(dpi=300, bbox_inches="tight")

def node_uncertainty_with_snr(model, quartiles=[2.5, 50, 97.5], N=1000,
    show_wg_limit=True, xlims=(1, 500), ylims=(0, 500), hide_node_ids=None, 
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
    ax.set_ylim(*ylims)
    ax.set_xlabel(r"${\rm Signal}$-${\rm to}$-${\rm noise}$ ${\rm ratio},$ $S/N$ $({\rm pixel}^{-1})$")

    default_ylabels = {
        "teff": r"${\rm Uncertainty}$ ${\rm in}$ ${\rm effective}$ ${\rm temperature},$ $\sigma_{T_{\rm eff}}$ $({\rm K})$",
        "logg": r"${\rm Uncertainty}$ ${\rm in}$ ${\rm surface}$ ${\rm gravity},$ $\sigma_\log{g}$ $({\rm dex})$",
        "feh": r"${\rm Uncertainty}$ ${\rm in}$ ${\rm metallicity},$ $\sigma_{\rm [Fe/H]}$ $({\rm dex})$",
    }
    ylabel = kwargs.get("ylabel", default_ylabels[model.parameter])
    ax.set_ylabel(ylabel)
    
    legend_kwds = dict(loc="upper right", frameon=False)
    legend_kwds.update(kwargs.get("legend_kwds", {}))
    plt.legend(**legend_kwds)
    
    ax.set(adjustable="box-forced", aspect="equal")
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