

from code.model import ensemble, plot



for wg in (10, ):

    for parameter in ("teff", "logg", "feh"):

        model = ensemble.EnsembleModel.read(
            "homogenisation-wg{}-{}-poly-sys.model".format(wg,parameter), None)

        # Plot the random uncertainty as a function of SNR
        fig = plot.node_uncertainty_with_snr(model)
        
        fig.savefig("figures/wg{wg}/wg{wg}-node-uncertainty-random-{param}.pdf"\
            .format(wg=wg, param=parameter))
        fig.savefig("figures/wg{wg}/wg{wg}-node-uncertainty-random-{param}.png"\
            .format(wg=wg, param=parameter))



        # Plot the distribution of biases for each node
        fig = plot.biases(model)
        fig.savefig("figures/wg{wg}/wg{wg}-biases-{param}.pdf".format(
            wg=wg, param=parameter))
        fig.savefig("figures/wg{wg}/wg{wg}-biases-{param}.png".format(
            wg=wg, param=parameter))
        



        # Plot the relative systematic uncertainty as a function of parameters
        fig = plot.node_relative_systematic_uncertainty(model)
        fig.savefig("figures/wg{wg}/wg{wg}-node-uncertainty-sys-relative-{param}.pdf"\
            .format(wg=wg, param=parameter))
        fig.savefig("figures/wg{wg}/wg{wg}-node-uncertainty-sys-relative-{param}.png"\
            .format(wg=wg, param=parameter))

        # Plot the baseline systematic uncertainty
        fig = plot.systematic_uncertainty(model)
        fig.savefig("figures/wg{wg}/wg{wg}-node-uncertainty-sys-constant-{param}.pdf"\
            .format(wg=wg, param=parameter))
        fig.savefig("figures/wg{wg}/wg{wg}-node-uncertainty-sys-constant-{param}.png"\
            .format(wg=wg, param=parameter))

        # Plot node correlations.
        fig = plot.node_correlations(model)
        fig.savefig("figures/wg{wg}/wg{wg}-correlations-{param}.pdf"\
            .format(wg=wg, param=parameter))
        fig.savefig("figures/wg{wg}/wg{wg}-correlations-{param}.png"\
            .format(wg=wg, param=parameter))