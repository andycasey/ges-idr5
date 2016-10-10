
import matplotlib.pyplot as plt

from code.model import ensemble, plot


model = ensemble.EnsembleModel.read("homogenisation-wg11-teff.model", None)


#plot.biases(model)

#plot.node_uncertainty_with_snr(model)

plot.node_systematic_uncertainty(model)