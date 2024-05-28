import numpy as np
from scipy.stats import chi2, t
#  wrapper for normfit to avoid using the matlab statistics toolbox.


def normfit_issm(x, alpha=None):
    if alpha is None:
        alpha = 0.05

    #  check for any NaN in any columns
    if not np.isnan(x).any():

        #  explicitly calculate the moments
        muhat = np.mean(x, 0)
    # numpy defaults to 0 delta degrees of freedom; matlab uses 1
        sigmahat = np.std(x, 0, ddof=1)

    # no way to ask this in python, assume 4 outputs
    #if (nargout > 2):
        prob = 1. - alpha / 2.

        if (np.size(x, 0) == 1):
            # operate like matlab normfit, mean, std, etc.
            n = np.size(x)
        else:
            n = np.size(x, 0)

        muci = np.zeros((2, np.size(muhat)))
        sigmaci = np.zeros((2, np.size(sigmahat)))

        try:
            muci[0, :] = muhat - t.ppf(prob, n - 1) * sigmahat / np.sqrt(n)
            muci[1, :] = muhat + t.ppf(prob, n - 1) * sigmahat / np.sqrt(n)
            sigmaci[0, :] = sigmahat * np.sqrt((n - 1) / chi2.ppf(prob, n - 1))
            sigmaci[1, :] = sigmahat * np.sqrt((n - 1) / chi2.ppf(1. - prob, n - 1))
        except:
            muci[0, :] = muhat
            muci[1, :] = muhat
            sigmaci[0, :] = sigmahat
            sigmaci[1, :] = sigmahat
    else:
        #  must loop over columns, since number of elements could be different
        muhat = np.zeros((1, np.size(x, 1)))
        sigmahat = np.zeros((1, np.size(x, 1)))
        muci = np.zeros((2, np.size(x, 1)))
        sigmaci = np.zeros((2, np.size(x, 1)))

    #  remove any NaN and recursively call column
        for j in range(np.shape(x, 1)):
            [muhat[j], sigmahat[j], muci[:, j], sigmaci[:, j]] = normfit_issm(x[not np.isnan(x[:, j]), j], alpha)

    return [muhat, sigmahat, muci, sigmaci]
