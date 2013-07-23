# -*- coding: utf-8 -*-
"""
analyze growth curves using the modified Gompertz
function from:
1. Zwietering MH, Jongenburger I, Rombouts FM, ’t Riet K van (1990)
Modeling of the bacterial growth curve. Applied and environmental
microbiology 56:1875–81.

y = carcap exp [-exp[linslope*exp[1]/carcap * (lagtime-time)+1]]

lagtime = lambda ~ 100 -
linslope = um ~ 0.05 - linear region slope
carcap = A ~ 1.2 - carrying capacity
"""
import sys
import numpy as np
import pymc as mc
import matplotlib.pyplot as plt
import subprocess
from scipy.stats.mstats import mquantiles
from pymc.Matplot import plot as mcplot


def loaddata(argv):
    od600 = np.loadtxt("data/20130608_40C.csv", delimiter=",",
                              skiprows=1, usecols=range(1,31))
    timemins = np.arange(0, 5*od600.shape[0],5)
    od600red = np.delete(od600,[9,19,29],1) - od600[:,[9,19,29]].mean()

    # add noise for testing
    #od600red = od600red + np.random.normal(0.0, 0.05, size=od600red.shape)

    flatod600 = np.reshape(od600red,od600red.size,order='F') #F = column-major
    flattime = np.resize(timemins,flatod600.size)
    return (od600, timemins, od600red, flatod600, flattime)

def gompertzmod():
    _, _, _, flatod600, flattime = loaddata(None)
    #od600red=np.delete(od600,[9,19,29],1)
    #od600red = od600[:,0:10]
    #flatod600 = np.reshape(od600red,od600red.size,order='F') #F = column-major
    #flattime = np.resize(timemins,flatod600.size)

    # normal hyperparameter priors
    # lagtime = mc.Normal("lagtime", 100, 0.001, value=100)
    # linslope = mc.Normal("linslope", 0.05, 0.001, value=0.05)
    # carcap = mc.Normal("carcap", 1.2, 0.001, value=1.2)

    # uninformative hyperparameter priors
    lagtime = mc.Uninformative("lagtime", value=300)
    linslope = mc.Uninformative("linslope", value=0.005)
    carcap = mc.Uninformative("carcap", value=1.2)

    @mc.deterministic
    def y_mean(lagtime=lagtime, linslope=linslope,
                carcap=carcap, time=flattime):
        return carcap*np.exp(-np.exp(linslope*np.exp(1)/carcap *
                                        (lagtime-time)+1))

    sigma = mc.Uniform('sigma', lower=0, upper=100, value=.05)

    # observation of a single growth curve
    # y_obs = mc.Normal("y_obs", value=od600[:,0].transpose(),
    #                   mu=y_mean, tau=sigma**-2, observed=True)
    # observation of multiple growth curves
    y_obs = mc.Normal("y_obs", value=flatod600,
                      mu=y_mean, tau=sigma**-2, observed=True)
    return vars()

def fit_gompertzmod():
    model = gompertzmod()

    mc.MAP(model).fit(method='fmin_powell')
    try:
        with open('gca.pickle'):
            m = mc.MCMC(model)
    except IOError:
        print "Saving simulation data to new pickle database"
        m = mc.MCMC(model, db='pickle', dbname='gca.pickle')
        #m = mc.MCMC(model, db='sqlite', dbname='gca.sqlite')

    #m = mc.MCMC(model, db='pickle', dbname='gca.pickle')
    #m = mc.MCMC(model, db='sqlite', dbname='gca.sqlite')

    # impose Adaptive Metropolis
    # m.use_step_method(mc.AdaptiveMetropolis,
    #     [m.lagtime, m.linslope, m.carcap, m.sigma])

    # low for testing
    #m.sample(iter=6000, burn=5000, thin=1)
    # medium
    #m.sample(iter=10000, burn=5000, thin=5)
    # long
    m.sample(iter=20000, burn=15000, thin=5)

    m.db.close()
    return m

def decorate_plot():
    pass

def plot_gompertzmod(m, ffname):

    _, time, od600red, _, _ = loaddata(None)
    od600redtrans = od600red.transpose()
    timemat = np.resize(time, od600redtrans.shape[::-1])

    y_mean = []
    fig1 = plt.figure(num=1, figsize=(12, 9), facecolor='w')
    ax1f1 = fig1.add_subplot(111)
    for lagtime, linslope, carcap in zip(m.lagtime.trace(),
                                                m.linslope.trace(),
                                                m.carcap.trace()):
        y = carcap*np.exp(-np.exp(linslope*np.exp(1)/carcap *
                            (lagtime-time)+1))
        # plot growth curve for each sampled parameter set
        #plt.plot(time, y, color="#7A68A6", alpha=.75, zorder=2)
        y_mean.append(y)

    # plot mean curve
    plt.plot(time, np.mean(y_mean, axis=0),
             color="k", ls="--", lw=1,
             label='modified Gompertz growth model')
    decorate_plot()

    # vectorized bottom and top 5% quantiles for "confidence interval"
    quanttest = mquantiles(y_mean, [0.05, 0.95], axis=0)
    plt.fill_between(time, *quanttest, alpha=0.7, color="#7A68A6")
    plt.plot(time, quanttest[0], label="95% CI", color = "#7A68A6", alpha =0.7)
    plt.plot(time, quanttest[1], color = "#7A68A6", alpha =0.7)

    plt.scatter( timemat, od600redtrans,
                    color = "k", s = 5, alpha = 0.5 )
    plt.legend(loc="lower right")
    plt.savefig(ffname, bbox_inches='tight',
        facecolor=fig1.get_facecolor(), edgecolor='none')
    plt.close()

def plot_pardists(m, ffname):
    fig2 = plt.figure(num=2, figsize=(12, 15), facecolor='w')

    ax1f2 = fig2.add_subplot(411)
    ax1f2.set_autoscaley_on(False)
    ax1f2.hist(m.lagtime.trace(), histtype='stepfilled',
        bins = 30, alpha = 0.85, label = "posterior of $\lambda$",
        color = "#467821", normed = True)
    ax1f2.legend(loc = "upper left")
    plt.title(r"Posterior distributions of the variables "
                "$\lambda,\;\mu_m,\;A$")
    #plt.xlim([15,30])
    plt.xlabel("$\lambda$ value")
    plt.ylabel("probability")

    ax2f2 = fig2.add_subplot(412)
    #ax2f2.set_autoscaley_on(False)
    ax2f2.hist(m.linslope.trace(), histtype='stepfilled',
        bins = 30, alpha = 0.85, label = "posterior of $\mu_m$",
        color = "#467821", normed = True)
    ax2f2.legend(loc = "upper left")
    #plt.xlim([15,30])
    plt.xlabel("$\mu_m$ value")
    plt.ylabel("probability")

    ax3f2 = fig2.add_subplot(413)
    #ax3f2.set_autoscaley_on(False)
    ax3f2.hist(m.carcap.trace(), histtype='stepfilled',
        bins = 30, alpha = 0.85, label = "posterior of $A$",
        color = "#467821", normed = True)
    ax3f2.legend(loc = "upper left")
    #plt.xlim([15,30])
    plt.xlabel("$A$ value")
    plt.ylabel("probability")

    ax4f2 = fig2.add_subplot(414)
    #ax4f2.set_autoscaley_on(False)
    ax4f2.hist(m.sigma.trace(), histtype='stepfilled',
        bins = 30, alpha = 0.85, label = "posterior of $\sigma$",
        color = "#467821", normed = True)
    ax4f2.legend(loc = "upper left")
    #plt.xlim([15,30])
    plt.xlabel("$\sigma$ value")
    plt.ylabel("probability")

    fig2.savefig(ffname,bbox_inches='tight',
        facecolor=fig2.get_facecolor(), edgecolor='none')
    plt.close()

def plot_parcorr(m):
    parnames = ["lagtime", "linslope", "carcap", "sigma"]
    for par in parnames:
        mcplot(m.trace(par), common_scale=False)
        plt.savefig(par + ".pdf", bbox_inches='tight', edgecolor='none')
        plt.close()

def gca(argv):
    """
    use MCMC to fit modified Gompertz model parameters to
    growth curve data and plot computed growth curves as well
    as inferred distributions for Gompertz model parameters
    """
    # run MCMC to fit modified Gompertz model
    mcmcoutput = fit_gompertzmod()
    print " "

    # save graphical representation of model
    try:
        with open('gcagraph.pdf'): pass
    except IOError:
        print "Saving graphical representation"
        mc.graph.graph(mcmcoutput, name="gcagraph", format="pdf",
                   prog="dot", legend=True, consts=True)

    # plot growth curves
    print "Plotting growth curves"
    growth_ffname = "growthfit.pdf"
    plot_gompertzmod(mcmcoutput, growth_ffname)
    #evince = subprocess.check_output(["evince",growth_ffname])

    # plot parameter distributions
    #parsplot_ffname = "parfit.pdf"
    #plot_pardists(mcmcoutput, parsplot_ffname)
    #evince = subprocess.check_output(["evince",parsplot_ffname])
    print "Plotting parameter distributions"
    plot_parcorr(mcmcoutput)

    # save parameter estimation data to csv file
    print "Saving parameter estimates"
    pars_ffname = "parests.csv"
    gompertzmodvars = ["lagtime", "linslope", "carcap", "sigma"]
    mcmcoutput.write_csv(pars_ffname, variables=gompertzmodvars)

    print "Combining figures"
    pdfcombine = subprocess.check_output(["pdftk", "growthfit.pdf", "lagtime.pdf", "linslope.pdf", "carcap.pdf", "sigma.pdf", "cat", "output", "allfigs.pdf"])

if __name__ == "__main__":
    #import doctest
    #doctest.testmod()
    gca(sys.argv[1:])
