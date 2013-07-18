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

def loaddata(argv):
    od600 = np.loadtxt("data/20130608_40C.csv", delimiter=",",
                              skiprows=1, usecols=range(1,31))
    timemins = np.arange(0, 5*od600.shape[0],5)
    return od600, timemins

def gompertzmod():
    od600, timemins = loaddata(None)

    # normal hyperparameter priors
    # lagtime = mc.Normal("lagtime", 100, 0.001, value=100)
    # linslope = mc.Normal("linslope", 0.05, 0.001, value=0.05)
    # carcap = mc.Normal("carcap", 1.2, 0.001, value=1.2)

    # uninformative hyperparameter priors
    lagtime = mc.Uninformative("lagtime", value=100)
    linslope = mc.Uninformative("linslope", value=0.05)
    carcap = mc.Uninformative("carcap", value=1.2)

    @mc.deterministic
    def y_mean(lagtime=lagtime, linslope=linslope,
                carcap=carcap, time=timemins):
        return carcap*np.exp(-np.exp(linslope*np.exp(1)/carcap *
                                        (lagtime-time)+1))

    sigma = mc.Uniform('sigma', lower=0, upper=100, value=1.)
    y_obs = mc.Normal("y_obs", value=od600[:,0].transpose(),
                      mu=y_mean, tau=sigma**-2, observed=True)
    return vars()

def fit_gompertzmod():
    model = gompertzmod()

    mc.MAP(model).fit(method='fmin_powell')
    m = mc.MCMC(model)
    #m.sample(iter=10000, burn=5000, thin=5)
    m.sample(iter=6000, burn=5000, thin=1)
    return m

def decorate_plot():
    pass

def plot_gompertzmod(m, ffname):
    od600, time = loaddata(None)
    y_mean = []
    fig1=plt.figure(num=1,figsize=(12,9))
    ax1f1 = fig1.add_subplot(111)
    for lagtime, linslope, carcap, sigma in zip(m.lagtime.trace(),
                                                m.linslope.trace(),
                                                m.carcap.trace(),
                                                m.sigma.trace()):
        y = carcap*np.exp(-np.exp(linslope*np.exp(1)/carcap *
                            (lagtime-time)+1))
        plt.plot(time, y, color='gray', alpha=.75, zorder=-1)
        y_mean.append(y)
    plt.plot(time, np.mean(y_mean, axis=0),
             color='green', label='modified Gompertz growth model')
    decorate_plot()
    plt.savefig(ffname,bbox_inches='tight',
        facecolor=fig1.get_facecolor(), edgecolor='none')
    plt.close()

def plot_pardists(m, ffname):
    fig2=plt.figure(num=2, figsize=(12,15), facecolor='w')

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

def gca(argv):
    mcmcoutput = fit_gompertzmod()

    growth_ffname = "growthfit.pdf"
    #plot_gompertzmod(mcmcoutput, growth_ffname)
    #evince = subprocess.check_output(["evince",growth_ffname])

    pars_ffname = "parfit.pdf"
    plot_pardists(mcmcoutput, pars_ffname)
    evince = subprocess.check_output(["evince",pars_ffname])

if __name__ == "__main__":
    #import doctest
    #doctest.testmod()
    gca(sys.argv[1:])
