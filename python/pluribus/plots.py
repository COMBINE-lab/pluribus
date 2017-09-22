import metrics
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

def plotGCBias(model):
    num_models = model.obs_.shape[0]

    maxn = float(num_models)
    labels = ["[{},{})".format(100.0*(i/maxn), 100.0*((i+1)/maxn)) for i in xrange(num_models)]
    print(labels)
    fig, ax = plt.subplots(num_models, 1, sharex=True, sharey=True)

    sns.set_style('white')

    for i in range(num_models):
        r = np.log2(model.obs_[i,:] / model.exp_[i,:])
        ax[i].plot(r,  label=labels[i])
        ax[i].axhline(y=0.0, ls='dashed', color='red')
        leg = ax[i].legend(loc='auto', shadow=True)
        leg.get_frame().set_facecolor('green')
        leg.get_frame().set_edgecolor('red')

    sns.despine()
    plt.suptitle('$\log_2$-odds ratio : observed / expected')
    plt.show()

def relDiffPlots(truthCol, predCol, df, bins=None):
    plt.cla()
    if bins:
        df['strata'] = pd.cut(df[truthCol], bins)
    else:
        df['strata'] = pd.cut(df[truthCol], 200)


    rd, ind = metrics.relDiff(truthCol, predCol, df)
    df['rd_{}'.format(predCol)] = rd.abs()
    sns.set(context='paper', style='ticks')
    df.groupby('strata').aggregate(np.mean)['rd_{}'.format(predCol)].dropna().plot(kind='bar', color='k', width=0.9, rot=90)
    sns.despine()
    #plt.show()
    plt.savefig(predCol+'_reldiff.pdf')

def relDiffPairPlot(truthCol, predCol1, predCol2, df, bins=None):
    plt.cla()
    if bins:
        df['strata'] = pd.cut(df[truthCol], bins)
    else:
        df['strata'] = pd.cut(df[truthCol], 200)


    rd1, ind = metrics.relDiff(truthCol, predCol1, df)
    rd2, ind = metrics.relDiff(truthCol, predCol2, df)
    df['rd_{}'.format(predCol1)] = rd1.abs()
    df['rd_{}'.format(predCol2)] = rd2.abs()
    df['rd_{}_vs_{}'.format(predCol1, predCol2)] = (rd1.abs() - rd2.abs())
    sns.set(context='paper', style='ticks')
    df.groupby('strata').aggregate(np.mean)['rd_{}_vs_{}'.format(predCol1, predCol2)].dropna().plot(kind='bar', color='k', width=0.9, rot=90)
    sns.despine()
    #plt.show()
    plt.savefig(predCol1+'_vs_'+predCol2+'_reldiff.pdf')


def relDiffPairScatter(truthCol, predCol1, predCol2, df):
    plt.cla()

    rd1, ind = metrics.relDiff(truthCol, predCol1, df)
    rd2, ind = metrics.relDiff(truthCol, predCol2, df)
    df['rd_{}'.format(predCol1)] = rd1.abs()
    df['rd_{}'.format(predCol2)] = rd2.abs()
    df['rd_{}_vs_{}'.format(predCol1, predCol2)] = (rd1.abs() - rd2.abs())
    sns.set(context='paper', style='ticks')
    plt.scatter(np.log(df[truthCol]+0.1), df['rd_{}_vs_{}'.format(predCol1,predCol2)], color='r', alpha=0.25, label="ard({}) - ard({})".format(predCol1, predCol2))
    #plt.scatter(np.log(df[truthCol]+0.1), df['rd_{}'.format(predCol2)], color='b', alpha=0.25, label=predCol2)
    sns.despine()
    plt.legend()
    plt.savefig(predCol1+'_vs_'+predCol2+'_rd_scatter.pdf')
