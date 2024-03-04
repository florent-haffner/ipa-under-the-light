import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme(style="ticks", palette='muted')


def lineplot_shapley(df, fontsize=30, palette=['#585153']):
    _, ax = plt.subplots(1, figsize=(18,10))
    sns.lineplot(df.T, legend=False, palette=palette)
    plt.xlabel("Wavenumbers cm-1", fontsize=fontsize)
    plt.ylabel('Absolute shapley values', fontsize=fontsize)
    plt.show()


def multiplot_shap_line(shap_values, spectra, fontsize=30, nbr_of_ticks=10, palette=['#212c3d', '#c05e31']):    
    _, ax1 = plt.subplots(1, figsize=(18,10))
    ax2 = ax1.twinx()
    
    sns.scatterplot(shap_values.T, legend=True, ax=ax1, alpha=.6, palette=palette)
    ax2 = sns.lineplot(spectra.T, legend=False, ax=ax2, alpha=1., dashes=False, palette=palette)
    
    ax1.set(xlabel='Wavenumbers cm-1', ylabel='Shapley values')
    ax2.set_ylabel('Absorbance')
    plt.show()


def plot_pls_coefficient(values, xlabel, fontsize=30, nbr_of_ticks=25):
    df = pd.DataFrame(values.T)
    df.columns = xlabel.T
    _, ax = plt.subplots(1, figsize=(18,10))

    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(fontsize-7)

    plt.axhline(y=0., color='b', linestyle='--')

    ax.plot(df.transpose())
    plt.xlabel("Wavenumbers cm-1", fontsize=fontsize)
    plt.ylabel('PLS coefficient', fontsize=fontsize)
    plt.show()