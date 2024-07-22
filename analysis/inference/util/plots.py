import os
os.environ['OPENBLAS_NUM_THREADS'] = '5'
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import warnings
warnings.filterwarnings("ignore")

def plot_mutations(mutdf, run):
    mutdf['rep'] = run['replicate']
    mutdf_obs = mutdf[mutdf['vaf']>0.01].reset_index(drop=True)
    polyp = mutdf_obs[mutdf_obs['group'] == 'Polyp']
    inutero = mutdf_obs[mutdf_obs['group'] == 'inutero']
    Normal = mutdf_obs[mutdf_obs['group'] == 'Normal']

    fig, ax = plt.subplots(1, 1, figsize=(4, 3))

    sns.histplot(data=mutdf_obs, x='vaf', hue='group', bins=np.linspace(0,1,51), alpha=.1, stat='density', common_norm=True, kde=True, ax=ax)
    # sns.kdeplot(data=polyp, x='vaf', hue='rep', fill=True, alpha=.3, common_norm=True)
    sns.despine()
    plt.tight_layout()
    plt.show()

def plot_results(df, mutdf, min_vaf=0.01, outpath=None):
    # Create a figure with 2 plots where the first is 75% of the width
    fig, ax = plt.subplots(1, 1, figsize=(8, 4))

    dfpolyp = df[df['Group'] == 'Polyp'].reset_index(drop=True)
    df = df[df['Group'] != 'Polyp'].reset_index(drop=True)

    # sns.lineplot(data=data, x="t", y=0, linewidth=1, color='black')
    sns.lineplot(data=df, x="t", y='Pop', hue='Group', linewidth=1.5, ax=ax)

    # Plot carry cap
    # plt.axhline(y=parameters['initSize'], color='black', linestyle='--', label='Carrying Capacity')
    ax.axvline(x=0, color='black', linestyle='--', alpha=0.5, zorder=0)

    # Add a legend
    plt.legend(loc='upper left')

    # Create the inset plot
    ax_inset = inset_axes(ax, width="25%", height="45%", loc='center')
    sns.lineplot(data=dfpolyp, x="t", y='Pop', ax=ax_inset, color='black')
    ax_inset.set_title('Polyp')
    ax_inset.set_xlabel('Time (years)')
    ax_inset.set_ylabel('')
    ax_inset.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    
    # Mutations
    mutdf_obs = mutdf[mutdf['vaf']>=min_vaf].reset_index(drop=True)
    mutdf_obs['ppVAF'] = mutdf_obs['vaf'] * 2
    mutdf_obs.loc[mutdf_obs['ppVAF'] > 1, 'ccf'] = 1
    polyp = mutdf_obs[mutdf_obs['group'] == 'Polyp']
    inutero = mutdf_obs[mutdf_obs['group'] == 'inutero']
    Normal = mutdf_obs[mutdf_obs['group'] == 'Normal']
    
    ax_inset2 = inset_axes(ax, width="25%", height="45%", loc='center right')
    # sns.histplot(data=mutdf_obs, x='vaf', hue='group', bins=np.linspace(0,1,51), alpha=.5, stat='density', common_norm=True, kde=True, ax=ax_inset2)
    sns.kdeplot(data=mutdf_obs, x='ppVAF', hue='group', alpha=.5, 
                common_norm=False, fill=True, linewidth=0, ax=ax_inset2,
                bw_adjust=0.5)
    ax_inset2.set_title('')
    ax_inset2.set_xlabel('ppVAF')
    ax_inset2.set_ylabel('Density')
    ax_inset2.set_xlim(0, 1)

    # Scientific notation for y-axis
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax.set_ylabel('Population Size')
    ax.set_xlabel('Time (years)')
    sns.despine()
    plt.tight_layout()
    if outpath is not None:
        plt.savefig(outpath + 'pop_plot.pdf')
    plt.show()