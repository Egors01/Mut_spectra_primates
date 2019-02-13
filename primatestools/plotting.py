import matplotlib
from itertools import combinations
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns



def plot_heatmaps(m96,sp_plotname_to_col_name,pictures_sets_list,title='insert title'):
    """ Function for plotting heatmaps.

        m96 :param df with 96-types.
        columns [SUBS, Context ,colname1, colname2,...] Row order is important
        sp_plotname_to_col_name: param {species name in plot : columnss name in m96}
        pictures_sets_list:param [ [ [pair00,pair01],[pair10,pair11] ],   ]
        title:param optional
    """
    ratios_df_single = create_ratio_table(m96=m96,
                                               sp_name_to_col_name=sp_plotname_to_col_name)

    for i, pic in enumerate(pictures_sets_list):
        print(i)
        plot_heatmaps_from_ratios_list(ratios_table=ratios_df_single,
                                       pairs_on_image_list=pictures_sets_list[i],
                                       kelley_notation=True,
                                       title=title)
    return

def create_ratio_table(m96, sp_name_to_col_name):
    pairs_to_compare = list(combinations(sp_name_to_col_name.keys(), 2))

    relations_df = m96[['SUBS', 'Context']]
    for i, pair in enumerate(pairs_to_compare):
        first = pair[0]
        second = pair[1]
        colname_first = sp_name_to_col_name[first]
        colname_second = sp_name_to_col_name[second]
        relations_df[first + ' vs ' + second] = float(m96[colname_first]) / float(m96[
            colname_second])
        relations_df[second + ' vs ' + first] = float(m96[colname_second]) / float(m96[
            colname_first])
    return relations_df

def create_heatmap_df_for_species( relations_df, species_column_name,
                                  kelley_notation=True):
    relations = relations_df
    if kelley_notation:
        iterables = [['C>A', 'C>G', 'C>T', 'A>G', 'A>C', 'A>T'],
                     ["T", "G", "C", "A"]]
        index = pd.MultiIndex.from_product(iterables,
                                           names=['subs', "5'- "])

    else:
        iterables = [['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G'],
                     ["A", "T", "G", "C"]]
        index = pd.MultiIndex.from_product(iterables,
                                           names=['subs', "5'- "])

    template_heatmap_df = pd.DataFrame(np.zeros(shape=(24, 4)) - 10,
                                       columns={"A", "C", "G", "T"},
                                       index=index)
    template_heatmap_df = template_heatmap_df[["A", "C", "G", "T"]]
    for idx in relations.index:
        subs = relations.loc[idx, 'SUBS']
        five_prime = relations.loc[idx, 'Context'][0]
        three_prime = relations.loc[idx, 'Context'][2]
        value = relations.loc[idx, species_column_name]
        template_heatmap_df.loc[subs, five_prime][three_prime] = value

    return template_heatmap_df


def plot_heatmaps_from_ratios_list(ratios_table, pairs_on_image_list, title='',
                                   kelley_notation=True, if_save=True,
                                   filename=''):
    colnames_to_plot = [pair[0] + ' vs ' + pair[1] for pair in
                        pairs_on_image_list]
    heatmaps_to_plot = [
        create_heatmap_df_for_species(relations_df=ratios_table,
                                           species_column_name=x,
                                           kelley_notation=kelley_notation)
        for x in colnames_to_plot]

    fig = plt.figure(figsize=(20, 15))
    fig.subplots_adjust(hspace=0.4, wspace=0.4)
    plt.title(title, y=1.08, fontsize=18)

    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    fig.axes[0].get_xaxis().set_visible(False)
    fig.axes[0].get_yaxis().set_visible(False)
    fig.axes[0].set_frame_on(False)
    for i, heatmap in enumerate(heatmaps_to_plot):
        ax = fig.add_subplot(1, len(heatmaps_to_plot), i + 1)
        # ax.text(0.5, 0.5, str((1, 4, i)), fontsize=18, ha='center')
        sns.heatmap(heatmap, annot=True, annot_kws={"size": 9}, linewidths=.6,
                    center=1, vmin=0.5, vmax=1.5, cmap='RdBu_r', ax=ax)
        fig.subplots_adjust(hspace=0.8)
        ax.set_title(pairs_on_image_list[i])

        # fig.subplots_adjust(top=0.5)
    # fig.subplots_adjust(top=.5)
    # fig.suptitle(title, size=12)
    # fig.subplots_adjust(hspace=0.4, wspace=0.4)
    # fig.subplots_adjust(top=.5)
    # plt.tight_layout()
    # plt.subplots_adjust(top=1.15)
    if if_save:
        plt.savefig(filename, dpi=300)
    plt.show()


def plot_mutational_spectra(matrix, species_colnames, title='',
                            kelley_notation=False,
                            save_to_pdf=False,
                            filename="spectrum.pdf",
                            ylim=None):
    if save_to_pdf:
        filename_full = os.path.join(filename)
        pp = matplotlib.backends.backend_pdf.PdfPages(filename_full)
    if kelley_notation:
        c_order = ['C>A', 'C>G', 'C>T', 'A>G',
                   "A>C",
                   "A>T"]
    else:
        c_order = ['C>A', 'C>G', 'C>T',
                   'T>A', "T>C",
                   "T>G"]
    # plt.figure(figsize=(20,20))
    # sns.set(rc={'figure.figsize':(20,20)})
    for species_name in species_colnames:
        if not kelley_notation:
            if len(matrix['SUBS'].values[:]) > 96:
                g = sns.FacetGrid(matrix, col="SUBS", hue="SUBS",
                                  col_wrap=6,
                                  col_order=['C>A', 'C>G', 'C>T',
                                             'T>A', "T>C",
                                             "T>G", 'G>T', 'G>T',
                                             'G>A', 'A>T',
                                             "A>G",
                                             "A>C"])
            else:
                g = sns.FacetGrid(matrix, col="SUBS", hue="SUBS",
                                  col_wrap=6,
                                  col_order=c_order)
        g.map(sns.barplot, "Context", species_name)
        # [x.set_ylim(0.05)
        #  for x in g.axes]
        if ylim:
            g.set(ylim=(0, ylim))
        g.set_xticklabels(rotation=90)
        g.set_titles("{col_name} ")
        g.fig.subplots_adjust(top=.4)
        g.fig.suptitle(title, size=14)
        g.fig.figsize = (15, 30)
        g.fig.set_size_inches(15, 10)
        if save_to_pdf:
            plt.savefig(pp, format='pdf')
        else:

            plt.show()
    if save_to_pdf:
        pp.close()
