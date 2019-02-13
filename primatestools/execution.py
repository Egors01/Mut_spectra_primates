import os

import matplotlib

from pandas import DataFrame

from primatestools.species_storage import SpeciesStorage
from primatestools import mspectools as mstools
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

chimp = SpeciesStorage(name='Chimpanzee',type='ensebml')
chimp.variant_df = pd.read_csv('pantro_sample.bed',sep='\t',header=None)
chimp.variant_df.columns = ['CHROM', 'chromStart','chromEnd','name']
chimp_bed_df = mstools.bed_to_columns(bed_df=chimp.variant_df,polym=False)
m192 = mstools.gen192_matrix_from_bed_df(bed_data=chimp_bed_df,name='Chimp_E')  # type: DataFrame
m96 = mstools.reduce_192_to_96_kelley_notation(m192=m192)
chimp.m96_common = mstools.reduce_192_to_96_common_notation(m192=m192)
#chimp.m96 = m96
context_ab_matrix = pd.read_csv('context_abundance.csv',sep='\t')
chimp.m192_mut_norm = \
    mstools.normalize_frequencies_192_context(m192=m192,
                                              context_abundance_matrix=context_ab_matrix,
                                              matrix_colname='panTro4')
print()


def plot_mutational_spectra( matrix, species_colnames, title='',
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
        if save_to_pdf:
            plt.savefig(pp, format='pdf')
        else:
            plt.show()
    if save_to_pdf:
        pp.close()

plot_mutational_spectra(matrix=m96,species_colnames=['Chimp_E'])