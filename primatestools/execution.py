import pandas as pd

from primatestools import mspectools as mstools
from primatestools.plotting import plot_heatmaps

# chimp = SpeciesStorage(name='Chimpanzee',type='ensebml')
# chimp.variant_df = pd.read_csv('pantro_sample.bed',sep='\t',header=None)
# chimp.variant_df.columns = ['CHROM', 'chromStart','chromEnd','name']
# chimp_bed_df = mstools.bed_to_columns(bed_df=chimp.variant_df,polym=False)
m192_dfs = {}
m96_dfs = {}
# [pd.DataFrame(),pd.DataFrame]
for i, name in enumerate(['pantro_sample.bed', 'pantro_sample2.bed']):
    sp_name = 'chimp' + str(i)
    loaded_df = pd.read_csv('pantro_sample.bed', sep='\t', header=None)
    loaded_df.columns = ['CHROM', 'chromStart', 'chromEnd', 'name']
    bed_df = mstools.bed_to_columns(loaded_df, polym=False)
    m192_dfs[sp_name] = mstools.gen192_matrix_from_bed_df(bed_data=loaded_df,
                                                          name=sp_name)

    m96_dfs[sp_name] = mstools.reduce_192_to_96_kelley_notation(
        m192=m192_dfs[sp_name])
print()
m96 = m96_dfs[sp_name][['SUBS', 'Context']]
for sp_name in m96_dfs.keys():
    m96 = pd.merge(m96, m96_dfs[sp_name])
sp_dict = {'chimpanzee': 'chimp0', 'chimpanzee2': 'chimp1'}
setlist = [
    [['chimpanzee', 'chimpanzee2'],
     ['chimpanzee2', 'chimpanzee']]
]
plot_heatmaps(m96=m96, sp_plotname_to_col_name=sp_dict,
              pictures_sets_list=setlist, title='test')
# chimp.m96_common = mstools.reduce_192_to_96_common_notation(m192=m192)
# #chimp.m96 = m96
# context_ab_matrix = pd.read_csv('context_abundance.csv',sep='\t')
# chimp.m192_mut_norm = \
#     mstools.normalize_frequencies_192_context(m192=m192,
#                                               context_abundance_matrix=context_ab_matrix,
#                                               matrix_colname='panTro4')
print()

# plot_mutational_spectra(matrix=m96[],species_colnames=['chimp1'])
