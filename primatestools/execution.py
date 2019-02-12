from pandas import DataFrame

from primatestools.species_storage import SpeciesStorage
from primatestools import mspectools as mstools
import pandas as pd

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
