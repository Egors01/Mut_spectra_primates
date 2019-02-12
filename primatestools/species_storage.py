import os

import matplotlib
import seaborn as sns
import pandas as pd
class SpeciesStorage():
    def __init__(self,name,assembly_name='',type='',ident=''):
        self.name = name
        self.assembly = assembly_name
        self.type = type
        self.ident = ident
        self.context_abundance = pd.DataFrame()
        self.variant_df = pd.DataFrame()
        self.m96 = pd.DataFrame()
        self.m192 = pd.DataFrame()
        self.m192_mut_norm = pd.DataFrame()
        self.m96_common = pd.DataFrame()

    def load_context_abundance(self,file_path,colname):
        context_abundance = pd.read_csv(file_path)
        species_context_ab_df = context_abundance[colname].to_frame()
        self.context_abundance = species_context_ab_df
        return species_context_ab_df


