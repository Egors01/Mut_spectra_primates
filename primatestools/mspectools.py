
import pandas as pd
import seaborn as sns
import re
from pyfaidx import Fasta
import os
from primatestools.constants import MutConstants
import json
import matplotlib.pyplot as plt

def bed_to_columns(bed_df, polym=False, rev_col=False):
    bed_data = bed_df
    bed_data['chromStart'] = bed_df['chromStart'].astype(int)
    bed_data['chromEnd'] = bed_df['chromEnd'].astype(int)
    bed_data['Ref'] = \
    bed_data['name'].str.split('ref:', expand=True)[1].str[0]
    bed_data['Alt'] = \
    bed_data['name'].str.split('alt:', expand=True)[1].str[0]
    bed_data['Context'] = bed_data['name'].str.split('cont:', expand=True)[
                              1].str[0:3]
    bed_data['SUBS'] = bed_data['Ref'].str[:] + '>' + bed_data['Alt'].str[
                                                      :]
    if polym:
        bed_data['AF'] = \
            bed_data['name'].str.split('AF:', expand=True)[1].str.split(
                ":").str[0]
        bed_data['AF'] = bed_data['AF'].astype(float)
        if rev_col:
            bed_data['Rev'] = \
                bed_data['name'].str.split(':arev', expand=True)[1].str[0]
    return bed_data


def gen192_matrix_from_bed_df(bed_data, name='species'):
    result = pd.DataFrame()
    df = bed_data.groupby(['SUBS', 'Context']).size().reset_index(
        name=name)
    df.sort_values(by=['SUBS', "Context"], inplace=True)
    df = df#df[df[name]<500]
    if df.shape[0]==192:
        result[name] = df[name]
        result['SUBS'] = df['SUBS']
        result['Context'] = df['Context']
    else:
        result['Context']= pd.Series(MutConstants().context192)
        result['SUBS'] = pd.Series(MutConstants().subs192)
        result = pd.merge(result, df, left_on=['SUBS', 'Context'],right_on=['SUBS', 'Context'], how='left').fillna(0)
        result = result.sort_values(by=['SUBS', "Context"])
        #result[name] = 0.0
    return result


def normalize_frequencies_192_context(m192, context_abundance_matrix, matrix_colname):

    m192_m = m192[m192.columns.values]
    m192_m['Full_context'] = m192_m['Context'].str[0] + m192_m['SUBS'].str[0] + \
                             m192_m['Context'].str[-1]
    name = matrix_colname
    sp_colname = [col for col in m192_m.columns.values if col != 'SUBS' and col != 'Context'][0]

    for idx in m192_m.index:
        full_context = m192_m.ix[idx, 'Full_context']
        denominator = \
            context_abundance_matrix.loc[
                context_abundance_matrix["Context"] == full_context][
                name].item()
        m192_m.ix[idx, sp_colname] = float(m192_m.ix[
                                       idx, sp_colname]) / float(denominator)

    m192_norm_sum = normalize_matrix_to_type_sum(m192_m)

    return m192_norm_sum


def get_complementary(nucleotide):
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'G':
        return "C"
    elif nucleotide == 'C':
        return "G"
    else:
        return '.'


def reverse_context(context):
    return "".join(
        [get_complementary(x) for x in ''.join(reversed(context))])


def normalize_matrix_to_type_sum( matrix):
    all_context_sum = {}
    matrix_norm_sum = matrix[['SUBS', 'Context']]
    for colname in matrix.columns.values:
        if colname != 'SUBS' and colname != 'Context' and colname != 'Full_context':
            matrix_norm_sum[colname] = matrix[colname] / matrix[
                colname].sum()
            all_context_sum[colname] = matrix[colname].sum()
    matrix_norm_sum = matrix_norm_sum.sort_values(by=['SUBS', 'Context'],
                                                  ascending=True)
    matrix_norm_sum = matrix_norm_sum.reset_index(drop=True)

    return matrix_norm_sum


def reduce_192_to_96_kelley_notation(m192):
    m192m = m192[m192.columns.values]
    #m192m = m192.copy(deep=True)
    colname_to_sp_name = dict(zip(
        [x for x in m192m.columns.values if x != 'SUBS' and x != 'Context'],
        [x.split('_')[0] for x in m192m.columns.values if
         x != 'SUBS' and x != 'Context' and x != 'Full_context']))

    for colname in colname_to_sp_name.keys():
        for idx in m192m.index:
            sbs = m192m.ix[idx, 'SUBS']
            context = m192m.ix[idx, 'Context']
            if sbs == 'T>A':
                m192m.ix[idx, 'SUBS'] = 'A>T'
                m192m.ix[idx, 'Context'] = reverse_context(context)
            elif sbs == 'T>G':
                m192m.ix[idx, 'SUBS'] = 'A>C'
                m192m.ix[idx, 'Context'] = reverse_context(context)
            elif sbs == 'T>C':
                m192m.ix[idx, 'SUBS'] = 'A>G'
                m192m.ix[idx, 'Context'] = reverse_context(context)
            elif sbs == 'G>C':
                m192m.ix[idx, 'SUBS'] = 'C>G'
                m192m.ix[idx, 'Context'] = reverse_context(context)
            elif sbs == 'G>A':
                m192m.ix[idx, 'SUBS'] = 'C>T'
                m192m.ix[idx, 'Context'] = reverse_context(context)
            elif sbs == 'G>T':
                m192m.ix[idx, 'SUBS'] = 'C>A'
                m192m.ix[idx, 'Context'] = reverse_context(context)

    m96 = m192m.groupby(['SUBS', 'Context'], as_index=False).sum()

    m96_res = normalize_matrix_to_type_sum(matrix=m96)

    return m96_res

def reduce_192_to_96_common_notation(m192):
    m192m = m192[m192.columns.values]
    #m192m = m192.copy(deep=True)
    colname_to_sp_name = dict(zip(
        [x for x in m192m.columns.values if x != 'SUBS' and x != 'Context'],
        [x.split('_')[0] for x in m192m.columns.values if
         x != 'SUBS' and x != 'Context' and x != 'Full_context']))

    for colname in colname_to_sp_name.keys():
        for idx in m192m.index:
            sbs = m192m.ix[idx, 'SUBS']
            context = m192m.ix[idx, 'Context']
            if sbs == 'A>T':
                m192m.ix[idx, 'SUBS'] = 'T>A'
                m192m.ix[idx, 'Context'] = reverse_context(context)
            elif sbs == 'G>C':
                m192m.ix[idx, 'SUBS'] = 'C>G'
                m192m.ix[idx, 'Context'] = reverse_context(context)
            elif sbs == 'G>A':
                m192m.ix[idx, 'SUBS'] = 'C>T'
                m192m.ix[idx, 'Context'] =  reverse_context(context)
            elif sbs == 'A>G':
                m192m.ix[idx, 'SUBS'] = 'T>C'
                m192m.ix[idx, 'Context'] = reverse_context(context)
            elif sbs == 'A>C':
                m192m.ix[idx, 'SUBS'] = 'T>G'
                m192m.ix[idx, 'Context'] = reverse_context(context)
            elif sbs == 'G>T':
                m192m.ix[idx, 'SUBS'] = 'C>A'
                m192m.ix[idx, 'Context'] = reverse_context(context)

    m96 = m192m.groupby(['SUBS', 'Context'], as_index=False).sum()

    m96_res = normalize_matrix_to_type_sum(matrix=m96)

    return m96_res

def get_ref_startpos(chrname):
    startpos_df = pd.read_csv('start_positions.csv',sep='\t')
    startpos = startpos_df[chrname].values[0]
    return startpos


def ucsc_get_window_alignment_data(species_name, chr_name, coord_start,
                                   coord_end=None, window_size=None,
                                   raw_fasta_path=None):
    pass
    window_dict = {}
    species_total_len = 0
    species_total_aligned = 0
    aligned_ref = 0
    window_start = coord_start
    if not raw_fasta_path:
        raw_fasta_path = '/uge_mnt/home/Egors01/find_sig_project/project_path/data/chr_raw_fasta_data/'
    ucsc20_sp_names = ['hg38', 'panTro4', 'panPan1', 'gorGor3', 'ponAbe2']

    if coord_end:
        window_end = coord_end
    elif window_size:
        window_size = int(window_size)
        window_end = window_start + window_size
    else:
        print('need to specify window size or coored end')
        raise ValueError

    target_species_name = species_name
    window_dict = {}

    for species_name in ucsc20_sp_names[1:]:
        species_fasta_filename = os.path.join(raw_fasta_path, chr_name,
                                              species_name + '.fasta')
        fasta_file = Fasta(species_fasta_filename, as_raw=True,
                           one_based_attributes=False)
        sequence = fasta_file[species_name][window_start:window_end].upper()
        atgc = re.compile('[ATGCatgc]')
        species_seq_len = len(sequence)
        aligned_len = len(atgc.findall(sequence))

        species_total_len += species_seq_len
        species_total_aligned += aligned_len

        window_dict[species_name] = {}
        if species_seq_len != 0:
            aligned_sp = float(aligned_len) / float(species_seq_len) * 100
        else:
            aligned_sp = 0
        window_dict[species_name]['aligned'] = round(aligned_sp, 3)
        # window_dict[species_name]['seq,aligned']=[species_seq_len,aligned_len]
        if species_name == target_species_name:
            window_dict['target_aligned'] = aligned_sp
    if species_total_len != 0:
        window_dict['total_aligned'] = round(
            float(species_total_aligned) / float(species_total_len) * 100, 4)
    else:
        window_dict['total_aligned'] = 0
    window_dict['sp_len,total_len'] = [species_total_aligned,
                                       species_total_len]
    window_dict['target_name'] = target_species_name
    window_dict['chr'] = chr_name
    window_dict['coord_start'] = window_start
    window_dict['coord_end'] = window_end
    window_dict['target_aligned'] = window_dict[target_species_name]['aligned']
    aligned_target = 0
    species_total_aligned = 0
    species_total_len = 0
    # print chr_name, 'window',window_start,window_end, 'aligned ',round(window_dict['total_aligned'],2)
    # pp.pprint(window_dict)
    return window_dict


def save_list_windows_to_json(windows_list, filename):
    species_chr_json = []
    for window_num, obj in enumerate(windows_list):
        species_chr_json.append(obj)
        for key, value in obj.items():
            if isinstance(value, pd.DataFrame):
                species_chr_json[window_num][key] = value.to_dict()
    with open(filename, 'w') as file:
        file.write(json.dumps(species_chr_json))
    return


def load_windows_for_chrjson(json_filename):
    with open(json_filename) as f:
        loaded = json.load(f, encoding='utf-8')
    loaded_perchr_list = []
    for i, item in enumerate(loaded):
        loaded_perchr_list.append(item)
        for key, value in item.items():
            if key[0] == 'm':
                df = pd.DataFrame(value)
                df.index = df.index.astype(int)
                df = df.sort_index()
                loaded_perchr_list[i][key] = df
    return loaded_perchr_list


def plot_mutational_spectra_vector(m96_vector, vector_name='vector', title='',
                                   subs=[], context=[],
                                   ylim=None):
    m96_context = ['A.A', 'A.C', 'A.G', 'A.T', 'C.A', 'C.C', 'C.G',
                   'C.T', 'G.A', 'G.C', 'G.G', 'G.T', 'T.A', 'T.C',
                   'T.G', 'T.T', 'A.A', 'A.C', 'A.G', 'A.T', 'C.A',
                   'C.C', 'C.G', 'C.T', 'G.A', 'G.C', 'G.G', 'G.T',
                   'T.A', 'T.C', 'T.G', 'T.T', 'A.A', 'A.C', 'A.G',
                   'A.T', 'C.A', 'C.C', 'C.G', 'C.T', 'G.A', 'G.C',
                   'G.G', 'G.T', 'T.A', 'T.C', 'T.G', 'T.T', 'A.A',
                   'A.C', 'A.G', 'A.T', 'C.A', 'C.C', 'C.G', 'C.T',
                   'G.A', 'G.C', 'G.G', 'G.T', 'T.A', 'T.C', 'T.G',
                   'T.T', 'A.A', 'A.C', 'A.G', 'A.T', 'C.A', 'C.C',
                   'C.G', 'C.T', 'G.A', 'G.C', 'G.G', 'G.T', 'T.A',
                   'T.C', 'T.G', 'T.T', 'A.A', 'A.C', 'A.G', 'A.T',
                   'C.A', 'C.C', 'C.G', 'C.T', 'G.A', 'G.C', 'G.G',
                   'G.T', 'T.A', 'T.C', 'T.G', 'T.T']
    m96_subs = ['C>A', 'C>A', 'C>A', 'C>A', 'C>A', 'C>A', 'C>A',
                'C>A', 'C>A', 'C>A', 'C>A', 'C>A', 'C>A', 'C>A',
                'C>A', 'C>A', 'C>G', 'C>G', 'C>G', 'C>G', 'C>G',
                'C>G', 'C>G', 'C>G', 'C>G', 'C>G', 'C>G', 'C>G',
                'C>G', 'C>G', 'C>G', 'C>G', 'C>T', 'C>T', 'C>T',
                'C>T', 'C>T', 'C>T', 'C>T', 'C>T', 'C>T', 'C>T',
                'C>T', 'C>T', 'C>T', 'C>T', 'C>T', 'C>T', 'T>A',
                'T>A', 'T>A', 'T>A', 'T>A', 'T>A', 'T>A', 'T>A',
                'T>A', 'T>A', 'T>A', 'T>A', 'T>A', 'T>A', 'T>A',
                'T>A', 'T>C', 'T>C', 'T>C', 'T>C', 'T>C', 'T>C',
                'T>C', 'T>C', 'T>C', 'T>C', 'T>C', 'T>C', 'T>C',
                'T>C', 'T>C', 'T>C', 'T>G', 'T>G', 'T>G', 'T>G',
                'T>G', 'T>G', 'T>G', 'T>G', 'T>G', 'T>G', 'T>G',
                'T>G', 'T>G', 'T>G', 'T>G', 'T>G']
    matrix = pd.DataFrame(
        {'SUBS': m96_subs, 'Context': m96_context, vector_name: m96_vector})
    c_order = ['C>A', 'C>G', 'C>T',
               'T>A', "T>C",
               "T>G"]

    g = sns.FacetGrid(matrix, col="SUBS", hue="SUBS",
                      col_wrap=6,
                      col_order=c_order)
    g.map(sns.barplot, "Context", vector_name)
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
    plt.show()