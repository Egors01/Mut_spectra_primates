
import pandas as pd

def bed_to_columns(bed_df, polym=True, rev_col=False):
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
    uniq_types = df['SUBS'] + df['Context']
    if len(uniq_types.unique()) != 96:
        Exception('SHAPE ERROR ')
    result[name] = df[name]
    result['SUBS'] = df['SUBS']
    result['Context'] = df['Context']
    return result


def normalize_frequencies_192_context(m192, context_abundance_matrix, matrix_colname):
    m192_m = m192
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
        m192_m.ix[idx, sp_colname] = m192_m.ix[
                                       idx, sp_colname] / denominator

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
    m192m = m192.copy(deep=True)
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

    m96 = normalize_matrix_to_type_sum(matrix=m96)

    return m96

def reduce_192_to_96_common_notation(m192):
    m192m = m192.copy(deep=True)
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

    m96 = normalize_matrix_to_type_sum(matrix=m96)

    return m96
