import pandas as pd
import os
import platform
import regex as re

import checkparents as cp

def clean_strs(list):
    '''
    fix poorly formatted entries
    :param list: str of vals with spaces and commas to be separated
    :return: clean list ints
    '''
    list = [''.join(entry.split()) for entry in list]
    list = [''.join(entry.split(',')) for entry in list]
    list = [int(entry) for entry in list]
    return list

def convert_to_rs_ranges(chromosomes, starts, ends, pos_intrs):
    '''
    make rs-like values to search snp space with
    :param chromosome: list of ints chromosome number
    :param start: list of strs that clean to ints start of introgression
    :param end: list of strs clean to ints start of introgression
    :return: list of lists of strs rs#s
    '''
    rs_s = []
    starts = clean_strs(starts)
    ends = clean_strs(ends)


    for i in range(len(starts)):
        rs = []
        for j in range(len(pos_intrs)):
            chrom = re.search('(?:S)(\d+)', pos_intrs[j]).group(1)
            if chromosomes[i] == int(chrom):
                snp_rs = pos_intrs[j]
                # check intr
                position = int(re.search('(?:_)(\d+)', snp_rs).group(1))
                if position >= starts[i] and position <= ends[i]:
                    rs.append(snp_rs)

        rs_s.append(rs)

    return rs_s

def clean_ints(intrs_df, pos_intrs):
    '''
    convert ints df from txt into usable information for script: fix formatting and create rs#s only for snps that exist
    :param ints_df: dataframe of raw ints data
    :return: new df with only useful info
    '''
    useful_intrs = intrs_df.loc[:, 'line':'end']
    useful_intrs['rs_strings'] = convert_to_rs_ranges(useful_intrs['chromosome'].tolist(), useful_intrs['start'].tolist(),
                                                     useful_intrs['end'].tolist(), pos_intrs)

    useful_intrs.to_csv('intrs_w_snps.csv')

    return useful_intrs

def percent_match(truth):
    return sum(truth.values()) / len(truth.values())

def check_ints(lines, founders, intrs):
    '''
    run analysis pipleline
    :param lines: dataframe of lines and their snps -- strs
    :param founders: dataframe of founders and their snps -- strs
    :param ints: dataframe of the introgression data
    :return: dataframe of summary stats for introgression parentage
    '''
    intrs_points = intrs['rs_strings'].tolist()
    parents = list(founders)[:-1]
    int_tops = []
    top_percents = []
    intr_percents = {}

    b73_snps = cp.col_to_snp_dict('B73', founders)

    for i, intr in enumerate(intrs_points):

        line = intrs.loc[i, 'line']
        line_snps = cp.col_to_snp_dict(line, lines.loc[lines['rs#'].isin(intr)])
        current_match_percent = 0.0
        current_match = ''

        parent_dict = {}
        for parent in parents:
            parent_snps = cp.col_to_snp_dict(parent, founders.loc[founders['rs#'].isin(intr)])
            truth = cp.dict_truth(line_snps, parent_snps, b73_snps, True)
            parent_match_percent = percent_match(truth)

            if parent_match_percent > current_match_percent:
                current_match = parent
                current_match_percent = parent_match_percent
            elif parent_match_percent == current_match_percent:
                current_match += ' '+str(parent)
            parent_dict[parent] = parent_match_percent

        int_tops.append(current_match)
        top_percents.append(current_match_percent)
        intr_percents[i] = parent_dict

    intr_percents = pd.DataFrame(intr_percents)
    intr_percents.to_csv('parent_percents_for_intrs.csv')


    return int_tops, top_percents

def add_compare_column(output, output_key, new_col_df, new_col_key, col_to_take, new_col_name):
    '''
    add columns to dataframe to make comparisons across scripts
    :param output: df output from this script
    :param output_key: str column name used to search other df
    :param new_col_df: df adding column from
    :param new_col_key: str column key name in that df
    :param col_to_take: str column name to add
    :param new_col_name: str name of column in output df
    :return: that df with the added column
    '''
    syn_par_dict = dict(zip(new_col_df[new_col_key], new_col_df[col_to_take]))
    output[new_col_name] = output[output_key].map(syn_par_dict)

    return output

def main():
    # open snp data
    # todo before final run, change to full analysis
    if 'Linux' in platform.platform():
        os.chdir('/local/workdir/parent_checker/inputs')
    elif 'Darwin' in platform.platform():
        os.chdir('/Users/kateharline/Desktop/nelson_lab/parent_checker/inputs')

    founders, lines = cp.open_founders_and_nils('10NN_CU_with_founders_full.hmp.txt', '10NN_CU_with_founders_full.hmp.txt')

    intrs_df = pd.read_csv('10NN_introgressions.txt', sep="\t")
    pos_ints = founders['rs#'].tolist()

    predictions_file = 'parentChecker_resultsSummary.csv'
    predictions_df = pd.read_csv(predictions_file)

    # todo select out possible rs#s

    if 'Linux' in platform.platform():
        os.chdir('/local/workdir/parent_checker/outputs')
    elif 'Darwin' in platform.platform():
        os.chdir('/Users/kateharline/Desktop/nelson_lab/parent_checker/outputs')

    intrs = clean_ints(intrs_df, pos_ints)
    output = pd.DataFrame()
    output['int_call'], output['percent_match'] = check_ints(lines, founders, intrs)

    # compare other calls
    output['line'] = intrs['line']
    output = add_compare_column(output, 'line', predictions_df, 'NIL line', 'Syngenta called Founder', 'Syngenta_call')
    output = add_compare_column(output, 'line', predictions_df, 'NIL line', 'tot_match', 'tot_call')
    output = add_compare_column(output, 'line', predictions_df, 'NIL line', 'chrom_match', 'chrom_call')

    # make in nice order
    output = output[['line', 'Syngenta_call', 'int_call', 'percent_match', 'chrom_call', 'tot_call']]

    # todo display/save output
    output.to_csv('intr_call_summary.csv')

    return None


if __name__ == '__main__':
    main()