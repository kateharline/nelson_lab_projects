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
            chrom = re.match('(?:S)(\d+)', pos_intrs[j])
            if chromosomes[j] == chrom:
                snp_rs = pos_intrs[j]
                # check intr
                position = re.match('(?:_)(\d+)', snp_rs)
                if position >= starts[i] and position <= ends[i]:
                    rs.append(snp_rs)

            # update chrom
            chrom = re.match('(?:S)(\d+)', pos_intrs[i])
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
        line_snps = cp.col_to_snp_dict(line, lines.loc[:, intr])
        current_match_percent = 0.0
        current_match = ''

        parent_dict = {}
        for parent in parents:
            print(intr.dtype)
            parent_snps = cp.col_to_snp_dict(parent, founders.loc[:, intr])
            truth = cp.dict_truth(line_snps, parent_snps, b73_snps, True)
            parent_match_percent = percent_match(truth)
            intr_percents[parent] = parent_match_percent

            if parent_match_percent > current_match_percent:
                current_match = parent
            elif parent_match_percent == current_match_percent:
                current_match += ' '+str(parent)
            parent_dict[parent] = parent_match_percent

        int_tops.append(current_match)
        top_percents.append(current_match_percent)
        intr_percents[i] = parent_dict

    intr_percents = pd.DataFrame(intr_percents)
    intr_percents.to_csv('parent_percents_for_intrs_.csv')


    return int_tops, top_percents


def main():
    # open snp data
    # todo before final run, change to full analysis
    if 'Linux' in platform.platform():
        os.chdir('/local/workdir/parent_checker/inputs')
    elif 'Darwin' in platform.platform():
        os.chdir('/Users/kateharline/Desktop/nelson_lab/parent_checker/inputs')

    founders, lines = cp.open_founders_and_nils('10NN_CML103_lines_with_founders_filtered.hmp.txt', '10NN_CML103_lines_with_founders_filtered.hmp.txt')

    intrs_df = pd.read_csv('int_test.txt', sep="\t")
    print(founders.head(10))
    pos_ints = founders['rs#'].tolist()

    # todo select out possible rs#s

    if 'Linux' in platform.platform():
        os.chdir('/local/workdir/parent_checker/outputs')
    elif 'Darwin' in platform.platform():
        os.chdir('/Users/kateharline/Desktop/nelson_lab/parent_checker/outputs')

    intrs = clean_ints(intrs_df, pos_ints)
    output = pd.DataFrame()
    output['top_parent(s)'], output['percent_match'] = check_ints(lines, founders, intrs)

    # compare syngenta calls
    pred_par_file = '10NN_CU_full_parent_matches.txt'
    predicted_par_df = pd.read_csv(pred_par_file, sep="\t")
    pred_par_dict = dict(zip(predicted_par_df['NIL line'], predicted_par_df['Syngenta called Founder']))
    intrs['Syngenta_call'] = intrs['line'].map(pred_par_dict)

    # todo display/save output
    output.to_csv('intr_call_summary.csv')

    return None


if __name__ == '__main__':
    main()