import pandas as pd
import os

import checkparents as cp

def clean_strs(list):
    return [''.join(entry.split()) for entry in list]

def convert_to_rs_ranges(chromosomes, starts, ends):
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

    for i in range(chromosomes):
        rs = []
        for pos in range(starts[i], ends[i]):
            rs.append('S'+str(chromosomes[i])+'_'+str(pos))
        rs_s.append(rs)

    return rs_s

def clean_ints(ints_df):
    '''
    convert ints df from txt into usable information for script
    :param ints_df: dataframe of raw ints data
    :return: new df with only useful info
    '''
    useful_ints = ints_df.loc[:, 'line':'end']
    useful_ints['rs_strings'] = convert_to_rs_range(useful_ints[['chromosome','start','end']])

    # todo convert
    # dict snp : val
    # do lookups of snps

    return useful_ints

def check_ints(lines, founders, ints):
    '''
    run analysis pipleline
    :param lines: dataframe of lines and their snps -- strs
    :param founders: dataframe of founders and their snps -- strs
    :param ints: dataframe of the introgression data
    :return: dataframe of summary stats for introgression parentage
    '''
    # todo pipeline
    return


def main():
    # open snp data
    # todo uncomment to run the script for analysis
    # founders, lines = cp.open_founders_and_nils('10NN_CU_with_founders_full.hmp.txt', '10NN_CU_with_founders_full.hmp.txt')

    os.chdir('/Users/kateharline/Desktop/nelson_lab/parent_checker/inputs')
    ints_df = pd.read_csv('int_test.txt', sep="\t")

    os.chdir('/Users/kateharline/Desktop/nelson_lab/parent_checker/outputs')
    ints = clean_ints(ints_df)
    int_analysis_df = check_ints(lines, founders, ints)

    # todo display/save output

    return None


if __name__ == '__main__':
    main()