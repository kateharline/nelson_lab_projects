import pandas as pd
import os

import checkparents as cp

def convert_to_rs(c_and_pos):
    '''
    convert chromosome and start position values into rs# in tassel format
    :param c_and_pos: dataframe of chromosome (int) and position (obj) values
    :return: one dataframe column with rs#s (str)
    '''
    # todo decide best concat

    return

def clean_ints(ints_df):
    '''
    convert ints df from txt into usable information for script
    :param ints_df: dataframe of raw ints data
    :return: new df with only useful info
    '''
    useful_ints = ints_df.loc[:, 'line':'end']
    useful_ints['start'] = convert_to_rs(useful_ints[['chromosome','start']])
    useful_ints['end'] = convert_to_rs(useful_ints['chromosome', 'end'])

    # todo convert

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
    founders, lines = cp.open_founders_and_nils('10NN_CU_with_founders_full.hmp.txt', '10NN_CU_with_founders_full.hmp.txt')

    os.chdir('/Users/kateharline/Desktop/nelson_lab/parent_checker/inputs')
    ints_df = pd.read_csv('10NN_introgressions.txt', sep="\t")

    os.chdir('/Users/kateharline/Desktop/nelson_lab/parent_checker/outputs')
    ints = clean_ints(ints_df)
    int_analysis_df = check_ints(lines, founders, ints)

    # todo display/save output

    return None


if __name__ == '__main__':
    main()