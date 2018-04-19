import pandas as pd
import datetime
import numpy as np
import os
import re

import parentvis as vis

import smtplib


### deprecated using lists

def make_truth(line, parent):
    '''
    make truth array of line and parent snp identities
    :param line: list of snps for a given line
    :param parent: list of snps for given parent
    :return: binary truth list
    '''
    return [line[i] == parent[i] for i in range(len(line))]



def b73_truth_out(truth, b73_truth):
    '''
    truth in parent but not in b73
    :param truth: truth matrix with another founder
    :param b73_truth: truth matrix with b73
    :return: new truth matrix with b73 snps removed
    '''
    return [truth[i] and not b73_truth[i] for i in range(len(truth))]

##### using dicts now

def col_to_snp_dict(line, lines):
    '''
    convert columns of lines dataframe into a snp dict
    :param line: str line to make dict of
    :param lines: dataframe of lines with snp info
    :return: dict of strs snp_id : snp
    '''
    return dict(zip(lines['rs#'], lines[line]))

def dict_truth(line_dict, par_dict, b73_dict, b73_out):
    truths_dict = {}

    for key, value in line_dict.items():
        if b73_out:
            truths_dict[key] = value is par_dict[key] and value is not b73_dict[key]
        else:
            truths_dict[key] = value is par_dict[key]

    return truths_dict

def percent_by_chrom(truth, parent, line):
    '''
    calculate the percent truths for each chrom
    :param truth: dict str snp : bool matches parent and not b73
    :param parent:
    :param line:
    :return: float percent value max for chromosome with max
    '''
    chrom_percents = np.zeros(11)
    snps_per_chrom = np.zeros(11)

    for key, value in truth.items():
        chrom = int(re.match('(?:S)(\d+)(?:_)', key)[1])
        snps_per_chrom[chrom] += 1
        if value:
            chrom_percents[chrom] += 1

    chrom_percents = chrom_percents / len(list(truth.values()))
    max_percent = chrom_percents.max()
    tots = chrom_percents.sum()

    return chrom_percents, max_percent, tots

def adjust_max(max_parents, percent_dict, p):
    '''
    helper function to track the parents with max matching values
    :param max_parents: str names of parents w given max percentage
    :param percent_dict: dict str parent : float percent match
    :param p: str current parent
    :return: str maximal parents
    '''
    if percent_dict[p] > percent_dict[max_parents.split(' ')[0]]:
        return str(p) + ' '
    elif percent_dict[p] == percent_dict[max_parents.split(' ')[0]]:
        return max_parents + str(p) + ' '
    else:
        return max_parents

def compare(line, line_name, parents, b73_out, by_chrom):
    '''
    for a given line compare its snps to each of a set of parents
    :param line: list or array str of snp calls
    :param line_name: str line
    :param parents: dataframe of str snp calls
    :param b73_out: bool take what isn't b73
    :param by_chrom: bool choose the max parent match by chromosome or total
    :return: dict of float percentage of snps for each parent { parent : % snps }
    '''
    tot_percent_dict = {'null':0}
    tot_max_parents = 'null '

    chrom_max_dict = {'null': 0}
    chrom_percent_dict = {'null': 0}
    chrom_max_parents = 'null '

    tots_from_chroms = {}

    pars = list(parents)[:-1]

    for p in pars:
        truth = dict_truth(line, col_to_snp_dict(p, parents), col_to_snp_dict('B73', parents), b73_out)

        chrom_percent_dict[p], chrom_max_dict[p], tots_from_chroms[p] = percent_by_chrom(truth, p, line_name)
        tot_percent_dict[p] = list(truth.values()).count(True) / len(list(truth.values()))

        chrom_max_parents = adjust_max(chrom_max_parents, chrom_max_dict, p)
        tot_max_parents = adjust_max(tot_max_parents, tot_percent_dict, p)

    del tot_percent_dict['null']
    del chrom_percent_dict['null']
    del chrom_max_dict['null']

    # percent by parent by chrom
    chrom_percent_df = pd.DataFrame(chrom_percent_dict)
    # chrom_percent_df.to_csv(line_name + '_by_chrom_percents.csv')

    return tot_percent_dict, tot_max_parents, chrom_max_dict, chrom_max_parents

def compares(lines, parents, b73_out=False, predicted_parents=None, by_chrom=True):
    '''
    make comparison between parent and line snps per lines
    :param lines: dataframe of lines and str snp calls
    :param parents: dataframe of str snp calls for parents
    :return: dict of dict of lines : { per parent (str) : percent snp calls (float) } }
    '''
    tot_parents_percents = {}
    tot_max_parents = {}
    chrom_parents_percents = {}
    chrom_max_parents = {}

    lins = list(lines)[:-2]

    for i, l in enumerate(lins):
        tot_parents_percents[l], tot_max_parents[l], chrom_parents_percents[l], chrom_max_parents[l] = compare(col_to_snp_dict(l, lines), l, parents, b73_out, by_chrom)
        print('running '+str(l)+' line '+str(i+1)+' of '+str(len(lins)))

    percents_df = pd.DataFrame.from_dict(tot_parents_percents)
    # percents_df.to_csv('all_chrom_percents.csv')

    return tot_max_parents, chrom_max_parents


def main():
    date = str(datetime.date.today())+'_'

    filename = '10NN_CU_with_founders_full'

    os.chdir('/Users/kateharline/Desktop/nelson_lab/parent_checker/inputs')

    lines_file = filename + '.hmp.txt'
    lines_df = pd.read_csv(lines_file, sep="\t")

    founders_file = '10NN_CU_with_founders_full.hmp.txt'
    founders_df = pd.read_csv(founders_file, sep="\t")

    lines = lines_df.loc[:, '10NN0001':'B73']
    take_max_percent = True
    lines['rs#'] = lines_df.loc[:, 'rs#']
    founders = founders_df.loc[:, 'B73':'Tzi8']
    founders['rs#'] = founders_df.loc[:, 'rs#']

    # parent predicted for each line : parent
    pred_par_file = '10NN_CU_full_parent_matches.txt'

    predicted_par_df = pd.read_csv(pred_par_file, sep="\t")
    pred_par_dict = dict(zip(predicted_par_df['NIL line'], predicted_par_df['Syngenta called Founder']))

    os.chdir('/Users/kateharline/Desktop/nelson_lab/parent_checker/outputs')

    tot_max_parents, chrom_max_parents = compares(lines, founders, True, pred_par_dict, take_max_percent)
    predicted_par_df['tot_match'] = predicted_par_df['NIL line'].map(tot_max_parents)
    predicted_par_df['chrom_match'] = predicted_par_df['NIL line'].map(chrom_max_parents)
    predicted_par_df.to_csv('parentChecker_resultsSummary.csv')

    # email when finished
    # msg = 'done'
    #
    # # Send the message
    # server = smtplib.SMTP('smtp.gmail.com', 587)  # port 465 or 587
    # server.ehlo()
    # server.starttls()s
    # server.ehlo()
    # server.login('harhark8@gmail.com', 'ubomujaxmhprtuwf')
    # server.sendmail('harhark8@gmail.com', 'kh694@cornell.edu', msg)
    # server.close()


if __name__ == '__main__':
    main()