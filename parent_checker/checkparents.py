import pandas as pd
import datetime
import numpy as np
import os

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


def compare(line, parents, b73_out):
    '''
    for a given line compare its snps to each of a set of parents
    :param line: list or array str of snp calls
    :param parents: dataframe of str snp calls
    :param b73_out: bool take what isn't b73
    :return: dict of float percentage of snps for each parent { parent : % snps }
    '''
    percent_dict = {'null':0}
    max_parents = 'null '

    pars = list(parents)

    for p in pars:
        truth = dict_truth(line, col_to_snp_dict(p, parents), col_to_snp_dict('B73', parents), b73_out)

        percent_parent = list(truth.values()).count(True) / len(list(truth.values()))
        percent_dict[p] = percent_parent

        if percent_parent > percent_dict[max_parents.split(' ')[0]]:
            max_parents = str(p)+' '
        elif percent_parent == percent_dict[max_parents.split(' ')[0]]:
            max_parents += str(p)+' '

    del percent_dict['null']

    return percent_dict, max_parents

def compares(lines, parents, b73_out=False, predicted_parents=None):
    '''
    make comparison between parent and line snps per lines
    :param lines: dataframe of lines and str snp calls
    :param parents: dataframe of str snp calls for parents
    :return: dict of dict of lines : { per parent (str) : percent snp calls (float) } }
    '''
    parents_percents = {}
    max_parents = {}

    lins = list(lines)[:-2]

    for i, l in enumerate(lins):
        line_percents, line_max_parents = compare(col_to_snp_dict(l, lines), parents, b73_out)

        parents_percents[l] = line_percents
        max_parents[l] = line_max_parents
        print('running '+str(l)+' line '+str(i)+' of '+str(len(lins)-1))

    return parents_percents, max_parents


def main():
    date = str(datetime.date.today())+'_'

    filename = '10NN_CU_with_founders_full'

    os.chdir('/Users/kateharline/Desktop/nelson_lab/parent_checker/inputs')

    lines_file = filename + '.hmp.txt'
    lines_df = pd.read_csv(lines_file, sep="\t")

    founders_file = '10NN_CU_with_founders_full.hmp.txt'
    founders_df = pd.read_csv(founders_file, sep="\t")

    lines = lines_df.loc[:, '10NN0001':'10NN0013']
    lines['rs#'] = lines_df.loc[:, 'rs#']
    founders = founders_df.loc[:, 'B73':'Tzi8']
    founders['rs#'] = founders_df.loc[:, 'rs#']

    # percent snp identities between founders and lines
    #
    # parent_percent_df = compares(lines, founders)
    # df = pd.DataFrame(parent_percent_df)
    # df.to_csv(date + filename + '_parent_percents.csv')

    # parent predicted for each line : parent
    pred_par_file = '10NN_CU_full_parent_matches.txt'

    predicted_par_df = pd.read_csv(pred_par_file, sep="\t")
    pred_par_dict = dict(zip(predicted_par_df['NIL line'], predicted_par_df['Syngenta called Founder']))

    os.chdir('/Users/kateharline/Desktop/nelson_lab/parent_checker/outputs')

    parent_not_b73, parent_maxs = compares(lines, founders, True, pred_par_dict)
    b73_out_df = pd.DataFrame.from_dict(parent_not_b73, orient='index')

    parent_maxs = pd.DataFrame.from_dict(parent_maxs, orient='index')

    b73_out_df['max_parents'] = parent_maxs[0]

    b73_out_df.to_csv('testing'+ date + filename + '_parent_percents_and_maxs_b73_out.csv')


    # email when finished
    msg = 'done'

    # Send the message
    server = smtplib.SMTP('smtp.gmail.com', 587)  # port 465 or 587
    server.ehlo()
    server.starttls()
    server.ehlo()
    server.login('harhark8@gmail.com', 'ubomujaxmhprtuwf')
    server.sendmail('harhark8@gmail.com', 'kh694@cornell.edu', msg)
    server.close()


if __name__ == '__main__':
    main()