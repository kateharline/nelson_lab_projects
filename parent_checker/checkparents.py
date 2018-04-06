import pandas as pd
import datetime
import numpy as np

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

def dict_truth(par_dict, line_dict, b73_dict, b73_out):
    truths_dict = {}

    for key, value in line_dict.items():


    return truths_dict


def compare(line, parents, b73_out):
    '''
    for a given line compare its snps to each of a set of parents
    :param line: list or array str of snp calls
    :param parents: list of lists of parent array of str snp calls
    :return: dict of float percentage of snps for each parent { parent : % snps }
    '''
    percent_dict = {}
    truth_array = []

    pars = list(parents)

    b73_truth = make_truth(line, col_to_snp_dict('B73', pars))

    for p in pars:
        truth = make_truth(line, col_to_snp_dict(p, pars))

        if b73_out:
            truth = b73_truth_out(truth, b73_truth)

        percent_parent = sum(truth) / len(truth)
        percent_dict[p] = percent_parent
        truth_array.append(truth)

    return percent_dict, truth_array

def compares(lines, parents, b73_out=False, predicted_parents=None):
    '''
    make comparison between parent and line snps per lines
    :param lines: dataframe of lines and str snp calls
    :param parents: dataframe of str snp calls for parents
    :return: dict of dict of lines : { per parent (str) : percent snp calls (float) } }
    '''
    parents_percents = {}
    max_parents = {}

    lins = list(lines)

    for l in lins:
        line_percents, line_truths = compare(col_to_snp_dict(l, lines), parents, b73_out)

        parents_percents[l] = line_percents
        percent_values = list(line_percents.values())
        max_percent_locs = np.argmax(np.array(percent_values))

        pars = list(line_percents.keys())

        max_parents[l] =

        # # make a graph?
        # takes too long, but would be cool one day
        # if predicted_parents[l] != max_parent:
        #     vis.subplotting(line_truths, pars, l)

    return parents_percents


def main():
    date = str(datetime.date.today())+'_'

    filename = 'unimputed_full'

    lines_file = filename + '.hmp.txt'
    lines_df = pd.read_csv(lines_file, sep="\t")

    founders_file = '10NN_CU_with_founders_full.hmp.txt'
    founders_df = pd.read_csv(founders_file, sep="\t")

    lines = lines_df.loc[:, '10NN0065:D25PAACXX:3:250184947':'10NN0805:C25A5ACXX:7:250185015']
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

    parent_not_b73 = compares(lines, founders, True, pred_par_dict)
    b73_out_df = pd.DataFrame(parent_not_b73)
    b73_out_df.to_csv(date + filename + '_parent_percents_b73_out.csv')

    # pick minimum from the dataframe
    max_ids = b73_out_df.idxmax()
    max_ids.to_csv(date + filename + '_parent_matches.csv')

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