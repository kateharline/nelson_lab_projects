import matplotlib.pyplot as plt
import numpy as np
import os
import datetime

## family colors
founder_colors = {
    'B73': '#0c0c0c',
    'B97': '#f47142',
    'CML103': '#b57a38',
    'CML228': '#ffb049',
    'CML247': '#ffec49',
    'CML277': '#aaff49',
    'CML322': '#49ff9d',
    'CML333': '#49ffdd',
    'CML52': '#49bfff',
    'CML69': '#4970ff',
    'HP301': '#6a49ff',
    'Il14H': '#ec49ff',
    'Ki11': '#ff4994',
    'Ki3': '#ff4964',
    'Ky21': '#961528',
    'M162W': '#961573',
    'M37W': '#681596',
    'MS71': '#301596',
    'Mo17': '#153d96',
    'Mo18W': '#157396',
    'NC350': '#15966f',
    'NC358': '#159637',
    'Oh43': '#629615',
    'Oh7B': '#939615',
    'P39': '#966615',
    'Tx303': '#964415',
    'Tzi8': '#961515'
}

def plot_snps(snps, i, fig, founders):
    '''
    make one subplot
    :param snps: list of x values
    :param i: int plot number
    :param fig: figure object adding subplots to
    :param founder: str founder truths based on
    :return: None
    '''


    plot_color = founder_colors[founders[i]]

    ax = fig.add_subplot(len(founders), 1, i+1)

    ax.set_ylabel(founders[i], color=plot_color, rotation=0, ha='right', size='small', position=(0,0))
    ax.set_yticks([])

    if not ax.is_last_row():
        # all but last
        ax.set_xticks([])

    if ax.is_last_row():
        ax.set_xlabel('snp')
        ax.get_xticklabels()


    ax.set_yticklabels(())

    xs = np.arange(len(snps))

    ax.bar(xs, np.array(snps), color=plot_color)

    return None


def subplotting(x_array, founders, line):
    '''
    make snp plots
    :param x_array: list of lists of ints truth for all founders and given line
    :param founders: list of str founder names
    :param line: str line being tested
    :return: none
    '''

    fig = plt.figure()
    fig.suptitle(line)



    for i, x in enumerate(x_array):
        plot_snps(x, i, fig, founders)
        print('plotting line ' + str(line) + ' parent ' + str(founders[i]))


    date = str(datetime.date.today())+'_'
    fig.savefig('/plots'+ date + str(line)+'_snp_dist.png')

    plt.clf()
    plt.close('all')

    return None

# x_array = np.random.randint(2, size=(27, 200))
# founders = list(founder_colors.keys())
# line = '10NN0049'
#
# subplotting(x_array, founders, line)