#!/usr/bin/env python

'''
Compute and plot sensitivity and specificity for various
variant effect prediction tools.

Usage: python graphs.py clinvar.summary.tsv

Note: Sift scores range from 1 to 0, but we convert the range
to 0-1 to make it easier to plot. We need to relabel the X axis

The input data set provides a score for each tool for each
pathogenic and benign variant from ClinVar.

For each tool we consider a range of scores as threshold values.
Variants scored >= the threshold are considered to be classified
by the tool as pathogenic, and the rest are considered to be classified
by the tool as benign.

We can then compute the true positive, false positive, true negative
and false negative rates for each threshold value. From there we
can compute sensitivity and specificity.

For a given threshold value V:

    - true positive (TP) is any ClinVar-pathogenic variant with score >= V
    - false positive (FP) is any ClinVar-benign variant with score >= V
    - true negative (TN) is any ClinVar-benign variant with score < V
    - false negative (FN) is any ClinVar-pathogenic variant with score < V

Percentage of Benign variants classified correctly:

    specificity = TN / (TN + FP)

Percentage of Pathogenic variants classied correctly:

    sensitivity = TP / (TP + FN)
'''

import csv
import sys
import matplotlib.pyplot as plt
from argparse import ArgumentParser

VERSION = '1.0'

def parse_args():
    'Parse the command line arguments for the program.'
    parser = ArgumentParser(
        description="Compute sensitivity and specificity for various variant effect prediction tools")
    parser.add_argument(
        '--version', action='version', version='%(prog)s ' + VERSION)
    parser.add_argument(
        '--dataset', type=str, help='Name of variant truth dataset', required=False)
    parser.add_argument(
        'variants', type=str, help='variants truth set')
    return parser.parse_args()

def specificity(true_negatives, false_positives):
    return float(true_negatives) / (true_negatives + false_positives)

def sensitivity(true_positives, false_negatives):
    return float(true_positives) / (true_positives + false_negatives)

def compute_sens_spec(total_num_pathogenic, total_num_benign, num_pathogenic, num_benign):
    # num_pathogenic and num_benign are counts of things we've seen below the threshold
    # total_num_pathogenic and total_num_benign are total counts of ClinVar classified variants
    true_positives = total_num_pathogenic - num_pathogenic
    true_negatives = num_benign
    false_positives = total_num_benign - num_benign 
    false_negatives = num_pathogenic 
    sens = sensitivity(true_positives, false_negatives)
    spec = specificity(true_negatives, false_positives)
    return sens, spec

def sens_spec(data, min_score=0.0, max_score=1.0, resolution=100):
    '''Compute sensitivity and specificity for a range of cutoff values'''
    # assumes input data is sorted

    total_num_pathogenic = len([item for item in data if item[1] == 'pathogenic'])
    total_num_benign = len([item for item in data if item[1] == 'benign'])

    range = max_score - min_score
    increment = float(range) / resolution

    threshold = min_score
    results = []
    num_pathogenic = 0
    num_benign = 0
    current = min_score 

    for next_score, next_classification in data:
        if next_score >= threshold:
            sens, spec = compute_sens_spec(total_num_pathogenic, total_num_benign, num_pathogenic, num_benign)
            results.append((threshold, sens, spec)) 
            threshold += increment
        if next_classification == 'pathogenic':
            num_pathogenic += 1
        elif next_classification == 'benign':
            num_benign += 1

    sens, spec = compute_sens_spec(total_num_pathogenic, total_num_benign, num_pathogenic, num_benign)
    results.append((threshold, sens, spec)) 
    return results

def get_score(row, name, action):
    try:
        score = float(row[name])
    except:
        pass
    else:
        action(score)
    

def read_input(clinvar_filename):
    sift = []
    polyphen = []
    cadd = []
    pri_ph_cons = []
    gerp = []
    fathmm = []
    with open(clinvar_filename) as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            # get the ClinVar classification of this variant
            classification = row['benign']
            # get the score for each tool (where available)
            get_score(row, 'SIFT', lambda score: sift.append((1.0 - score, classification)))
            get_score(row, 'POLYPHEN', lambda score: polyphen.append((score, classification)))
            get_score(row, 'maxCADD', lambda score: cadd.append((score, classification)))
            get_score(row, 'priPhCons', lambda score: pri_ph_cons.append((score, classification)))
            get_score(row, 'GerpRS', lambda score: gerp.append((score, classification)))
            get_score(row, 'FATHMM', lambda score: fathmm.append((score, classification)))

    # sort in ascending order by score
    sift.sort()
    polyphen.sort()
    cadd.sort()
    pri_ph_cons.sort()
    gerp.sort()
    fathmm.sort()

    return sift, polyphen, cadd, pri_ph_cons, gerp, fathmm


class Plot(object):

    def __init__(self, dataset_name):
        self.dataset_name = dataset_name

    def roc(self, values, label):
        sens = [v[1] for v in values]
        spec = [1.0 - v[2] for v in values]
        title = "ROC"
        if self.dataset_name is not None:
            title += " (" + self.dataset_name + ")"
        plt.title(title)
        plt.plot(spec, sens, label=label)
    
    def plot_roc_curves(self, sift_results, polyphen_results, cadd_results, pri_ph_cons_results, gerp_results, fathmm_results):
        plt.xlabel("1-specificity")
        plt.ylabel("sensitivity")
    
        self.roc(sift_results, "sift")
        self.roc(polyphen_results, "polyphen")
        self.roc(cadd_results, "cadd")
        self.roc(pri_ph_cons_results, "pri_ph_cons")
        self.roc(gerp_results, "gerp")
        self.roc(fathmm_results, "fathmm")
        plt.legend(loc="lower right")
    
        plt.savefig("roc.png")
        plt.close()

    def plot_sens_spec_graph(self, algorithm, values, legend_loc='lower center', min_x=0.0, max_x=1.0):
        plt.xlabel("score threshold")
        thresholds = [v[0] for v in values]
        sens = [v[1] for v in values]
        spec = [v[2] for v in values]
        plt.plot(thresholds, sens, label='sensitivity: TP/(TP + FN)')
        plt.plot(thresholds, spec, label='specificity: TN/(TN + FP)')
        plt.legend(loc=legend_loc)
        title = algorithm + " sensitivity and specificity"
        if self.dataset_name is not None:
            title += " (" + self.dataset_name + ")"
        plt.title(title)
        plt.xlim([min_x, max_x]) 
        plt.ylim([0, 1]) 
        plt.savefig(algorithm + '.png')
        plt.close()


    

def main():
    args = parse_args()
    sift, polyphen, cadd, pri_ph_cons, gerp, fathmm = read_input(args.variants)

    max_cadd = cadd[-1][0]
    max_gerp = gerp[-1][0]

    sift_results = sens_spec(sift)
    polyphen_results = sens_spec(polyphen)
    cadd_results = sens_spec(cadd, min_score=0.0, max_score=max_cadd)
    pri_ph_cons_results = sens_spec(pri_ph_cons)
    gerp_results = sens_spec(gerp, max_score=max_gerp)
    fathmm_results = sens_spec(fathmm)

    plotter = Plot(args.dataset)

    plotter.plot_sens_spec_graph('sift', sift_results, legend_loc='center left')
    plotter.plot_sens_spec_graph('polyphen', polyphen_results)
    plotter.plot_sens_spec_graph('cadd', cadd_results, legend_loc='center right', max_x=max_cadd)
    plotter.plot_sens_spec_graph('pri_ph_cons', pri_ph_cons_results)
    plotter.plot_sens_spec_graph('gerp', gerp_results, legend_loc='center right', max_x=max_gerp)
    plotter.plot_sens_spec_graph('fathmm', fathmm_results)

    plotter.plot_roc_curves(sift_results, polyphen_results, cadd_results, pri_ph_cons_results, gerp_results, fathmm_results)

if __name__ == '__main__':
    main()
