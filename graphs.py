#!/usr/bin/env python

'''
Compute and plot sensitivity and specificity for various
variant effect prediction tools.

Usage: python graphs.py clinvar.summary.tsv

Note: Sift scores range from 1 to 0, but we convert the range
to 0-1 to make it easier to plot. We need to relabel the X axis
'''

import csv
import sys
import matplotlib.pyplot as plt

def specificity(true_negatives, false_positives):
    return float(true_negatives) / (true_negatives + false_positives)

def sensitivity(true_positives, false_negatives):
    return float(true_positives) / (true_positives + false_negatives)

def compute_sens_spec(total_num_pathogenic, total_num_benign, num_pathogenic, num_benign):
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

    threshold = min_score + increment 
    results = []
    num_pathogenic = 0
    num_benign = 0
    current = min_score 

    for next_score, next_classification in data:
        if next_score > threshold:
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
            try:
                sift_score = 1.0 - float(row['SIFT'])
                polyphen_score = float(row['POLYPHEN'])
                cadd_score = float(row['maxCADD'])
                pri_ph_cons_score = float(row['priPhCons'])
                gerp_score = float(row['GerpRS'])
                fathmm_score = float(row['FATHMM'])
                classification = row['benign']
            except:
                # ignore any row where we are missing scores 
                pass
            else:
                sift.append((sift_score, classification))
                polyphen.append((polyphen_score, classification))
                cadd.append((cadd_score, classification))
                pri_ph_cons.append((pri_ph_cons_score, classification))
                gerp.append((gerp_score, classification))
                fathmm.append((fathmm_score, classification))

    # sort by score
    sift.sort()
    polyphen.sort()
    cadd.sort()
    pri_ph_cons.sort()
    gerp.sort()
    fathmm.sort()

    return sift, polyphen, cadd, pri_ph_cons, gerp, fathmm

def plot_graph(algorithm, values, legend_loc='lower center', min_x=0.0, max_x=1.0):
    plt.xlabel("Score Threshold")
    plt.title(algorithm)
    thresholds = [v[0] for v in values]
    sens = [v[1] for v in values]
    spec = [v[2] for v in values]
    plt.plot(thresholds, sens, label='sensitivity: TP/(TP + FN)')
    plt.plot(thresholds, spec, label='specificity: TN/(TN + FP)')
    plt.legend(loc=legend_loc)
    plt.xlim([min_x, max_x]) 
    plt.ylim([0,1]) 
    plt.savefig(algorithm + '.png')
    plt.close()

def main():
    if len(sys.argv) < 2:
        exit("usage: graphs.py <clinvar_stats.tsv>")
    sift, polyphen, cadd, pri_ph_cons, gerp, fathmm = read_input(sys.argv[1])

    max_cadd = cadd[-1][0]
    max_gerp = gerp[-1][0]

    sift_results = sens_spec(sift)
    polyphen_results = sens_spec(polyphen)
    cadd_results = sens_spec(cadd, min_score=0.0, max_score=max_cadd)
    pri_ph_cons_results = sens_spec(pri_ph_cons)
    gerp_results = sens_spec(gerp, max_score=max_gerp)
    fathmm_results = sens_spec(fathmm)

    plot_graph('sift', sift_results, legend_loc='center left')
    plot_graph('polyphen', polyphen_results)
    plot_graph('cadd', cadd_results, legend_loc='center right', max_x=max_cadd)
    plot_graph('pri_ph_cons', pri_ph_cons_results)
    plot_graph('gerp', gerp_results, legend_loc='center right', max_x=max_gerp)
    plot_graph('fathmm', fathmm_results)

if __name__ == '__main__':
    main()
