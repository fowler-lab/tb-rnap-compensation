#! /usr/bin/env python

import argparse

import pandas
import xlsxwriter

import tb_rnap_compensation

if __name__ == "__main__":

    # collect the command line arguments using argparse
    parser = argparse.ArgumentParser()

    parser.add_argument("--results_path", required=False, default='/Users/viktoriabrunner/Documents/Studium/PhD/Project1_rev/repository/tb-rnap-compensation/', help="the path to the folder containing the results.csv file from calculate-fisher-tests.py and the reference file.")
    parser.add_argument("--p_value", default=0.01, type=float, help="the cut-off for significantly associated mutations.")
    parser.add_argument("--method", default='inclusive', type=str, help="the type of reference list to be used for comparison. Either 'inclusive' or 'conservative'.")
    parser.add_argument("--debug", action='store_true', help="print info on running analysis to terminal.")
    options = parser.parse_args()

    # don't let people give you negative numbers
    assert options.p_value > 0, "cannot have a negative number!"

    # load csv with results from fishers exact test
    results=pandas.read_csv(options.results_path + 'results.csv')

    # set parameters
    p_value = options.p_value/len(results)
    n_tests = len(results)
    method = options.method
    if method == 'inclusive':
        sheet_name='described_CMs_binary'
    if method == 'conservative':
        sheet_name='conservative_CMs'

    if options.debug:
        print("There were %i tests performed and the analysis looks for other mutations that are significantly correlated to resistance mutations on a %i percent level." %(n_tests, options.p_value*100))

    # load known CMs from reference papers
    reference=pandas.read_excel(options.results_path + 'Ref_known_CMs.xlsx', sheet_name=sheet_name)
    reference.drop([0,1,2,3], axis=0, inplace=True)
    reference.drop(index=reference.index[-2:], axis=0, inplace=True)

    # determine which of the reference mutations appear in ALL tested other_mutations
    unique_results=pandas.Series(results['other_mutation'].unique())

    # determine which other_mutations are significantly correlated with resistance mutations on a 1% level including bonferroni multiple testing correction
    hits = results[(results.p_right_tail<(p_value))].other_mutation
    unique_hits = pandas.Series(hits.unique())

    if options.debug:
        print("There are %i other mutations significantly correlated on a %i percent level" %(len(unique_hits), options.p_value*100))

    # determine which percentage of reference CMs are significantly correlated in our dataset
    found_ref_sign = len(unique_hits[unique_hits.isin(reference['mutation'])])/len(unique_results[unique_results.isin(reference['mutation'])])

    # determine which percentage of reference CMs are found in our significant hits
    found_ref = len(unique_hits[unique_hits.isin(reference['mutation'])])/len(reference)

    # store reference hits that are found with their corresponding resistance mutations and p-values in dataframe
    ref_hits = unique_hits[unique_hits.isin(reference['mutation'])]
    ref_hits = results[results['other_mutation'].isin(ref_hits)&(results.p_right_tail<(p_value/n_tests))]

    # store previoulsy not described hits with their corresponding resistance mutation and p-vlaues in dataframe
    new_hits = unique_hits[~unique_hits.isin(reference['mutation'])]
    new_hits = results[results['other_mutation'].isin(new_hits)&(results.p_right_tail<(p_value/n_tests))]

    # write xlsx file with metadata of analysis and significantly correlated reference and new hits
    rows = []
    rows.append([p_value, n_tests, method, found_ref, found_ref_sign])
    analysis = pandas.DataFrame(rows, columns = ['p_value', 'n_tests', 'method', 'found_ref', 'sign_found_ref'])

    rows1 = []
    rows1 = ref_hits

    rows2 = []
    rows2 = new_hits

    writer = pandas.ExcelWriter('Fisher_hits_analysis.xlsx', engine='xlsxwriter')
    analysis.to_excel(writer, sheet_name='parameters_statistics', index=False)
    rows1.to_excel(writer, sheet_name='reference_hits', index=False)
    rows2.to_excel(writer, sheet_name='new_hits', index=False)
    writer.save()   
