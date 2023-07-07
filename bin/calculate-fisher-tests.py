#! /usr/bin/env python

import argparse
import copy
import pkg_resources

import pandas
from fisher import pvalue
from tqdm.autonotebook import tqdm

import tb_rnap_compensation

if __name__ == "__main__":

    # collect the command line arguments using argparse
    parser = argparse.ArgumentParser()

    parser.add_argument("--table_path", required=False, default=pkg_resources.resource_filename("tb_rnap_compensation", 'tables/'), help="the path to the folder containing all the CRyPTIC tables.")
    parser.add_argument("--test_method", default='numerical', help="the statistical association test used. Options are 'numerical' or 'fisher' test.")
    parser.add_argument("--n_resistant", default=50, type=int, help="the minimum number of samples to allow in a resistant.")
    parser.add_argument("--n_other", default=50, type=int, help="the minimum number of samples to allow in an other.")
    parser.add_argument("--include_syn", action='store_true', help="whether to include synonymous mutations in the analysis or not.")
    parser.add_argument("--debug", action='store_true', help="the minimum number of samples to allow in an other.")
    parser.add_argument("--outfile", default='results.csv', help="the minimum number of samples to allow in an other.")
    options = parser.parse_args()

    # don't let people give you negative numbers
    assert options.n_resistant > 0, "cannot have a negative number!"
    assert options.n_other > 0, "cannot have a negative number!"

    # read in the MUTATIONS data table
    MUTATIONS = pandas.read_pickle(options.table_path + '/MUTATIONS.pkl.gz')
    MUTATIONS.reset_index(inplace=True)
    MUTATIONS = MUTATIONS.loc[(MUTATIONS.IS_FILTER_PASS) & (~MUTATIONS.IS_HET) & (~MUTATIONS.IS_NULL)]
    MUTATIONS.set_index('UNIQUEID', inplace=True)

    #add column indicating synonymous mutations
    boolean = MUTATIONS['MUTATION'].apply(lambda x: x[0]==x[-1])
    MUTATIONS['SYNONYMOUS'] = boolean

    # convert GENE from categorical
    MUTATIONS = MUTATIONS.astype({'GENE':'str'})
    MUTATIONS['GENE_MUTATION'] = MUTATIONS['GENE'] + '_' + MUTATIONS['MUTATION']

    # read in the EFFECTS data table
    EFFECTS = pandas.read_pickle(options.table_path + '/EFFECTS.pkl.gz')
    EFFECTS.reset_index(inplace=True)
    EFFECTS['GENE_MUTATION'] = EFFECTS['GENE'] + '_' + EFFECTS['MUTATION']
    EFFECTS = EFFECTS.loc[(EFFECTS.DRUG=='RIF') & (EFFECTS.PREDICTION=='R') & (~EFFECTS.MUTATION.str[-1].isin(['O','X']))]
    EFFECTS.set_index('UNIQUEID', inplace=True)

    RESISTANT_MUTATIONS = pandas.DataFrame(EFFECTS.GENE_MUTATION.unique(), columns = ['GENE_MUTATION'])

    if not options.include_syn:
        MUTATIONS=MUTATIONS[~MUTATIONS.SYNONYMOUS]

    OTHER_MUTATIONS = pandas.DataFrame(list(set(MUTATIONS.GENE_MUTATION.unique()) - set(RESISTANT_MUTATIONS.GENE_MUTATION.unique())), columns=['GENE_MUTATION'])
    OTHER_MUTATIONS = OTHER_MUTATIONS[['GENE_MUTATION']]

    # now that we've defined both sets of mutations, let's subset down based on the minimum number of samples with each mutation
    df = EFFECTS.GENE_MUTATION.value_counts()>options.n_resistant
    df = pandas.DataFrame(df[df])
    df.reset_index(inplace=True)
    # HIGH_COUNTS = df[['index']]
    # HIGH_COUNTS.rename(columns={'index': 'GENE_MUTATION'}, inplace=True)
    HIGH_COUNTS = df[['GENE_MUTATION']]
    HIGH_COUNTS.set_index('GENE_MUTATION', inplace=True)
    RESISTANT_MUTATIONS.set_index('GENE_MUTATION', inplace=True)
    RESISTANT_MUTATIONS = RESISTANT_MUTATIONS.join(HIGH_COUNTS, how='inner')
    RESISTANT_MUTATIONS.reset_index(inplace=True)

    if options.debug:
        print("There are %i resistant mutations using a threshold of %i samples" %(len(RESISTANT_MUTATIONS), options.n_resistant))

    # Remove all mutations from OTHER_MUTATIONS df that appear less than n times in all samples
    df = MUTATIONS.GENE_MUTATION.value_counts()>options.n_other
    df = pandas.DataFrame(df[df])
    df.reset_index(inplace=True)
    # HIGH_COUNTS = df[['index']]
    # HIGH_COUNTS.rename(columns={'index': 'GENE_MUTATION'}, inplace=True)
    HIGH_COUNTS = df[['GENE_MUTATION']]
    HIGH_COUNTS.set_index('GENE_MUTATION', inplace=True)
    OTHER_MUTATIONS.set_index('GENE_MUTATION', inplace=True)
    OTHER_MUTATIONS = OTHER_MUTATIONS.join(HIGH_COUNTS, how='inner')
    OTHER_MUTATIONS.reset_index(inplace=True)

    if options.debug:
        print("There are %i non-resistant mutations using a threshold of %i samples" %(len(OTHER_MUTATIONS), options.n_other))

    # finally create a table of SAMPLES
    SAMPLES = pandas.DataFrame(MUTATIONS.index.unique(), columns=['UNIQUEID'])
    SAMPLES.set_index('UNIQUEID', inplace=True)

    if options.debug:
        print("There are %i samples with mutations in RNAP" %(len(SAMPLES)))

    rows=[]
    n_tests = 0

    # iterate through the non-resistant mutations first since there are more of them
    for other_mutation in tqdm(OTHER_MUTATIONS.GENE_MUTATION):

        OTHER_SAMPLES = copy.deepcopy(MUTATIONS.loc[MUTATIONS.GENE_MUTATION == other_mutation])
        OTHER_SAMPLES['IS_OTHER'] = True
        OTHER_SAMPLES = OTHER_SAMPLES[['IS_OTHER']]

        for resistant_mutation in tqdm(RESISTANT_MUTATIONS.GENE_MUTATION):

            RESISTANT_SAMPLES = copy.deepcopy(EFFECTS.loc[EFFECTS.GENE_MUTATION == resistant_mutation])
            RESISTANT_SAMPLES['IS_RESISTANT'] = True
            RESISTANT_SAMPLES=RESISTANT_SAMPLES[['IS_RESISTANT']]

            # create a copy and join to both the lists
            df = copy.deepcopy(SAMPLES)
            df = df.join(RESISTANT_SAMPLES, how='left')
            df = df.join(OTHER_SAMPLES, how='left')

            # since matching rows will have a True, and non-matching will have a NaN, fill the NaNs with False
            df.fillna(False,inplace=True)

            # now we can do the crosstab
            test_set = pandas.crosstab(df.IS_RESISTANT, df.IS_OTHER)
            test_set = test_set.to_numpy()

            #now we use one of the two methods to calculate the p-value
            if options.test_method == 'numerical':
                
                if test_set[1,1] > 0:
                    p =  tb_rnap_compensation.numerical_test(100000, 10000, len(RESISTANT_SAMPLES), len(OTHER_SAMPLES), test_set[1,1])
                    n_tests = n_tests + 1

                else:
                    p = 1

                rows.append([resistant_mutation, other_mutation, p, test_set[0,0], test_set[0,1], test_set[1,0], test_set[1,1], len(RESISTANT_SAMPLES), len(OTHER_SAMPLES)])

            if options.test_method == 'fisher':
                
                if test_set[1,1] > 0:
                    p = tb_rnap_compensation.calculate_fisher_pvalue(test_set)
                    n_tests = n_tests + 1

                    rows.append([resistant_mutation, other_mutation, p.right_tail, test_set[0,0], test_set[0,1], test_set[1,0], test_set[1,1], len(RESISTANT_SAMPLES), len(OTHER_SAMPLES)])


                else:
                    p = 1

                    rows.append([resistant_mutation, other_mutation, p, test_set[0,0], test_set[0,1], test_set[1,0], test_set[1,1], len(RESISTANT_SAMPLES), len(OTHER_SAMPLES)])

            if options.test_method == 'chi-square':

                if test_set[1,1] > 0:
                    stat, p, dof, expected = tb_rnap_compensation.calculate_chi_square_pvalue(test_set)
                    n_tests = n_tests + 1

                else:
                    stat, p, dof, expected = 1, 1, 1, 1

                rows.append([resistant_mutation, other_mutation, stat, p, dof, expected, test_set[0,0], test_set[0,1], test_set[1,0], test_set[1,1], len(RESISTANT_SAMPLES), len(OTHER_SAMPLES)])
    
    if options.test_method == 'chi-square':
        rows.append(['number', 'of tests', 'performed:', n_tests, 5, 6, 7, 8, 9, 10, 11, 12])
        # now convert back to a DataFrame and save to disc
        results = pandas.DataFrame(rows,columns=['resistant_mutation', 'other_mutation','chi-square statistic','p_value','dof','expected','None','other','resistant','both', 'n_resistant', 'n_other'])

    else:
        rows.append(['number', 'of tests', 'performed:', n_tests, 5, 6, 7, 8, 9])
        # now convert back to a DataFrame and save to disc
        results = pandas.DataFrame(rows,columns=['resistant_mutation', 'other_mutation','p_value','None','other','resistant','both', 'n_resistant', 'n_other'])

    results.to_csv(options.outfile, index=False)
