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
    parser.add_argument("--n_resistant", default=50, type=int, help="the minimum number of samples to allow in a resistant.")
    parser.add_argument("--n_other", default=50, type=int, help="the minimum number of samples to allow in an other.")
    parser.add_argument("--syn", default=True, type=bool, help="whether to include synonymous mutations in the analysis or not.")
    parser.add_argument("--debug", action='store_true', help="the minimum number of samples to allow in an other.")
    parser.add_argument("--outfile", default='results.csv', help="the minimum number of samples to allow in an other.")
    options = parser.parse_args()

    # don't let people give you negative numbers
    assert options.n_resistant > 0, "cannot have a negative number!"
    assert options.n_other > 0, "cannot have a negative number!"
    assert type(options.syn) == bool, "type must be a boolean!"

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

    if options.syn!=True:
        MUTATIONS=MUTATIONS[MUTATIONS.SYNONYMOUS==False]

    OTHER_MUTATIONS = pandas.DataFrame(list(set(MUTATIONS.GENE_MUTATION.unique()) - set(RESISTANT_MUTATIONS.GENE_MUTATION.unique())), columns=['GENE_MUTATION'])
    OTHER_MUTATIONS = OTHER_MUTATIONS[['GENE_MUTATION']]

    # now that we've defined both sets of mutations, let's subset down based on the minimum number of samples with each mutation
    df = EFFECTS.GENE_MUTATION.value_counts()>options.n_resistant
    df = pandas.DataFrame(df[df])
    df.reset_index(inplace=True)
    HIGH_COUNTS = df[['index']]
    HIGH_COUNTS.rename(columns={'index': 'GENE_MUTATION'}, inplace=True)
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
    HIGH_COUNTS = df[['index']]
    HIGH_COUNTS.rename(columns={'index': 'GENE_MUTATION'}, inplace=True)
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
        print("There are %i samples" %(len(SAMPLES)))

    rows=[]

    # iterate through the non-resistant mutations first since there are more of them
    for other_mutation in tqdm(OTHER_MUTATIONS.GENE_MUTATION):

        OTHER_SAMPLES = copy.deepcopy(MUTATIONS.loc[MUTATIONS.GENE_MUTATION == other_mutation])
        OTHER_SAMPLES['IS_OTHER'] = True
        OTHER_SAMPLES = OTHER_SAMPLES[['IS_OTHER']]

        for resistant_mutation in RESISTANT_MUTATIONS.GENE_MUTATION:

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
            fisher_test_set = pandas.crosstab(df.IS_RESISTANT, df.IS_OTHER)
            fisher_test_set = fisher_test_set.to_numpy()

            p = tb_rnap_compensation.calculate_fisher_pvalue(fisher_test_set)

            # should be a bit faster to build a list and convert to a DataFrame once finished
            rows.append([resistant_mutation, other_mutation, p.right_tail, p.left_tail])

            # consider using one-sided test, here resultsGreater -> means we expect a larger odds ratio
            # print(resistant_mutation, other_mutation, oddsr, p)

    # now convert back to a DataFrame and save to disc
    results = pandas.DataFrame(rows,columns=['resistant_mutation', 'other_mutation','p_right_tail','p_left_tail'])

    results.to_csv(options.outfile, index=False)
