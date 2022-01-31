import argparse
import copy

import pandas
from fisher import pvalue

from tqdm import tqdm

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--table_path", required=True, help="the path to the folder containing all the CRyPTIC tables.")
    parser.add_argument("--n_resistant", default=50, type=int, help="the minimum number of samples to allow in a resistant.")
    parser.add_argument("--n_other", default=50, type=int, help="the minimum number of samples to allow in an other.")
    options = parser.parse_args()

    assert options.n_resistant > 0, "cannot have a negative number!"

    MUTATIONS = pandas.read_pickle(options.table_path + '/MUTATIONS.pkl.gz')
    MUTATIONS.reset_index(inplace=True)
    MUTATIONS = MUTATIONS[(MUTATIONS.IS_FILTER_PASS) & (~MUTATIONS.IS_HET) & (~MUTATIONS.IS_NULL)]
    MUTATIONS.set_index('UNIQUEID', inplace=True)
    MUTATIONS = MUTATIONS.astype({'GENE':'str'})
    MUTATIONS['GENE_MUTATION'] = MUTATIONS['GENE'] + '_' + MUTATIONS['MUTATION']

    EFFECTS = pandas.read_pickle(options.table_path + '/EFFECTS.pkl.gz')
    EFFECTS.reset_index(inplace=True)
    EFFECTS['GENE_MUTATION'] = EFFECTS['GENE'] + '_' + EFFECTS['MUTATION']
    EFFECTS = EFFECTS[(EFFECTS.DRUG=='RIF') & (EFFECTS.PREDICTION=='R') & (~EFFECTS.MUTATION.str[-1].isin(['O','X']))]
    EFFECTS.set_index('UNIQUEID', inplace=True)

    # create a table with one row per sample that is resistant to RIF
    RESISTANT_SAMPLES = copy.deepcopy(EFFECTS)
    RESISTANT_SAMPLES['IS_RESISTANT'] = True
    RESISTANT_SAMPLES=RESISTANT_SAMPLES[['IS_RESISTANT']]

    RESISTANT_MUTATIONS = pandas.DataFrame(EFFECTS.GENE_MUTATION.unique(), columns = ['GENE_MUTATION'])
    RESISTANT_MUTATIONS['IS_RESISTANT'] = True


    OTHER_MUTATIONS = pandas.DataFrame(list(set(MUTATIONS.GENE_MUTATION.unique()) - set(RESISTANT_MUTATIONS.GENE_MUTATION.unique())), columns=['GENE_MUTATION'])
    OTHER_MUTATIONS['OTHER_MUTATIONS'] = True
    OTHER_MUTATIONS = OTHER_MUTATIONS[['GENE_MUTATION']]

    # Remove all mutations from OTHER_MUTATIONS df that appear less than n times in all samples
    a = MUTATIONS.GENE_MUTATION.value_counts()>options.n_other
    df = pandas.DataFrame(a[a])
    df.reset_index(inplace=True)
    HIGH_COUNTS = df[['index']]
    HIGH_COUNTS.rename(columns={'index': 'GENE_MUTATION'}, inplace=True)
    HIGH_COUNTS.set_index('GENE_MUTATION', inplace=True)

    OTHER_MUTATIONS.set_index('GENE_MUTATION', inplace=True)
    print(len(OTHER_MUTATIONS))
    OTHER_MUTATIONS = OTHER_MUTATIONS.join(HIGH_COUNTS, how='inner')
    # OTHER_MUTATIONS = OTHER_MUTATIONS.loc[OTHER_MUTATIONS.GENE_MUTATION, :]
    # OTHER_MUTATIONS = OTHER_MUTATIONS.drop(columns=['GENE_MUTATION'])
    OTHER_MUTATIONS.reset_index(inplace=True)
    # print(OTHER_MUTATIONS[:3])

    print(len(OTHER_MUTATIONS))

    a = EFFECTS.GENE_MUTATION.value_counts()>options.n_resistant
    df = pandas.DataFrame(a[a])
    df.reset_index(inplace=True)
    HIGH_COUNTS = df[['index']]
    HIGH_COUNTS.rename(columns={'index': 'GENE_MUTATION'}, inplace=True)

    RESISTANT_MUTATIONS = RESISTANT_MUTATIONS.loc[RESISTANT_MUTATIONS.GENE_MUTATION.isin(HIGH_COUNTS.GENE_MUTATION)]


    SAMPLES = pandas.DataFrame(MUTATIONS.index.unique(), columns=['UNIQUEID'])
    SAMPLES.set_index('UNIQUEID', inplace=True)

    results = pandas.DataFrame(columns=('resistant_mutation', 'other_mutation','p_right_tail','p_left_tail'))
    counter = 0

    for resistant_mutation in tqdm(RESISTANT_MUTATIONS.GENE_MUTATION):

        RESISTANT_SAMPLES = copy.deepcopy(EFFECTS.loc[EFFECTS.GENE_MUTATION == resistant_mutation])
        RESISTANT_SAMPLES['IS_RESISTANT'] = True
        RESISTANT_SAMPLES=RESISTANT_SAMPLES[['IS_RESISTANT']]

        for other_mutation in OTHER_MUTATIONS.GENE_MUTATION:
            OTHER_SAMPLES = copy.deepcopy(MUTATIONS.loc[MUTATIONS.GENE_MUTATION == other_mutation])
            OTHER_SAMPLES['IS_OTHER'] = True
            OTHER_SAMPLES = OTHER_SAMPLES[['IS_OTHER']]

            df = copy.deepcopy(SAMPLES)
            df = df.join(RESISTANT_SAMPLES, how='left')
            # print(df)
            df = df.join(OTHER_SAMPLES, how='left')
            df.fillna(False,inplace=True)
            # print(df[df.IS_RESISTANT==True])
            fisher_test_set = pandas.crosstab(df.IS_RESISTANT, df.IS_OTHER)
            fisher_test_set = fisher_test_set.to_numpy()
            p = pvalue(fisher_test_set[0,0],fisher_test_set[0,1],fisher_test_set[1,0],fisher_test_set[1,1])
            # print(fisher_test_set)
            # oddsr, p = fisher_exact(fisher_test_set, alternative='two-sided')

            results.loc[counter]=[resistant_mutation, other_mutation, p.right_tail, p.left_tail]
            counter = counter+1
            #consider using one-sided test, here resultsGreater -> means we expect a larger odds ratio
            # print(resistant_mutation, other_mutation, oddsr, p)

    print(results)
