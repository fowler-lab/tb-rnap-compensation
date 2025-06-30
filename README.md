[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](
https://colab.research.google.com/github/fowler-lab/tb-rnap-compensation/blob/main/Recreate_figures_paper.ipynb
)
# tb-rnap-compensation
Brunner V, Fowler PW

Compensatory mutations are associated with increased in vitro growth in resistant clinical samples of *Mycobacterium tuberculosis*.

mGen 10:001187 https://doi.org/10.1099/mgen.0.001187

## Idea

Rifampicin resistance is introduced by mutations in the RNA polymerase of *M. tuberculosis*. This resistance also introduces a fitness cost in M. tuberculosis bacteria. This fitness cost can be alleviated by compensatory mutations (CMs) in other parts of the same protein. 
The direct causal link between each resistance mutation and compensatory mutation should allow identifying a comprehensive set of compensatory mutations through statistical association testing, if conducted with a large enough dataset.

We have around 70,000 sequenced *M. tuberculosis* genomes on hand, which can be employed to identify a high-confidence set of CMs. Pairing this with *in vitro* growth data (available for about 15,000 samples) allows evaluation of the genotype-phenotype relationship. 

## Way of working

### Recreate figures

This repository is the basis for the paper *"Compensatory Mutations are associated with increased in vitro Growth Levels in resistant Clinical Samples of Mycobacterium tuberculosis"* Clicking "Open In Colab" above the readme will take you to a google colab that you can use to run the jupyter notebook `Recreate_figures_paper.ipynb`. This notebook allows to reproduce most of the figures as seen in the paper. A few exceptions are figures based on PyMOL mappings, which were produced outside of the notebook environment. 

Important notes:
- The first cell of the notebook is only to be executed when using the google colab environment. It mounts the repository on the instance of google colab that is running.
- the medians and confidence intervals of the boxplots are slightly different from the ones shown in the paper figures, since `seaborn` stopped support of the `usermedians` command. The medians are hence not bootstrapped but calculated using seaborns inbuilt function for boxplot medians.

### Run statistical association tests

The statistical test used (fisher test for association of resistance mutations with potential CMs) can be reproduced as well by using cloning the github repository to your local machine and using `setup.py` to install the module. The command `calculate-fisher-tests.py` should then allow to run pairwise testing of resistance mutations with other mutations detected in the RNA polymerase. It is possible to specify several parameters.

There is also the option to run `results-evaluation.py`, which will automatically take the list of hits resulting from `calculate-fisher-tests.py` (as long as the output file is called 'results.csv') and apply the specified p-value cutoff with/ without multiple testing correction. There is also the option to automatically exclude synonymous mutations, flag mutations known to be lineage associated, and to compare the resulting hits to known literature CMs.

## Useful resources

* Several important references are included `references/`
* The subset of the CRyPTIC datatables used in this project are in `tables/` and include a Data Schema. The dataset is hetergeneous, especially on the phenotype/growth side.
* [This website](https://mycobrowser.epfl.ch/) is useful for looking up genes in the *M. tuberculosis* H37Rv v3 genome.
