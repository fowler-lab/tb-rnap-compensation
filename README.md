# tb-rnap-compensation
Investigation of the role of compensatory mutations in rifamycin resistance in M. tuberculosis

## Idea

I looked at this back in 2018 (!) but haven't been able to take it any further. The fundamental idea is that resistance to rifamipicin arises through mutations in the amino acids surrounding the drug which are located in the *rpoB* gene. These mutations have a fitness cost; this can be tolerated as long as the bacterium is exposed to antibiotics.

Then one of two things happen. Upon prolonged exposure, another mutation in the RNA polymerase arises that compensates for this fitness cost, thereby "locking in" the resistance-conferring mutation. Or, if exposure is removed (e.g. if treatment is stopped or the strain infects another untreated individual), the original mutation reverts. Most of the time this is silent i.e. the strain resembles wild-type again. Occasionally, however, it falls a different way on the codon table resulting in an apparent mutation at an amino acid where mutation is known to confer resistance but with 2 or 3 of the bases in the codon mutated. (e.g. S450W).

One can therefore look for compensatory mutations since they should only occur in conjunction with mutations known to confer resistance and never/rarely on their own. To date, it has been assumed that they all occur in *rpoC* but there is no good reason why they cannot occur in any protein in the RNAP complex, including *rpoB*.

One can also look for sites where reversion is probably happening by looking for mutations involving more than one base in the codon.

We now have at least twice as much data, including much more accurate and comprehensive phenotype, growth and lineage information, and there is a good story here, so I hope we can rapidly make progress and, if it goes well, think about writing a preprint.

## Way of working

The idea is this repository contains some information and data to help you get started. If you could e.g. create Jupyter notebooks within this repo and `git push` back to GitHub at the end of every day (or when you send me a question) I can have a look at your code.

When using a Jupyter notebook I'd strongly recommend including cells with Markdown explaining what you are thinking and doing, so that it becomes, in effect, an electronic lab notebook.

## Useful resources

* I've included the PDFs of a few papers in `references/` to get you started.
* A subset of the CRyPTIC datatables are in `tables/` and include a Data Schema. The dataset is hetergeneous, especially on the phenotype/growth side.
* [This website](https://mycobrowser.epfl.ch/) is useful for looking up genes in the *M. tuberculosis* H37Rv v3 genome.
