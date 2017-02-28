# Gondi-Dialect-Analysis

* This paper is about analyzing 46 Gondi dialects dialects spoken in Central India.

* The paper is available as [pdf] (paper/gondi_dialects.pdf) here and is accepted for publishing at [VarDial 2017](http://ttg.uni-saarland.de/vardial2017/index.html)

> Please cite the paper as:
```
@inproceedings{rama-coltekin-sofroniev2017gondi,
  title = {A computational analysis fo Gondi dialects},
  author = {Rama, Taraka and \c{C}\"{o}ltekin, \c{C}a\u{g}r\i{} and Sofroniev, Pavel},
  date = {2017},
  booktitle = {Proceedings of the Fourth Workshop on {NLP} for Similar Languages, Varieties and Dialects},
  pages = {(to appear)},
  location = {Valencia, Spain},
}
```
# Cognate clustering and MrBayes

The program online_pmi.py produces the nexus file which can be fed into MrBayes for the purpose of tree building.

The program takes three inputs: a file with a seed list of probable cognates, a data file containing the data and the coding option for processing the file. The program then processes the file by computing the PMI scores for sound segments and then uses the scoring matrix to cluster the words.

For example: `python3 online_pmi.py gondi_ldn_ipa1.txt data/gondi_combined.tsv IPA outputfilename`
The program outputs the cognate judgments for each word belonging to a concept and the number of clusters found for each concept.


The data folder contains `gondi_combined.tsv` file that contains the word lists in IPA, ASJP, and SCA format.
Another file `gondi_combined_cognates.csv` contains the cognate information given by Taraka Rama (the lead author of the paper).

Finally the `mrbayes/run.mb` file consists of [MrBayes](http://mrbayes.sourceforge.net/) commands that produce the consensus tree that can be visualized using [FigTree](http://beast.bio.ed.ac.uk/figtree). We provide a consensus tree for visualization. The tree is a rooted tree and uses Independent Gamma Branch Rates for the purpose of inferring trees. The `.nexus` file is also provided in the `mrbayes` folder. The `.tre` file provides the consensus tree from our analysis.

## Maps
The `gondi.kml` file is in the `maps` folder that is useful for the purpose of visualization.
- [] Add a pdf to the map

## Autoencoders

## Requirements
