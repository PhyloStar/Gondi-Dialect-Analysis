# Gondi-Dialect-Analysis

This paper is about analyzing 46 Gondi dialects dialects spoken in Central India.

# Cognate clustering and MrBayes

The program online_pmi.py produces the nexus file which can be fed into MrBayes for the purpose of tree building.

The program takes three inputs: a file with a seed list of probable cognates, a data file containing the data and the coding option for processing the file. The program then processes the file by computing the PMI scores for sound segments and then uses the scoring matrix to cluster the words.

For example: `python3 online_pmi.py gondi_ldn_ipa1.txt data/gondi_combined.tsv IPA <outputfilename>`

