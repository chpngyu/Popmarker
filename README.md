# Popmarker
Popmarker: identifying phylogenetic markers at the population level

As phylogenomic approach become a common practice for constructing a tree reflecting true bacterial phylogeny, it has become apparent that single molecular markers such as 16s rDNA often lead to misclassification of species. In this study, we present a program: Popmarker that utilize the true species phylogeny and identify a minimum set of proteins/genes reflecting the bacterial evolution history and phylogenetic relationship at the resolution of populations. 
  Popmarker rank the proteome according to the correlation of whole species tree or sub-tree branch length against orthologous sequence distances. We demonstrate five proteins of two top ranks achieve the same resolution as concatenation of 2,203 single copy orthologous genes and the right species classification as well as correct split of the two groups of Vibrio campbellii. The top ranking genes selected by Popmarker are candidates that lead to speciation and are useful in distinguishing close related species in microbiome study.

Popmarker is implemented by Python and there is no need to install it. This pipeline includes the following five steps, where the root of directory is assumed in $HOME, the source codes are in the $HOME/source, and the example data used in this pipeline are given in the $HOME/Example.

# Prerequisites Packages

  - Numpy (1.11.1 or higher)
  - Scipy (0.17.1 or higher)
  - Matplotlib (1.5.1 or higher)
  - Biopython (1.67 or higher)

and an additional package StatsModels (0.6.1 or later) is recommended for calculating the False Discovery Rate (FDR).

### Step1 Finding Orthoship

First, go to the folder of Step1
```sh
$ cd $HOME/Example/Step1/
```
and obtain 1-to-1 orthologous relationships among species (48 species in sp_list.txt) based on the output of OrthoFinder (OrthologousGroups.txt)
```sh
$ python $HOME/source/orthoships.py -o $HOME/Example/InputFile/OrthologousGroups.txt -s $HOME/ Example /InputFile/sp_list.txt -r orthoships_result.txt
```
The orthoships_result.txt will give 2287(/14265 in total) 1-to-1 orthologous groups. For the detailed options, -h (or --help) is available
```sh
$ python $HOME/source/orthoships.py -h
```

### Step2 Constructing Species Tree
Go to the folder of Step2
```sh
$ cd $HOME/Example/Step2.SpeciesTree/
```
Construct a reference tree using above all orthologous groups
```sh
$ python $HOME/source/maketree.py --accuracy localpair --gap 10 --wag --gamma --orthogrp $HOME/Example/Step1/orthoships_result.txt -s sp_file_path.txt 
```

The sp_file_path.txt tells maketree.py where the files of the orthologous sequences are located (see $HOME/Example/sp_file_path.txt). This script will produce concatenated sequences for all species (concatenated_seqs_2206.aln) and a constructed tree (concatenated_seqs_2206.tree). For the detailed options, type

```sh
$ python $HOME/source/maketree.py -h
```

### Step3. Ranking Gene (Rank1 and/or Rank2)
Go to the folder of Step3
```sh
$ cd $HOME/Example/Step3.RankGene/
```
Copy (symbolically link) the species tree (concatenated_seqs_2206.tree) and 1-to-1 orthologous groups (orthoship_result.txt) in the current folder by
```sh
$ cp -s $HOME/Example/Step2.SpeciesTree/concatenated_seqs_2206.tree ref.sp.tree
$ cp -s $HOME/Example/Step1/orthoships_result.txt .
```
and copy the result of all-against-all blastp for each pair of orthologous protein sequences among species by
```sh
$ cp $HOME/Example/goodProteins.blast .
```
It has an option for comparing species in only a clade as the following format (see $HOME/ Example/Vca_clade.txt)
Species_name1
Species_name2
Species_name3
…

To rank genes with respect to the reference tree or a clade of species if –clade is given as following:
```sh
$ python $HOME/source/rankgene.py -b goodProteins.blast -t  ref.sp.tree -o orthoships_result.txt --clade Vca_clade.txt --corr 0 -r rankgene_PCC_result.txt
```
This command uses Pearson correlation coefficient (--corr 0) as the measurement and also tell program to rank genes in accord with the Vca clade (--clade Vca_clade.txt). For the other options, type
```sh
$ python $HOME/source/rankgene.py -h
```

### Step4. Constructing Gene Tree
Go to the folder of Step4
```sh
$ cd $HOME/Example/Step4.GeneTree/
```
and copy three previously obtained files for later use
```sh
$ cp -s $HOME/Example/Step2.SpeciesTree/sp_file_path.txt .
$ cp -s $HOME/Example/Step3.RankGene/rankgene_PCC_result.txt .
$ cp -s $HOME/Example/Step3.RankGene/rankgene_Kendall_result.txt .
```

To construct a gene tree based on the top N genes of the ranking gene list
```sh
$ python $HOME/source/maketree.py --wag --gamma --accuracy localpair -g rankgene_PCC_result.txt -s sp_file_path.txt -t 1
```
This tells maketree.py to construct a gene tree using top one gene (-t 1) in the rank of rankgene_PCC_result.txt and to yield two outputs, ID000001.aln and ID000001.tree.
You might rename the output files to meaningful names as
```sh
$ mv ID000001.tree PCC.top1.tree
$ mv ID000001.aln PCC.top1.aln
```
Or, to construct a gene tree based on the concatenated multiple sequences in the top 2 genes (-t 2) as
```sh
$ python $HOME/source/maketree.py --wag --gamma --accuracy localpair -g rankgene_PCC_result.txt -s sp_file_path.txt -t 2 --con
```
It will produce concatenated_seqs_2.aln and concatenated_seqs_2.tree 
And you might rename the output files to
```sh
$ mv concatenated_seqs_2.tree PCC.top2.tree
$ mv concatenated_seqs_2.aln PCC.top2.aln
```

The maketree.py can construct the gene tree using the concatenated gene sequences from N top genes in the Species Rank and also can specify more top N-th genes in order to include top genes un the Clade Rank. For example, using Top4 genes (from 1st gene to 4th gene) and additional top 15th genes, which is the top 12th gene in the Clade Rank, in the Species Rank to construct a gene tree by
```sh
$ python $HOME/source/maketree.py --wag --gamma --accuracy localpair -g rankgene_PCC_result.txt -s sp_file_path.txt -t 4 -n 15 --con
```
It will produce concatenated_seqs_5.aln and concatenated_seqs_5.tree 
You might rename the output files to
```sh
$ mv concatenated_seqs_5.tree PCC.top4+15.tree
$ mv concatenated_seqs_5.aln PCC.top4+15.aln
```

The other example shows using the Top6 genes (1st gene to 6th gene) in the Species Rank and more top 8th gene, which is the 38th in Clade Rank, based on the correlation method of Kendall tau, as
```sh
$ python $HOME/source/maketree.py --wag --gamma --accuracy localpair -g rankgene_Kendall_result.txt -s sp_file_path.txt -t 6 -n 8 --con
```
It will produce concatenated_seqs_7.aln and concatenated_seqs_7.tree 
You might rename the output files to
```sh
$ mv concatenated_seqs_7.tree Kendall.top6+8.tree
$ mv concatenated_seqs_7.aln Kendall.top6+8.aln
```
For the detailed options, type
```sh
$ python $HOME/source/maketree.py -h
```

### Step5. Calculating Distance Between Two Trees
Go to the folder of Step5
```sh
$ cd $HOME/Example/Step5.TreeTreeDist/
```
Copy the species tree and the gene trees to the current folder for the comparisons:
```sh
$ cp –s $HOME/Example/Step2.SpeciesTree/concatenated_seqs_2206.tree ref.sp.tree
$ cp -s $HOME/Example/Step4.GeneTree/PCC.top4+15.tree .
$ cp -s $HOME/Example/Step4.GeneTree/PCC.top1.tree .
$ cp -s $HOME/Example/Step4.GeneTree/PCC.top2.tree .
$ cp -s $HOME/Example/Step4.GeneTree/PCC.top4+15.tree .
```
For calculating the distance between species tree and gene tree, it is done by
```sh
$ python $HOME/source/treedist.py -r ref.sp.tree -q *.tree -o whole_tree_distance.txt
```

For calculating the distance between only Vca clades in species tree and gene tree, it is
```sh
$ python $HOME/source/treedist.py -r ref.sp.tree -q *.tree -t Vca051011E Vca051011F Vca051011G Vca1114GL Vca1116 Vca151112c Vca200612B VcaCCS02 VcaDS40M4 VcaHY01 VcaKC13 VcaNBRC15631 VcaUMTGB204 -o VcaClade_tree_distance.txt
```
For the detailed options, type
```sh
$ python $HOME/source/treedist.py -h
```

[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)


   [dill]: <https://github.com/joemccann/dillinger>
   [git-repo-url]: <https://github.com/joemccann/dillinger.git>
   [john gruber]: <http://daringfireball.net>
   [df1]: <http://daringfireball.net/projects/markdown/>
   [markdown-it]: <https://github.com/markdown-it/markdown-it>
   [Ace Editor]: <http://ace.ajax.org>
   [node.js]: <http://nodejs.org>
   [Twitter Bootstrap]: <http://twitter.github.com/bootstrap/>
   [jQuery]: <http://jquery.com>
   [@tjholowaychuk]: <http://twitter.com/tjholowaychuk>
   [express]: <http://expressjs.com>
   [AngularJS]: <http://angularjs.org>
   [Gulp]: <http://gulpjs.com>

   [PlDb]: <https://github.com/joemccann/dillinger/tree/master/plugins/dropbox/README.md>
   [PlGh]: <https://github.com/joemccann/dillinger/tree/master/plugins/github/README.md>
   [PlGd]: <https://github.com/joemccann/dillinger/tree/master/plugins/googledrive/README.md>
   [PlOd]: <https://github.com/joemccann/dillinger/tree/master/plugins/onedrive/README.md>
   [PlMe]: <https://github.com/joemccann/dillinger/tree/master/plugins/medium/README.md>
   [PlGa]: <https://github.com/RahulHP/dillinger/blob/master/plugins/googleanalytics/README.md>
