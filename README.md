# Comprehensive Guide on generating Admixture Graphs using Admixtools2 in R
Admixture graphs allow researchers to model historical admixture events, which is particularly useful in hybridization studies to calculate the relative contributions of parental lineages to hybrids. They provide additional insights beyond classic Treemix and f-branch plots. After an extensive search for a straightforward tutorial on generating admixture graphs from genetic variant data yielded little success, I decided to write a comprehensive guide to help integrate this powerful analytical tool into your research.
Find the full tutorial in https://www.timlheller/admixturegraphs


## Input Format

Admixtools requires an atypical bouquet of file formats to run, referred to as EIGENSTRAT, which must be converted accordingly: 
- A .ind file which is a tab-separated table without header containing the columns  ``` SampleID  Sex  Population ``` . For Sex, we can write either M (male) or F (female) if this information is accessible, otherwise we will declare it as U (undefined).
- A .geno file with all SNP variants and their relative states in the samples
- A .snp highlighting the positions of all SNP saved in the .geno files.

Whilst the conversion of VCFs to this assortment is complex, there exist scripts to do the work for you. We just need to provide a VCF, which should be LD pruned for unbiased results, and a population file inclosing two columns one: Sample	Population
The associated Linux package Admixtools2 provides a function, convertf, for conversion of Plink to EIGENSTRAT. However, at the time of writing this tutorial, I was unsuccessful with this application.

## Reading Eigenstrat and calculating f2

We now provide the file name stem provided to the conversion script as a file prefix and generate our f2 statistics. I will use the simulated data I provided in this repository, you can follow the tutorial.

``` library(admixtools)
prefix = "heller_sim"
f2_blocks = f2_from_geno(prefix)
``` 

## Finding the optimal tree topology 

Afterwards, we have two options to find the optimal trees to represent our data in an ideal way. We can either provide a tree and fit the gene flow accordingly or use this tree as a startpoint topology for the find_graph() function. The generation of the initial tree can be omitted, when the topology is unknown.
The first argument is the previously calculated f2 statistics as an object, then the outgroup as specified in the population file is declared, then we define the number of admixture (hybridization) events and get to choose further iteration parameters. Here, we can specify a tree with initgraph = as an R-object 
```
candidate_graphs = find_graphs(
  f2_blocks, outpop="Outgroup", #numadmix = 1
  initgraph = tree_hybrid,
  stop_gen2 = 5,
  numgraphs = 10,
)
 ```

## Fitting data to tree topology
To fit the data now, we run the command qpgraph and specify the f2 statistics as the first and the tree topology as the second argument. The column $score displays the model fit for the data on the topology. We can plot these admixture and drift values using the $edges column. 

If you have all elements within the candidate_graphs, we can automate the fitting and evaluation of every graph as follows: 

``` library(tidyverse)
best_graph = NULL
best_score = Inf
for (graph in candidate_graphs$graph) {
  fit_find = qpgraph(f2_blocks, graph)
  if (fit_find$score < best_score) {
    best_score = fit_find$score
    best_run = graph
    best_tree = fit_find$edges    
  }
}
```
## Bootstrapping your data 
We are getting very close to the final product, the only step left is bootstrapping by applying variations to your dataset. We use the following command. The second element is a list with all fitted graphs as R-objects and the third element denounces the number of bootstrap iteration. We can also simply put one graph into the list() function.
```
bootstrap = qpgraph_resample_multi(f2_blocks, list(best_tree,fit_graph), nboot = 100)
compare_fits(fits[[1]]$score, fits[[2]]$score)
bootstrap[[2]] %>% summarize_fits() %>% plotly_graph(print_highlow = TRUE)
```
If you encounter complications even after thoroughly reading this article, feel free to contact me.
