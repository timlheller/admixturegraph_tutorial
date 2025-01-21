#Loading relevant libraries
library(admixtools)
library(plotly)
library(tidyverse)

#Loading genotype data, ensure to be in the correct working directory
prefix = "heller_sim"
f2_blocks = f2_from_geno(prefix)

#OPTIONAL: Programming tree topologies
tree_simple = edges_to_igraph(matrix(c(
  'R','Outgroup',
  'R','n0',
  'n0','Taxon1',
  'n0', 'Taxon2'
  ), , 2, byrow = T))

plotly_graph(tree_simple, fix=T)

tree_hybrid = edges_to_igraph(matrix(c(
  'R','Outgroup',
  'R','n0',
  'n0','n1',
  'n0','n2',
  'n1','Taxon1',
  'n1','a1',
  'n2','Taxon2',
  'n2','a1',
  'a1','Hybrid'
  ), , 2, byrow = T)) 
plotly_graph(tree_hybrid,fix=T)

#Running the graph finding algorithm
candidate_graphs = find_graphs(
  f2_blocks, outpop="Outgroup", #numadmix=1, 
  initgraph = tree_hybrid,
  stop_gen2 = 5,
  numgraphs = 10,
)

#OPTIONAL: Fitting the data to a predefined topology 
fit_function = qpgraph(f2_blocks,tree_hybrid)
fit_graph = fit_function$edges
plotly_graph(fit_graph,fix=T)
fit_function$score

#Iterating through the generated tree topologies
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
best_tree %>% plotly_graph(fix=T)

#Bootstrapping the admixture and drift parameters
bootstrap = qpgraph_resample_multi(f2_blocks, list(best_tree,fit_graph), nboot = 100)
compare_fits(bootstrap[[1]]$score, bootstrap[[2]]$score)
bootstrap[[2]] %>% summarize_fits() %>% plotly_graph(print_highlow = TRUE, fix=T)

#OPTIONAL: Exporting the graphs to SVG. The resulting file will be spawned in the working directory. 
#You need python to be installed and configured for this to work.
install.packages("reticulate")
library(reticulate)
reticulate::py_config()
reticulate::py_install("kaleido")

final_plot = bootstrap[[2]] %>% summarize_fits() %>% plotly_graph(print_highlow = TRUE, fix=T)
save_image(final_plot, file = "Bootstrap_Tree_Fit.svg")
  
