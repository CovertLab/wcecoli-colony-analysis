library(argparse)
library(phytools)

parser <- ArgumentParser()
parser$add_argument("phylogeny", help = "Newick file with phylogenetic tree")
parser$add_argument("survival",
                    help = "CSV with each agent's survival outcome")
args <- parser$parse_args()

tree <- read.newick(file = args$phylogeny)
survival <- read.csv(args$survival)

x <- survival[["survival"]]
names(x) <- survival[["agents"]]

phylosig(tree, x, method = "K", test = TRUE, nsim = 1e5)
phylosig(tree, x, method = "lambda", test = TRUE)
