suppressMessages(library(tidyverse))
library(readxl)
library(ggridges)
library(argparse)

NUM_GAMMA_SAMPLES <- 1000
VOLUME_KEY <- "volume"
MAX_PROTEINS <- 100

# Set up Command-Line Interface
parser <- ArgumentParser()
parser$add_argument("expression_csv", help = "Input expression data file")
parser$add_argument("chosen_csv",
                    help = "CSV of chosen genes and associated wcEcoli keys")
parser$add_argument("output", help = "Output folder")
parser$add_argument("experiment_id", help = "ID of experiment")
args <- parser$parse_args()

# Load Chosen Genes and Proteins to Compare
chosen <- suppressMessages(read_csv(args$chosen_csv))

# Load Reference Data
taniguchi_s6_path <- file.path("assets", "TableS6.xls")
taniguchi_s6 <- read_excel(taniguchi_s6_path)
gamma_parameters <- taniguchi_s6[, c("Gene Name", "A_Protein",
                                     "B_Protein")]

# Load Simulation Data
data <- suppressMessages(read_csv(args$expression_csv))

# Set up vectors we will add data to
proteins <- c()
expressions <- c()
trials <- c()
data_sources <- c()
ks_stats <- c()
mean_diffs <- c()
ks_proteins <- c()

# Fill vectors with data
num_proteins <- min(nrow(chosen), MAX_PROTEINS)
for (i in seq_len(num_proteins)) {
    # Get gene and protein names
    gene <- chosen["taniguchi_gene_name"][[i, 1]]
    protein_key <- chosen["wcecoli_protein_key"][[i, 1]]

    # Get expression levels
    protein_concentrations <- data[protein_key][[1]]  # Units of counts/fL
    volumes <- data[VOLUME_KEY][[1]]  # Units of fL/cell
    stopifnot(length(protein_concentrations) == length(volumes))
    protein_counts <- protein_concentrations * volumes  # Units of counts/cell
    num_cells <- length(protein_concentrations)

    proteins <- c(proteins, rep(protein_key, num_cells))
    expressions <- c(expressions, protein_counts)
    trials <- c(trials, rep(1, num_cells))
    data_sources <- c(data_sources, rep("simulated", num_cells))

    # Get Taniguchi Parameters
    gene_parameters_map <- gamma_parameters$"Gene Name" == gene
    a <- gamma_parameters[gene_parameters_map, "A_Protein"][[1]]
    b <- gamma_parameters[gene_parameters_map, "B_Protein"][[1]]

    # Sample from Taniguchi Distribution
    taniguchi_samples <- rgamma(n = NUM_GAMMA_SAMPLES, shape = a,
                                scale = b)
    proteins <- c(proteins, rep(protein_key, NUM_GAMMA_SAMPLES))
    expressions <- c(expressions, taniguchi_samples)
    trials <- c(trials, rep(1, NUM_GAMMA_SAMPLES))
    data_sources <- c(data_sources, rep("expected", NUM_GAMMA_SAMPLES))

    # Run KS Test
    ks_stat_obj <- ks.test(protein_counts, "pgamma", shape = a, scale = b)
    ks_stat <- ks_stat_obj$statistic
    ks_stats <- c(ks_stats, ks_stat)
    ks_proteins <- c(ks_proteins, protein_key)

    # Find distance between distribution means
    mean_diff <- abs(mean(taniguchi_samples) - mean(protein_counts))
    mean_diffs <- c(mean_diffs, mean_diff)
}

# Plot KS statistics
stat_sort <- tibble(protein = ks_proteins, value = ks_stats + mean_diffs)
stat_sort <- arrange(stat_sort, desc(value))
protein_order = stat_sort$protein

stat_data <- tibble(protein = c(ks_proteins, ks_proteins),
                    ks = c(ks_stats, mean_diffs / 10000),
                    stat = c(rep("KS Statistic", length(ks_proteins)),
                             rep("Difference Between Means", length(ks_proteins))))
#stat_data <- arrange(stat_data, desc(ks))
stat_data$protein <- factor(stat_data$protein, levels = protein_order)
stats_plot <- ggplot(data = stat_data) +
    geom_bar(mapping = aes(x = protein, y = ks, fill = stat), stat = "identity") +
    labs(title = "KS Test Between Experimental and Simulated Protein Counts",
         subtitle = paste("From Experiment ", args$experiment_id)) +
    xlab("Protein") + ylab("KS Statistic + Difference Between Means / 10,000") + coord_flip()
ggsave(file.path(args$output, "expression_ks_stats.pdf"), stats_plot,
       height = num_proteins * 0.5, units = "cm", limitsize = FALSE)

# Plot Distributions
transformed <- tibble(protein = proteins, expression = expressions,
                      trial = trials, data_source = data_sources)
transformed$protein <- factor(transformed$protein, levels = protein_order)
expression_plot <- ggplot(transformed,
                          aes(x = expression, y = protein,
                              fill = data_sources,
                              linetype = data_sources,
                              point_color = data_sources)) +
    geom_density_ridges(
        jittered_points = TRUE,
        position = position_raincloud(height = 0.3),
        point_alpha = 0.3,
        alpha = 0.5, scale = 0.7, rel_min_height = 0,
    ) +
    labs(title = "Protein Concentration Distributions",
         subtitle = paste("From Experiment ", args$experiment_id)) +
    xlab("Protein Counts per Cell") + ylab("Protein")

expression_range <- max(expressions) - min(expressions)
suppressMessages(ggsave(file.path(args$output,
                                  "expression_distributions.pdf"),
                        expression_plot, height = num_proteins * 10,
                        width = expression_range * 0.001,
                        units = "cm", limitsize = FALSE))

# Plot Expression Data for AmpC and AcrAB-TolC
protein <- c()
expression <- c()

for (value in data$"EG10040-MONOMER[p]") {
    protein <- c(protein, "AmpC")
    expression <- c(expression, value)
}
for (value in data$"TRANS-CPLX-201[s]") {
    protein <- c(protein, "AcrAB-TolC")
    expression <- c(expression, value)
}
transformed <- tibble(protein = protein, expression = expression)
expression_plot <- ggplot(transformed, aes(x = expression, y = protein)) +
    geom_density_ridges(
        jittered_points = TRUE,
        position = position_points_jitter(width = 0.05, height = 0),
        point_shape = "|", point_size = 3, point_alpha = 1,
        alpha = 0.7,
    ) +
    labs(title = "Protein Concentration Distributions",
         subtitle = paste("From Experiment ", args$experiment_id)) +
    xlab("Protein Concentration (counts/fL)") + ylab("Protein")

suppressMessages(ggsave(file.path(args$output,
                                  "antibiotics_expression.pdf"),
                        expression_plot))
