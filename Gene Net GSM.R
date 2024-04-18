# Load required libraries
library(GeneNet)
library(graph)
library(Rgraphviz)
library(readr)  # Load the readr library for reading CSV files

# Read gene expression data from CSV file
Gene_Expression <- read_csv("GeneNet/Gene Expression.csv")

# Create a numeric matrix from Gene_Expression (excluding the first column)
numeric_matrix <- as.matrix(Gene_Expression[, -1])

# Assuming 'numeric_matrix' is your gene expression matrix
# Check the dimensions of the numeric matrix
dim(numeric_matrix)

# Inspect summary of gene expression data
summary(numeric_matrix)

# Assuming 'numeric_matrix' is your gene expression matrix
boxplot(numeric_matrix, names = colnames(numeric_matrix), las = 2, main = "Gene Expression Boxplot", xlab = "Samples", ylab = "Expression Values")

# Inspect pairwise scatter plots
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste(prefix, txt, sep = "")
  if (missing(cex.cor)) cex = 0.8 / strwidth(txt)
  text(0.5, 0.5, txt, cex = cex * r)
}

# Pairs plot with smoothing
pairs(numeric_matrix[, 1:10], lower.panel = panel.smooth, upper.panel = panel.cor)

# Compute Partial Correlations and Select Relevant Edges
pcor.dyn <- ggm.estimate.pcor(numeric_matrix, method = "dynamic")
gene_net <- network.test.edges(pcor.dyn, direct = TRUE)

# Use the strongest 150 edges
gene_net <- extract.network(gene_net, method.ggm = "number", cutoff.ggm = 150)

# Check the structure of gene_net and node.labels
str(gene_net)
str(node.labels)

# Construct Graph
# Extracting Gene Symbols from node labels
gene_symbols <- node.labels

# Create a graph object
gene_graph <- network.make.graph(gene_net, gene_symbols, drop.singles = TRUE)

# Some information about the graph
cat("Number of nodes:", num.nodes(gene_graph), "\n")
cat("Correlations:\n", edge.info(gene_graph)$weight, "\n")

# Number of directed ("forward") and undirected ("none") edges
cat("Number of directed ('forward') and undirected ('none') edges:\n", table(edge.info(gene_graph)$dir), "\n")

# Well-Connected Nodes
well_connected_nodes <- sort(node.degree(gene_graph), decreasing = TRUE)[1:20]
cat("Nodes connected with many edges:\n", well_connected_nodes, "\n")

# Print descriptions of top nodes (Gene Symbol)
cat("Node Descriptions:\n", Gene_Expression$`Gene Symbol`[well_connected_nodes], "\n")

# Print descriptions of top nodes
cat("Node Descriptions:\n", gene_symbols[well_connected_nodes], "\n")

# Set attributes of some particular nodes
nodeAttrs <- list()
nodeAttrs$fillcolor <- ifelse(gene_symbols %in% gene_symbols[well_connected_nodes], "red", gray(0.95))  # highlight hub nodes

# Ensure the nodeAttrs$fillcolor vector has names
names(nodeAttrs$fillcolor) <- gene_symbols

# Set edge attributes
edi <- edge.info(gene_graph)  # edge directions and correlations
edgeAttrs <- list()
edgeAttrs$dir <- edi$dir  # set edge directions 
cutoff <- quantile(abs(edi$weight), c(0.2, 0.8))  # thresholds for line width / coloring
edgeAttrs$lty <- ifelse(edi$weight < 0, "dotted", "solid")  # negative correlation
edgeAttrs$color <- ifelse(abs(edi$weight) <= cutoff[1], "grey", "black")  # lower 20% quantile
edgeAttrs$lwd <- ifelse(abs(edi$weight) >= cutoff[2], 2, 1)  # upper 20% quantile

# Set global node and edge attributes
globalAttrs <- list()
globalAttrs$edge <- list(color = "black", lty = "solid", lwd = 1, arrowsize = 1)
globalAttrs$node <- list(fillcolor = gray(0.95), shape = "ellipse", width = 1.5, height = 1, fixedsize = FALSE)

# Plot the network
plot(gene_graph, attrs = globalAttrs, nodeAttrs = nodeAttrs, edgeAttrs = edgeAttrs, "fdp")

# Install and load the igraph package
# install.packages("igraph")
library(igraph)

# Create a graph object from the gene_net data frame
gene_graph <- graph.data.frame(gene_net, directed = TRUE)

# Plot the graph
plot(gene_graph, edge.arrow.size = 0.5, layout = layout_with_fr(gene_graph))
























