library(igraph)
library(ggraph)
library(dplyr)
library(scales)
library(pheatmap)
library(reshape2)

coexpression_network <- function(data, filename, annotation) {
    cor_matrix <- cor(data, method = "pearson", use = "pairwise.complete.obs")
    threshold <- 0.5
    graph <- graph_from_adjacency_matrix(abs(cor_matrix) > threshold, weighted = TRUE, mode = "undirected", diag = FALSE)
    edges <- get.edgelist(graph)
    weights <- sapply(1:nrow(edges), function(i) cor_matrix[edges[i, 1], edges[i, 2]])
    E(graph)$weight <- weights
    E(graph)$layout_weight <- abs(weights)  # Absolute weights for layout

    # Add positive/negative correlation attribute
    E(graph)$color <- ifelse(E(graph)$weight < 0, "blue", "red")

    # Calculate node size based on degree
    V(graph)$size <- degree(graph)

    # Create layout for the graph using weights for the Fruchterman-Reingold algorithm
    layout <- create_layout(graph, layout = "fr", weights = E(graph)$layout_weight)

    # Draw the network graph, displaying positive and negative correlations
    network_plot <- ggraph(layout) +
    geom_edge_link(aes(width = abs(weight), color = color), edge_alpha = 0.8, lineend = 'round') +
    geom_node_point(aes(size = size), color = "darkred", shape = 21, fill = "white") +
    geom_node_text(aes(label = name), repel = TRUE, size = 4) +
    theme_graph() +
    scale_edge_width(range = c(0.2, 3)) +
    scale_edge_color_manual(values = c("blue" = "blue", "red" = "red"), labels = c("red" = "Positive", "blue" = "Negative"))

    ggsave(filename, plot = network_plot, width = 12, height = 8, dpi = 600)
}

coexpression_network_with_specified_edges <- function(data, filename, only_keep = NULL) {
    cor_matrix <- cor(data, method = "pearson", use = "pairwise.complete.obs")

    if (!is.null(only_keep)) {
        # Ensure that only_keep is a dataframe with columns 'source' and 'target'
        if (!("source" %in% names(only_keep) && "target" %in% names(only_keep))) {
            stop("only_keep must be a dataframe with 'source' and 'target' columns.")
        }

        filtered_cor_matrix <- matrix(0, ncol = ncol(cor_matrix), nrow = nrow(cor_matrix),
                                      dimnames = list(rownames(cor_matrix), colnames(cor_matrix)))

        for (row in rownames(cor_matrix)) {
            for (col in colnames(cor_matrix)) {
                if (dim(only_keep[only_keep$source == row & only_keep$target == col, ])[1] > 0) {
                    filtered_cor_matrix[row, col] <- cor_matrix[row, col]
                }

                if (dim(only_keep[only_keep$source == col & only_keep$target == row, ])[1] > 0) {
                    filtered_cor_matrix[row, col] <- cor_matrix[row, col]
                }
            }
        }   
    } else {
        filtered_cor_matrix <- cor_matrix
    }

    threshold <- 0.5
    graph <- graph_from_adjacency_matrix(abs(filtered_cor_matrix) > threshold, weighted = TRUE, mode = "undirected", diag = FALSE)
    
    edges <- get.edgelist(graph)
    weights <- sapply(1:nrow(edges), function(i) filtered_cor_matrix[edges[i, 1], edges[i, 2]])
    E(graph)$weight <- weights
    E(graph)$layout_weight <- abs(weights)  # Absolute weights for layout

    # Add positive/negative correlation attribute
    E(graph)$color <- ifelse(E(graph)$weight < 0, "blue", "red")

    # Calculate node size based on degree
    V(graph)$size <- degree(graph)

    # Create layout for the graph using weights for the Fruchterman-Reingold algorithm
    layout <- create_layout(graph, layout = "fr", weights = E(graph)$layout_weight)

    # Draw the network graph, displaying positive and negative correlations
    network_plot <- ggraph(layout) +
    geom_edge_link(aes(width = abs(weight), color = color), edge_alpha = 0.8, lineend = 'round') +
    geom_node_point(aes(size = size), color = "darkred", shape = 21, fill = "white") +
    geom_node_text(aes(label = name), repel = TRUE, size = 4) +
    theme_graph() +
    scale_edge_width(range = c(0.2, 3)) +
    scale_edge_color_manual(values = c("blue" = "blue", "red" = "red"), labels = c("red" = "Positive", "blue" = "Negative"))

    ggsave(filename, plot = network_plot, width = 12, height = 8, dpi = 600)
}

coexpression_network_with_specified_edges_layout <- function(data, filename, only_keep = NULL, annotation = NULL, threshold = 0.5) {
    cor_matrix <- cor(data, method = "pearson", use = "pairwise.complete.obs")

    if (!is.null(only_keep)) {
        if (!("source" %in% names(only_keep) && "target" %in% names(only_keep))) {
            stop("only_keep must be a dataframe with 'source' and 'target' columns.")
        }

        filtered_cor_matrix <- matrix(0, ncol = ncol(cor_matrix), nrow = nrow(cor_matrix),
                                      dimnames = list(rownames(cor_matrix), colnames(cor_matrix)))

        for (row in rownames(cor_matrix)) {
            for (col in colnames(cor_matrix)) {
                if (dim(only_keep[only_keep$source == row & only_keep$target == col, ])[1] > 0 ||
                    dim(only_keep[only_keep$source == col & only_keep$target == row, ])[1] > 0) {
                    filtered_cor_matrix[row, col] <- cor_matrix[row, col]
                }
            }
        }
    } else {
        filtered_cor_matrix <- cor_matrix
    }

    graph <- graph_from_adjacency_matrix(abs(filtered_cor_matrix) > threshold, weighted = TRUE, mode = "undirected", diag = FALSE)
    
    # Remove isolated nodes (nodes with no edges)
    graph <- delete_vertices(graph, V(graph)[degree(graph) == 0])
    
    edges <- get.edgelist(graph)
    weights <- sapply(1:nrow(edges), function(i) filtered_cor_matrix[edges[i, 1], edges[i, 2]])
    E(graph)$weight <- weights
    E(graph)$layout_weight <- abs(weights)

    E(graph)$color <- ifelse(E(graph)$weight < 0, "blue", "red")
    V(graph)$size <- sqrt(degree(graph)) * 2

    if (!is.null(annotation)) {
        if (!("feature" %in% names(annotation) && "direction" %in% names(annotation))) {
            stop("annotation must be a dataframe with 'feature' and 'direction' columns.")
        }

        # Merge annotation to set node color based on direction
        node_colors <- setNames(as.character(annotation$direction), annotation$feature)
        default_color <- "grey"  # Default color for nodes not in the annotation dataframe
        node_color_values <- c("Up" = "red", "Down" = "blue", "NoDiff" = "grey")
        # V(graph)$color <- ifelse(names(V(graph)) %in% names(node_colors), node_colors[names(V(graph))], default_color)
        V(graph)$color <- ifelse(names(V(graph)) %in% names(node_colors), node_color_values[node_colors[names(V(graph))]], default_color)
    } else {
        V(graph)$color <- "white"
    }

    layout <- create_layout(graph, layout = "kk", weights = E(graph)$layout_weight)

    network_plot <- ggraph(layout) +
    geom_edge_link(aes(width = abs(weight), color = color), edge_alpha = 0.6, lineend = 'round') +
    geom_node_point(aes(size = size), color = "darkred", shape = 21, fill = V(graph)$color) +
    geom_node_text(aes(label = name, y = y - 0.05), repel = TRUE, size = 4, segment.size = 0.2, point.padding = 0.3, vjust = 2) +
    theme_graph() +
    scale_edge_width(range = c(0.1, 2), name = "Weight (abs)") +
    scale_edge_color_manual(values = c("blue" = "blue", "red" = "red"), labels = c("red" = "Positive", "blue" = "Negative")) +
    labs(edge_color = "Direction", size = "Degree") +
    theme(legend.title = element_text(size = 12),
          legend.text = element_text(size = 12)) 

    ggsave(filename, plot = network_plot, width = 12, height = 8, dpi = 600)
}



data <- read.table("data.tsv", header=TRUE, sep="\t")
rownames(data) <- data$feature

## Annotate Direction
control_columns <- grep("^F", names(data), value = TRUE)
treatment_columns <- grep("^C", names(data), value = TRUE)
control_means <- rowMeans(data[, control_columns])
treatment_means <- rowMeans(data[, treatment_columns])
fold_changes <- treatment_means / control_means
results <- data.frame(Feature = data$feature, FoldChange = fold_changes)
results$log2FC <- log2(fold_changes)
data$direction <- ifelse(results$log2FC > 0.4, "Up",
                         ifelse(results$log2FC < -0.4, "Down", "NoDiff"))

data_matrix <- data[, -c(1,2,3)]
rownames(data_matrix) <- rownames(data)
data_matrix_normalized <- data_matrix

## Sample Correlation
# data_matrix_normalized <- data_matrix %>% mutate(across(-c(1, 2, 3), scale))
sample_cor_matrix <- cor(data_matrix_normalized, method = "pearson", use = "pairwise.complete.obs")
# png("sample_heatmap.png", width = 150, height = 150, units = "cm", res = 600)
pheatmap(sample_cor_matrix, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", clustering_method = "complete", display_numbers = TRUE, fontsize_number = 8)
# dev.off()

## Feature Correlation
feature_cor_matrix <- cor(t(data_matrix_normalized), method = "pearson", use = "pairwise.complete.obs")
annotation_col <- data.frame(Category = data$category)
rownames(annotation_col) <- rownames(data)
# png("feature_heatmap.png", width = 150, height = 150, units = "cm", res = 600)
pheatmap(feature_cor_matrix, annotation_col = annotation_col, clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", clustering_method = "complete", display_numbers = TRUE, fontsize_number = 6)
# dev.off()

## Co-expression Network for Gene
gene <- data[data$category == "gene",]
gene_data <- gene[, -c(1,2,3)]
annotation <- data.frame(feature = rownames(gene), direction = gene$direction)
coexpression_network_with_specified_edges_layout(t(gene_data), "coexpression_network_gene.png", NULL, annotation)


## Co-expression Network for Gene-Microbe
gene_microbe <- data[data$category != "metabolite",]
gene_microbe_data <- gene_microbe[, -c(1,2,3)]
kept_gene_microbe <-gene_microbe[gene_microbe$direction == "Up" | gene_microbe$direction == "Down", ]
combinations <- expand.grid(source = rownames(kept_gene_microbe), target = rownames(kept_gene_microbe))
combinations <- merge(combinations, kept_gene_microbe, by.x = "source", by.y = "row.names")
combinations <- merge(combinations, kept_gene_microbe, by.x = "target", by.y = "row.names", suffixes = c(".source", ".target"))
gene_to_microbe <- combinations[combinations$category.source == "gene" & combinations$category.target == "microbe", ]
only_keep <- gene_to_microbe[, c("source", "target")]
annotation <- kept_gene_microbe[,c("feature", "direction")]
coexpression_network_with_specified_edges_layout(t(gene_microbe_data), "coexpression_network_gene_microbe.png", only_keep, annotation)


## Microbe Correlation 
microbe_metabolite <- data[data$category == "microbe" | data$category == "metabolite",]
microbe_metabolite_data <- microbe_metabolite[, -c(1,2,3)]
kept_microbe_metabolite <- microbe_metabolite[microbe_metabolite$direction == "Up" | microbe_metabolite$direction == "Down", ]
combinations <- expand.grid(source = rownames(kept_microbe_metabolite), target = rownames(kept_microbe_metabolite))
combinations <- merge(combinations, kept_microbe_metabolite, by.x = "source", by.y = "row.names")
combinations <- merge(combinations, kept_microbe_metabolite, by.x = "target", by.y = "row.names", suffixes = c(".source", ".target"))
microbe_to_metabolite <- combinations[combinations$category.source == "microbe" & combinations$category.target == "metabolite", ]
only_keep <- microbe_to_metabolite[, c("source", "target")]
annotation <- kept_microbe_metabolite[,c("feature", "direction")]
coexpression_network_with_specified_edges_layout(t(microbe_metabolite_data), "coexpression_network_microbe_metabolite.png", only_keep, annotation = annotation, threshold = 0.7)
