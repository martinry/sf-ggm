compute_network_metrics <- function(graph) {
    G_list <- graph$G_list
    G_common <- graph$CommonG
    
    # Convert G matrices to boolean
    G_common_bool <- G_common > 0
    
    # Create a list to store the metrics
    metrics_list <- lapply(names(G_list), function(x) {
        G <- graph_from_adjacency_matrix(G_list[[x]], mode = "undirected")
        
        # Compute the metrics
        degree_centrality <- degree(G)
        closeness_centrality <- closeness(G)
        betweenness_centrality <- betweenness(G)
        eigenvector_centrality <- eigen_centrality(G)$vector
        clustering_coef <- transitivity(G, type="local")
        
        # Compute the number of unique connections
        G_bool <- G_list[[x]] > 0
        unique_bool <- G_bool & !G_common_bool
        unique_connections <- rowSums(unique_bool)
        
        # Compute the number of shared connections
        shared_bool = G_bool & G_common_bool
        shared_connections <- rowSums(shared_bool)
        
        # Combine the metrics into a data frame
        metrics_df <- data.table(
            id = names(degree_centrality),
            degree = degree_centrality,
            closeness = closeness_centrality,
            betweenness = betweenness_centrality,
            eigenvector = eigenvector_centrality,
            clustering_coef = clustering_coef,
            unique_connections = unique_connections,
            shared_connections = shared_connections
        )
        
        return(metrics_df)
    })
    
    # Assign names to the list elements
    names(metrics_list) <- names(G_list)
    
    return(metrics_list)
}

# Function to create igraph object and remove vertices with non-shared edges
jewel_to_igraph <- function(res){
    # Get adjacency matrices from model output
    G_list <- res$G_list
    G <- res$CommonG
    
    # Create igraph object for each graph
    G_list_graph <- lapply(G_list, function(x) 
        graph_from_adjacency_matrix(x, mode = "undirected"))
    G_graph <- graph_from_adjacency_matrix(G, mode = "undirected")
    
    # Find vertices by number of shared edges between graphs
    un <- do.call(igraph::union, G_list_graph)
    
    # Get isolated vertices
    isolated <- which(igraph::degree(un) == 0)
    
    # Remove isolated vertices
    G_list_plot <- lapply(G_list_graph, function(x) delete.vertices(x, isolated))
    G_plot <- delete.vertices(G_graph, isolated)
    
    # Set edge color to gray
    E(G_plot)$color <- "gray"
    E(G_plot)$width = .3

    # Construct the difference between each graph and the intersection and color accordingly.
    G_specific <- lapply(G_list_plot, function(x) difference(x, G_plot))
    colours_k <- c("#138D75", "#E67E22")
    K <- length(colours_k)
    for (k in 1:K) {
        E(G_specific[[k]])$color <- colours_k[k]
        E(G_specific[[k]])$width <- 1.2
    }
    
    # Reassemble the intersection and class-specific edges
    G_list_coloured <- lapply(G_specific, function(x) igraph::union(x, G_plot))
    for (k in 1:K) {
        edge_colors <- edge_attr(G_list_coloured[[k]], "color_1")
        specific <- which(is.na(edge_attr(G_list_coloured[[k]], "color_1")))
        edge_colors[specific] <- edge_attr(G_list_coloured[[k]], "color_2")[specific]
        G_list_coloured[[k]] <- set_edge_attr(G_list_coloured[[k]], "color", 
                                              value = edge_colors)
    }
    
    # Our list of graphs
    return(G_list_coloured)
}

plot_graphs <- function(G,
                        min_cluster_size = 2,
                        identical_layouts = F,
                        community_detection = T,
                        vlc = 0.5,
                        highlight_ids = NULL) {
    set.seed(123)
    # Used if we want a identical layout between each plot
    if(identical_layouts) {
        l <- layout_nicely(G[[1]])
        coords <- tkplot(G[[1]], layout = l, vertex.size = 0)
        l <- tkplot.getcoords(coords)
    } else {
        l <- layout_nicely
    }
    
    # Set font and margins
    font_add_google("Poppins", "myfont")
    showtext_auto()
    par(family = "myfont")
    par(mar = c(.5, .5, 2, .5) + .1)
    par(mfrow = c(1, 2))
    
    # Plot each graph
    for(i in seq_along(G)) {
        m <- names(G[i])
        g = G[[m]]
        
        set.seed(123)
        
        if(!is.null(highlight_ids)){
            V(g)$name = ifelse(grepl(paste(highlight_ids, collapse = "|"), V(g)$name), V(g)$name, "")
        }
        
        if(community_detection){
            # Detect communities
            lou <- cluster_louvain(g)
            
            # Count community sizes and 'remember' communities with at least min_cluster_size members
            c_keep_ids <- as.numeric(names(sizes(lou)[sizes(lou) >= min_cluster_size]))
            
            # Create color vector
            cols = alpha(rep("black", uniqueN(lou$membership)), .1)
            
            # Color only communities we want to keep
            cols[c_keep_ids] <- alpha(brewer.pal(length(c_keep_ids), "Set1")[1:length(c_keep_ids)], .4)
            
            
            # Plot the graph
            plot(lou, g, layout = l,
                 mark.border = cols,
                 mark.col = cols,
                 vertex.shape = 'none',
                 vertex.size = 0,
                 vertex.frame.color = alpha("white", 0),
                 vertex.label.cex = vlc, 
                 vertex.label.color = "black",
                 vertex.label.family = "myfont",
                 main = m)
        } else {
            # Plot the graph (without mark.border and mark.col)
            plot(g, layout = l,
                 vertex.shape = 'none',
                 vertex.size = 0,
                 vertex.frame.color = alpha("white", 0),
                 vertex.label.cex = vlc, 
                 vertex.label.color = "black",
                 vertex.label.family = "myfont",
                 main = m)
            
        }
    }
}

jewel_unique_edges <- function(res){
    # Get adjacency matrices from model output
    G_list <- res$G_list
    G <- res$CommonG
    
    # Create igraph object for each graph
    G_list_graph <- lapply(G_list, function(x) 
        graph_from_adjacency_matrix(x, mode = "undirected"))
    G_graph <- graph_from_adjacency_matrix(G, mode = "undirected")
    
    # Find vertices by number of shared edges between graphs
    un <- do.call(igraph::union, G_list_graph)
    
    # Get isolated vertices
    isolated <- which(degree(union) == 0)
    
    # Remove isolated vertices
    G_list_plot <- lapply(G_list_graph, function(x) delete.vertices(x, isolated))
    G_plot <- delete.vertices(G_graph, isolated)
    
    # Set edge color to gray
    E(G_plot)$color <- "gray"
    E(G_plot)$width = .3
    
    # Construct the difference between each graph and the intersection.
    G_shared = lapply(G_list_plot, function(x) intersect(x, G_plot))
    G_specific <- lapply(G_list_plot, function(x) difference(x, G_plot))
    colours_k <- c("#138D75", "#E67E22")
    K <- length(colours_k)
    for (k in 1:K) {
        E(G_specific[[k]])$color <- colours_k[k]
        E(G_specific[[k]])$width <- 1.2
    }
    
    # Return the list of unique graphs
    return(G_specific)
}

get_node_ids <- function(res, option = "unique") {
    # Get adjacency matrices from model output
    G_list <- res$G_list
    G <- res$CommonG
    
    # Convert G matrices to boolean
    G_common <- G > 0
    
    # Depending on the option selected, return either shared or unique node ids
    if (option == "unique") {
        # Compute unique ids for each group
        unique_ids <- lapply(names(G_list), function(x) {
            G_group <- G_list[[x]] > 0
            unique_bool <- G_group & !G_common
            rownames(G_group)[rowSums(unique_bool) > 0]
        })
        return(unique_ids)
    } else if (option == "shared") {
        # Compute shared ids
        shared_bool <- Reduce(`&`, lapply(G_list, function(x) x > 0))
        shared_ids <- rownames(G_common)[rowSums(shared_bool) > 0]
        return(shared_ids)
    } else {
        stop("Invalid option selected. Choose either 'unique' or 'shared'.")
    }
}

map_node_names <- function(graph, mapping_df) {
    G_list <- graph$G_list
    
    # Replace node names in adjacency matrices
    G_list <- lapply(G_list, function(g) {
        row_map_df <- merge(data.frame(SeqId = rownames(g)), mapping_df, by = "SeqId", all.x = TRUE)
        col_map_df <- merge(data.frame(SeqId = colnames(g)), mapping_df, by = "SeqId", all.x = TRUE)
        
        # Replace missing mappings with "NA"
        row_map_df[is.na(row_map_df$EntrezGeneSymbol), "EntrezGeneSymbol"] <- "NA"
        col_map_df[is.na(col_map_df$EntrezGeneSymbol), "EntrezGeneSymbol"] <- "NA"
        
        rownames(g) <- row_map_df$EntrezGeneSymbol
        colnames(g) <- col_map_df$EntrezGeneSymbol
        return(g)
    })
    
    # Replace node names in the CommonG matrix
    CommonG <- graph$CommonG
    row_map_df <- merge(data.frame(SeqId = rownames(CommonG)), mapping_df, by = "SeqId", all.x = TRUE)
    col_map_df <- merge(data.frame(SeqId = colnames(CommonG)), mapping_df, by = "SeqId", all.x = TRUE)
    
    # Replace missing mappings with "NA"
    row_map_df[is.na(row_map_df$EntrezGeneSymbol), "EntrezGeneSymbol"] <- "NA"
    col_map_df[is.na(col_map_df$EntrezGeneSymbol), "EntrezGeneSymbol"] <- "NA"
    
    rownames(CommonG) <- row_map_df$EntrezGeneSymbol
    colnames(CommonG) <- col_map_df$EntrezGeneSymbol
    
    return(list(G_list = G_list, CommonG = CommonG))
}

map_node_names_combined <- function(graph, mapping_df, idtype = "EntrezGeneSymbol") {
    G_list <- graph$G_list
    
    # Replace node names in adjacency matrices
    G_list <- lapply(G_list, function(g) {
        row_map_df <- merge(data.frame(SeqId = rownames(g)), mapping_df, by = "SeqId", all.x = TRUE)
        col_map_df <- merge(data.frame(SeqId = colnames(g)), mapping_df, by = "SeqId", all.x = TRUE)
        
        # Replace missing mappings with "NA"
        row_map_df[is.na(row_map_df[, idtype]), idtype] <- "NA"
        col_map_df[is.na(col_map_df[, idtype]), idtype] <- "NA"
        
        # Combine EntrezGeneSymbol and SeqId
        rownames(g) <- paste(row_map_df[, idtype], row_map_df$SeqId, sep = "\n")
        colnames(g) <- paste(col_map_df[, idtype], col_map_df$SeqId, sep = "\n")
        return(g)
    })
    
    # Replace node names in the CommonG matrix
    CommonG <- graph$CommonG
    row_map_df <- merge(data.frame(SeqId = rownames(CommonG)), mapping_df, by = "SeqId", all.x = TRUE)
    col_map_df <- merge(data.frame(SeqId = colnames(CommonG)), mapping_df, by = "SeqId", all.x = TRUE)
    
    # Replace missing mappings with "NA"
    row_map_df[is.na(row_map_df[, idtype]), idtype] <- "NA"
    col_map_df[is.na(col_map_df[, idtype]), idtype] <- "NA"
    
    # Combine EntrezGeneSymbol and SeqId
    rownames(CommonG) <- paste(row_map_df[, idtype], row_map_df$SeqId, sep = "\n")
    colnames(CommonG) <- paste(col_map_df[, idtype], col_map_df$SeqId, sep = "\n")
    
    return(list(G_list = G_list, CommonG = CommonG))
}



plot_graphs_ext <- function(G,
                        min_cluster_size = 2,
                        identical_layouts = F,
                        community_detection = T,
                        vlc = 0.5,
                        centrality = "degree",
                        highlight_ids = NULL) {
    set.seed(123)
    # Used if we want a identical layout between each plot
    if(identical_layouts) {
        l <- layout_nicely(G[[1]])
        coords <- tkplot(G[[1]], layout = l, vertex.size = 0)
        l <- tkplot.getcoords(coords)
    } else {
        l <- layout_nicely
    }
    
    # Set font and margins
    font_add_google("Poppins", "myfont")
    showtext_auto()
    par(family = "myfont")
    par(mar = c(.5, .5, 2, .5) + .1)
    par(mfrow = c(1, 2))
    
    # Plot each graph
    for(i in seq_along(G)) {
        m <- names(G[i])
        g = G[[m]]
        
        set.seed(123)

        if(!is.null(highlight_ids)){
            V(g)$name = ifelse(grepl(paste(highlight_ids, collapse = "|"), V(g)$name), V(g)$name, "")
        }
        

        
        if(centrality == "eigenvector") {
            vs = eigen_centrality(g)$vector
            vs[is.nan(vs)] = 0
            vs = (0.5 + vs) * 5
        } else if(centrality == "closeness") {
            vs = closeness(g)
            vs[is.nan(vs)] = 0
            vs = vs * 1000
            vs = log2(vs)
            } else if(centrality == "clustering_coef") {
            vs = transitivity(g, type="local")
            vs[is.nan(vs)] = 0
            vs = (0.5 + vs) * 5
        } else {
            vs = get(centrality)(g)
            vs = log2(1 + vs)
        }

        
        if(community_detection){
            # Detect communities
            lou <- cluster_louvain(g)
            
            # Count community sizes and 'remember' communities with at least min_cluster_size members
            c_keep_ids <- as.numeric(names(sizes(lou)[sizes(lou) >= min_cluster_size]))
            
            # Create color vector
            cols = alpha(rep("black", uniqueN(lou$membership)), .1)
            
            # Color only communities we want to keep
            cols[c_keep_ids] <- alpha(brewer.pal(length(c_keep_ids), "Set1")[1:length(c_keep_ids)], .4)
            
            
            # Plot the graph
            set.seed(123)
            plot(lou, g, layout = l,
                 mark.border = cols,
                 mark.col = cols,
                 #vertex.shape = 'none',
                 vertex.size= vs,
                 #vertex.frame.color = alpha("white", 0),
                 vertex.label.cex = vlc, 
                 vertex.label.color = "black",
                 vertex.label.family = "myfont",
                 main = m)
        } else {
            set.seed(123)
            # Plot the graph (without mark.border and mark.col)
            plot(g, layout = l,
                 #vertex.shape = 'none',
                 vertex.size= vs,
                 #vertex.frame.color = alpha("white", 0),
                 vertex.label.cex = vlc, 
                 vertex.label.color = "black",
                 vertex.label.family = "myfont",
                 main = m)
            
        }
    }
}
