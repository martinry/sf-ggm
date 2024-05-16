library(igraph)
library(stringr)
library(clipr)
library(dplyr)
library(data.table)
library(magrittr)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(plotly)
library(showtext)

rstudioapi::getActiveDocumentContext()[["path"]] %>%
    sub(basename(.), "", .) %>% setwd %>% paste0(., " -> ", getwd()) %>% cat

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

    
    
    ccols = list(brewer.pal(12, "Paired")[1:6], brewer.pal(12, "Paired")[6:12])
    
    # Plot each graph
    for(i in seq_along(G)) {
        m <- names(G[i])
        g = G[[m]]
        
        set.seed(123)
        
        # Set node frame colors based on if the protein is up- or downregulated
        frame_colors <- rep("lightgray", vcount(g))
        upregulated <- est1$SeqId[est1$estimate > 0]
        downregulated <- est1$SeqId[est1$estimate <= 0]
        frame_colors[match(upregulated, V(g)$name)] <- "blue"
        frame_colors[match(downregulated, V(g)$name)] <- "red"
        
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
            cols[c_keep_ids] <- alpha(ccols[[i]][1:length(c_keep_ids)], .7)
        
            V(g)$label <- NA
            V(g)$color = "lightgray"
            
            # Plot the graph
            
            set.seed(123)
            lout <- layout.fruchterman.reingold(g, niter=10000)

            community_membership <- membership(lou)
            node_community_mapping <- data.table(node = names(community_membership), community = community_membership)
            
            set.seed(123)
            plot(lou, g, layout = lout,
                 mark.border = cols,
                 mark.col = cols,
                 col = "lightgray",
                 vertex.color = "black",
                 vertex.size = 4,
                 vertex.frame.color = frame_colors,
                 vertex.frame.size = 2,
                 vertex.label.cex = vlc, 
                 vertex.label.color = "black",
                 vertex.label.family = "myfont",
                 main = m)
            
            if(i == 1) {
                text(community_coords.1$x, community_coords.1$y, labels =  LETTERS[1:5])
            } else if(i == 2) {
                text(community_coords.2$x, community_coords.2$y, labels =  LETTERS[6:11])
            }
            
            
        } else {
            V(g)$label <- NA
            V(g)$color = "lightgray"
            # Plot the graph (without mark.border and mark.col)
            set.seed(123)
            plot(g, layout = layout.fruchterman.reingold(g, niter=10000),
                 vertex.size = 4,
                 vertex.frame.color = frame_colors,
                 vertex.frame.size = 2,
                 vertex.label.cex = vlc, 
                 vertex.label.color = "black",
                 vertex.label.family = "myfont",
                 main = m)
            
        }
    }
}

# > community_coords.1
# $x
# [1] -0.400219376 -0.599562361  0.213143653  0.627163697 -0.001533408
# 
# $y
# [1]  0.7774376  0.4170869  0.2484120 -0.3496169 -1.1546559
# 
# > community_coords.2
# $x
# [1] -1.13165479 -0.64863140 -0.97064699 -0.67163252 -0.07360356  0.02606793
# 
# $y
# [1] -0.1119388  0.1487405 -0.5642940 -0.3036147 -0.3879521 -1.1469889

