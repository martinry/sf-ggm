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
        l <- layout.fruchterman.reingold(G[[1]], niter=10000)
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
            
            #V(g)$label <- NA
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
            plot(g, layout = l,
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


est1 = fread("../data/estimates.EARLYOA_vs_HEALTHY.UniProt.Entrez.2023-06-01.txt")

m4 = readRDS("../data/m4.L1_0.1_L2_0.001_ITER_10K_STABILITYSS_1K.rds")

g = m4 %>% jewel_to_igraph()
#g = m4 %>% jewel_unique_edges()




pdf("../results/m4.identical-layouts.pdf", width = 6, height = 7, pointsize = 6, compress = F)
g %>% plot_graphs(community_detection = F, identical_layouts = T)
dev.off()





G = m4 %>% jewel_to_igraph()



m <- names(G[1])
g = G[[m]]

pclasses <- data.table(SeqId = V(g)$name)


pclasses$UniProt <- est1$UniProt[match(pclasses$SeqId, est1$SeqId)]
pclasses$GeneSymbol <- est1$EntrezGeneSymbol[match(pclasses$SeqId, est1$SeqId)]

pclasses$Classification = pclasses$UniProt %>% panther.classify()

pclasses.healthy = pclasses[Group == "Healthy"]
pclasses.healthy = merge(pclasses.healthy, healthy.metrics, by.x = "SeqId", by.y = "id")

pclasses.oa = pclasses[Group == "Mild OA"]
pclasses.oa = merge(pclasses.oa, oa.metrics, by.x = "SeqId", by.y = "id")

pclasses = rbindlist(list(pclasses.healthy, pclasses.oa))

merge(pclasses, est1, by = "SeqId", all.x = T)

pclasses2 = merge(pclasses, est1, by = "SeqId", all.x = T)

pclasses2$Regulation = "Downregulated"
pclasses2[estimate > 0]$Regulation = "Upregulated"

pclasses2 = pclasses2[, c("Community", "SeqId", "UniProt.x", "GeneSymbol", "Classification", "degree_centrality", "betweenness_centrality", "estimate", "lower.CL", "upper.CL")]

names(pclasses2) = c("Community", "SeqId", "UniProt", "GeneSymbol", "Classification", "Degree", "Betweenness", "Log2 fold-change", "lower.CL", "upper.CL")

fwrite(pclasses2, "../results/m4.unique_edges.communities.classification.txt", sep = "\t")








# Read Jewel output
m4 = readRDS("m4.L1_0.1_L2_0.001_ITER_10K_STABILITYSS_1K.rds")

# Read mapping file
up = fread(up, "up.txt")

# Source mapping + node id script
source("get_node_ids.R")

# Map node names (from SeqId to SeqId+Entrez)
m.entrez = map_node_names_combined(m4, up)

# Get nodes for healthy with unique edges
m.healthy = get_node_ids(m.entrez, option = "unique")
m.healthy = m.healthy[[1]]

# Get nodes for OA with unique edges
m.oa = get_node_ids(m.entrez, option = "unique")
m.oa = m.oa[[2]]

# Get nodes with edges in healthy and OA
shared = get_node_ids(m.entrez, option = "shared")

c(m.oa, m.healthy, shared) %>% uniqueN


m.healthy %>% sub("\n", "\t", .) %>% clipr::write_clip()
m.oa  %>% sub("\n", "\t", .) %>% clipr::write_clip()
shared  %>% sub("\n", "\t", .) %>% clipr::write_clip()

shared[shared %in% intersect(m.healthy, m.oa)]

shared %>% uniqueN




metrics = compute_network_metrics(m4)

healthy.metrics = metrics$HEALTHY

oa.metrics = metrics$EARLYOA

healthy.metrics = merge(healthy.metrics, est1, by.x = "id", by.y = "SeqId")
oa.metrics = merge(oa.metrics, est1, by.x = "id", by.y = "SeqId")


# Compute the metrics
compute_metrics = function(g) {
    degree_centrality <- degree(g)
    betweenness_centrality <- betweenness(g)
    data.table(id = names(degree_centrality),
               degree_centrality = degree_centrality,
               betweenness_centrality = betweenness_centrality)
}


g = m4 %>% jewel_unique_edges()

healthy.metrics = compute_metrics(g$HEALTHY)
oa.metrics = compute_metrics(g$EARLYOA)



oa.classes = pclasses2[Community %in% LETTERS[6:11]]

merge(oa.classes, oa.metrics, by.x = "SeqId", by.y = "id") %>% View
