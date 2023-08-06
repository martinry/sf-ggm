

m4 = readRDS("../data/m4.L1_0.1_L2_0.001_ITER_10K_STABILITYSS_1K.rds")

G = m4 %>% jewel_unique_edges()

ccols = list(brewer.pal(12, "Paired")[1:6], brewer.pal(12, "Paired")[6:12])

par(family = "myfont")
par(mar = c(.5, .5, 2, .5) + .1)
par(mfrow = c(1, 2))

m <- names(G[1])
g = G[[m]]

set.seed(123)

# Set node frame colors based on if the protein is up- or downregulated
frame_colors <- rep("lightgray", vcount(g))
upregulated <- est1$SeqId[est1$estimate > 0]
downregulated <- est1$SeqId[est1$estimate <= 0]
frame_colors[match(upregulated, V(g)$name)] <- "blue"
frame_colors[match(downregulated, V(g)$name)] <- "red"

# Detect communities
lou <- cluster_louvain(g)

# Count community sizes and 'remember' communities with at least min_cluster_size members
c_keep_ids <- as.numeric(names(sizes(lou)[sizes(lou) > min_cluster_size]))

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


community_coords = locator()

community_coords.1 = community_coords

text(community_coords.1$x, community_coords.2$y, labels =  LETTERS[1:5])


community_coords.2 = community_coords

text(community_coords.2$x, community_coords.2$y, labels =  LETTERS[6:12])







soma = fread("../data/normSMP.txt")
soma = soma[Organism == "Human" & nchar(UniProt) > 0]

soma = soma$UniProt

soma = unique(soma)

soma.dt = data.table("UniProt" = soma, "Classification" = NA)

soma.dt$Classification = soma.dt$UniProt %>% panther.classify(.)

soma.df$Classification = sapply(soma.df$Classification, unlist)
soma.df$Classification %<>% unlist

tmp1 = soma.df$Classification %>% unlist


fwrite(file = "../results/all.somalogic.uniprot.classification.txt", soma.dt, sep = "\t")
