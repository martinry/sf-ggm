


est1 = fread("../data/estimates.EARLYOA_vs_HEALTHY.UniProt.Entrez.2023-06-01.txt")

est1$Selected = F

est1[order(abs(estimate), decreasing = T)][1:800]$Selected = T


fwrite(est1, file = "../results/S1.txt", sep = "\t")


selected = est1[Selected == T]



netw = clipr::read_clip_tbl() %>% as.data.table
netw[id %in% V(g)$name] %>% clipr::write_clip()
