
# Libraries ---------------------------------------------------------------

# Data wrangling
library(dplyr)
library(data.table)
library(magrittr)
library(clipr)

# Statistics
library(lme4)
library(emmeans)

# Set work directory ------------------------------------------------------

rstudioapi::getActiveDocumentContext()[["path"]] %>%
    sub(basename(.), "", .) %>% setwd %>%
    paste0(., " -> ", getwd()) %>% cat


# Pre-processing -----------------------------------------------------------

# Read sample annotation data
sdat <- fread("../lib/sampledata.txt")

# Read sample mapping data (MENIX - Somascan IDs)
smap <- fread("../lib/smap.txt")

# R does not like only numbers as column names
smap$Matched <- paste0("SID_", smap$Somalogic)

# Add the missing initial 0's for patients...
sdat <- sdat[, c("study_id", "Studygroup", "sex", "Age")]
sdat$study_id %<>% as.character
sdat[nchar(study_id) == 3]$study_id <- paste0("000", sdat[nchar(study_id) == 3]$study_id)

# Map ids (Somalogic ids - sample annotation data)
setkey(sdat, "study_id")
setkey(smap, "Samples")
smap <- sdat[smap]

# Standardize group names and remove empty
smap$Studygroup %<>% toupper()
smap <- smap[!is.na(Studygroup) & nchar(Studygroup) > 0]

smap <- smap[, c("study_id", "Matched", "Studygroup", "sex", "Age")]

# Rename study groups
smap[Studygroup == "EARLYOA"]$Studygroup = "MILD.DEGENERATION"
smap[Studygroup == "LATEOA"]$Studygroup = "LATE.OA"

# Read LOD adjusted data
lfq <- fread("../data/lfq-LOD.adjusted.txt")

# Add remaining columns
lfq <- merge(lfq, smap, by.x = "variable", by.y = "Matched", all.x = TRUE)

# Set data types
lfq$SeqId %<>% as.factor
lfq$study_id %<>% as.factor
lfq$Studygroup %<>% as.factor
lfq$sex %<>% toupper() %>%  as.factor

lfq$Studygroup <- factor(lfq$Studygroup, levels = c("HEALTHY", "MILD.DEGENERATION", "LATE.OA") )

lfq <- lfq[!is.na(value)]

lfq = dcast(lfq, SeqId ~ variable, value.var = "value")

lfq = na.omit(lfq)

lfq = melt(lfq, id.vars = "SeqId")

lfq <- merge(lfq, smap, by.x = "variable", by.y = "Matched", all.x = TRUE)


# Analysis ----------------------------------------------------------------

set.seed(123)

unique_proteins <- lfq$SeqId %>% unique

for(i in 20:60) {
    n = i
    
    # Determine the number of chunks
    num_chunks <- ceiling(length(unique_proteins) / n)
    
    # Create a factor which will be used to split the vector
    f <- sample(rep(1:num_chunks, each = n, len = length(unique_proteins)))
    
    # Split the vector into chunks
    inds <- split(unique_proteins, f)
    
    
    print(i)
    print(inds %>% lengths)
}

n = 51

# Determine the number of chunks
num_chunks <- ceiling(length(unique_proteins) / n)

# Create a factor which will be used to split the vector
f <- sample(rep(1:num_chunks, each = n, len = length(unique_proteins)))

# Split the vector into chunks
inds <- split(unique_proteins, f)

inds %>% lengths

estimates_output_list <- list()

counter = 0

for(i in 1:length(inds)) {
    
    print(counter)
    
    #if(i > 2) break
    
    p = inds[i][[1]]
    
    lmm1 <- tryCatch(
        expr = {
            lm1 <- lmer(value ~ SeqId + Studygroup + SeqId:Studygroup + sex + Age + (1 | study_id),
                        data = lfq[SeqId %in% p], REML = T)
        },
        error = function(e){ 
            print(e)
            e
        }
    )
    
    
    if(inherits(lmm1, "error")) {
        counter = counter + 1
        next
    }
    
    emm <- tryCatch(
        expr = {
            pairs(emmeans(lmm1, ~ Studygroup | SeqId), reverse = T)
        },
        error = function(e){ 
            print(e)
            e
        }
    )
    
    if(inherits(emm, "error")) next
    
    emm2 <- as.data.table(emm)
    emm3 <- as.data.table(confint(emm))
    emm3$p.value <- emm2$p.value
    emm3$Model = i
    
    estimates_output_list <- rbindlist(list(estimates_output_list, emm3), fill = T)
    
    counter = counter + 1
    
}

