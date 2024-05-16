
# Libraries ---------------------------------------------------------------

library(dplyr)
library(data.table)
library(magrittr)

# Set work directory ------------------------------------------------------

rstudioapi::getActiveDocumentContext()[["path"]] %>%
    sub(basename(.), "", .) %>% setwd %>%
    paste0(., " -> ", getwd()) %>% cat

# Pre-processing ----------------------------------------------------------

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

# Read normalized quant data
lfq <- fread("../data/normSMP.txt")

# Remove proteins missing accession
lfq <- lfq[nchar(UniProt) > 0]

# Keep only human entries
lfq <- lfq[Organism == "Human"]

somalogic.samples <- fread("../lib/somalogic.sampledata.txt")
somalogic.samples <- somalogic.samples[, c("SampleId", "SampleType")]
somalogic.samples$SampleId <- paste0("SID_", somalogic.samples$SampleId)

somalogic.samples <- somalogic.samples[SampleType == "Buffer"]

cols <- c("SeqId", somalogic.samples$SampleId)

# mean_of_buffer+5*SD_of_buffer

# Get buffer samples
buffer <- lfq[, grep("SID_220114", names(lfq)), with = F]
buffer_names <- paste0("Buffer", seq(1:ncol(buffer)))
names(buffer) <- buffer_names

buffer %<>% log2

buffer$SeqId <- lfq$SeqId

my_dt_rowmean <- buffer[, rowMeans(.SD), .SDcols = buffer_names]
my_dt_rowsd <- buffer[, matrixStats::rowSds(as.matrix(.SD)), .SDcols = buffer_names]

buffer$LOD <- my_dt_rowmean + 5 * my_dt_rowsd

# Make lfq < LoD missing

# Uniprot accessions + Somascan sample ids will be our column names
cols_to_keep <- c("SeqId", smap$Matched)
lfq <- lfq[, ..cols_to_keep]

lfq[, 2:ncol(lfq)] %<>% log2

lfq <- melt(lfq, id.vars = "SeqId")

buffer <- buffer[, c("SeqId", "LOD")]

lfq <- merge(lfq, buffer, by = "SeqId")

lfq[value < LOD]$value <- NA

lfq <- lfq[, -"LOD"]

# Write to file
fwrite(lfq, file = "../data/lfq-LOD.adjusted.txt", sep = "\t")