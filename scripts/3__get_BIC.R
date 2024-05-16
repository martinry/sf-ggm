# Libraries ---------------------------------------------------------------

# Data wrangling
library(dplyr)
library(data.table)
library(magrittr)
library(clipr)

# Plotting
library(ggplot2)
library(ggtext)
library(viridis)
library(RColorBrewer)
library(plotly)
library(stringr)
library(igraph)
library(showtext)

# Statistics
library(lme4)
library(emmeans)
library(jewel)

# SetWD -------------------------------------------------------------------

rstudioapi::getActiveDocumentContext()[["path"]] %>%
    sub(basename(.), "", .) %>% setwd %>% paste0(., " -> ", getwd()) %>% cat

est = fread("../results/selected.aptamers.txt")

apt = est$SeqId


# ID Mapping

up = fread("../data/normSMP.txt")

up = up[, c("SeqId", "UniProt", "EntrezGeneSymbol")]

up[UniProt == "P56192"]$EntrezGeneSymbol = "MARS1"


# Preprocessing -----------------------------------------------------------

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

lfq <- fread("../data/lfq-LOD.adjusted.txt")

# Add remaining columns
lfq <- merge(lfq, smap, by.x = "variable", by.y = "Matched", all.x = TRUE)

# Set data types
lfq$SeqId %<>% as.factor
lfq$study_id %<>% as.factor
lfq$Studygroup %<>% as.factor
lfq$sex %<>% toupper() %>%  as.factor

lfq$Studygroup <- factor(lfq$Studygroup, levels = c("HEALTHY", "EARLYOA", "LATEOA") )

lfq <- lfq[!is.na(value)]

lfq = lfq[SeqId %in% apt]

lfq <- dcast(lfq, formula = SeqId ~ variable, value.var = "value")

# Create a list which will fit our two groups, Early and Healthy
X <- vector(mode = "list", length = 2)

# Name vectors within list according to our groups
names(X) <- c("HEALTHY", "EARLYOA")

# Find sample names by group (according to our mapping variable) and assign
# corresponding quant data to each vector in the list
X[[1]] <- t(lfq[, names(lfq) %in% smap$Matched[smap$Studygroup == "HEALTHY"], with = F])
X[[2]] <- t(lfq[, names(lfq) %in% smap$Matched[smap$Studygroup == "EARLYOA"], with = F])

# Add Uniprot accessions as column name (as we have transposed the matrix)
colnames(X[[1]]) <- colnames(X[[2]]) <- lfq$SeqId

# Number of samples (rows), proteins (columns) in each vector?
sapply(X, dim)




model_eval = function(l1, l2) {
    myjewel(X, lambda1 = l1, lambda2 = l2, verbose = T, maxIter = 10000)
}

jewels = list()
l1 = seq(from = 0.01, to = 0.25, by = 0.02)

for(i in 1:length(l1)) {
    jewels[[i]] = model_eval(l1 = l1[i], l2 = 0.001)
}

bics1 <- data.table(
    Model = paste("Model ", seq(1:length(l1))),
    Lambda1 = l1,
    Lambda2 = rep(0.001, length(l1)),
    Iterations = rep(10000, length(l1)),
    BIC = sapply(jewels, function(x) x$BIC)
)

m1 = jewel(X, lambda1 = 0.0001, lambda2 = 0.0001, maxIter = 10000, stability = T, stability_nsubsets = 1000)
m2 = jewel(X, lambda1 = 0.01, lambda2 = 0.0001, maxIter = 10000, stability = T, stability_nsubsets = 1000)
m3 = jewel(X, lambda1 = 0.1, lambda2 = 0.001, maxIter = 10000, stability = T, stability_nsubsets = 1000)


library(RcppEigen)
library(SMUT)

myjewel = function (X, lambda1, lambda2 = NULL, Theta = NULL, W = NULL, 
                    tol = 0.01, maxIter = 10000, verbose = TRUE) 
{
    X <- mapply(function(y) scale(y), X, SIMPLIFY = FALSE)
    K <- length(X)
    n_k <- sapply(X, function(x) dim(x)[1])
    p <- dim(X[[1]])[2]
    vars <- colnames(X[[1]])
    if (is.null(vars)) {
        warning("Colnames of some datasets are empty. Be aware that without variable names there is no check if their order (columns order) is the same between datasets. That can lead to wrong estimation.")
    }
    else {
        if (sum(sapply(2:K, function(k) sum(colnames(X[[k]]) %in% 
                                            vars) == p)) != (K - 1) | sum(sapply(2:K, function(k) sum(vars %in% 
                                                                                                      colnames(X[[k]])) == p)) != (K - 1)) {
            stop("Variables don't match across classes. Check colnames of each element of X.")
        }
        else {
            for (k in 1:K) {
                X[[k]] <- X[[k]][, vars]
            }
        }
    }
    Xl <- do.call(rbind, X)
    nindex <- rep(1:K, n_k)
    index <- c(1:K)
    if (is.null(Theta)) {
        Theta <- Xi <- Gamma <- matrix(0, nrow = (p - 1) * K, 
                                       ncol = p)
        Active <- matrix(TRUE, nrow = (p - 1), ncol = p)
        Active_K_long <- matrix(TRUE, nrow = (p - 1) * K, ncol = p)
        R_specific <- R_common <- Xl
        Gamma_list <- rep(list(NA), K)
        Xi_list <- rep(list(NA), K)
    }
    else {
        r <- mapply(function(x, y) x - eigenMapMatMult(x, y), 
                    X, Theta, SIMPLIFY = FALSE)
        R_specific <- R_common <- do.call(rbind, r)
        remove(r)
        Active <- lapply(Theta, function(x) x != 0)
        Active <- Reduce("+", Active)
        Active <- (Active == 3)
        Xi <- lapply(Theta, function(x) x[Active])
        Gamma <- mapply(function(x, y) x - y, Theta, Xi)
        Active_K <- lapply(Gamma, function(x) jewel:::removeDiagonal(x != 
                                                                         0))
        Active_K_long <- do.call(rbind, Active_K)
        remove(Active_K)
        Gamma <- lapply(Gamma, jewel:::removeDiagonal)
        Gamma <- do.call(rbind, Gamma)
        Xi <- lapply(Xi, jewel:::removeDiagonal)
        Xi <- do.call(rbind, Xi)
        Gamma_list <- rep(list(NA), K)
        Xi_list <- rep(list(NA), K)
    }
    if (is.null(W)) {
        W <- matrix(1, nrow = (p - 1) * K, ncol = p)
        colnames(W) <- vars
    }
    else {
        if (sum(sapply(W, function(x) !is.null(colnames(x)))) != 
            K) {
            warning("Colnames of some weight matrices are empty. Be aware that without variable names there is no check if variable order matches with the one in datasets provided. That can lead to wrong estimation.")
            W <- lapply(W, jewel:::removeDiagonal)
            W <- do.call(rbind, W)
        }
        else if (sum(sapply(1:K, function(k) sum(colnames(W[[k]]) %in% 
                                                 vars) == p)) != K | sum(sapply(1:K, function(k) sum(vars %in% 
                                                                                                     colnames(W[[k]])) == p)) != K) {
            stop("There are some colnames that do not match between provided X and W. Please check.")
        }
        else {
            W <- lapply(W, function(x) x[vars, vars])
            W <- lapply(W, jewel:::removeDiagonal)
            W <- do.call(rbind, W)
        }
    }
    if (is.null(lambda2)) {
        lambda2 <- lambda1 * 1.4
    }
    eps <- 2.220446e-16
    numIter <- 1
    check_conv <- 10000
    if (verbose) 
        message("1/3 Initialization completed. Starting iterations...")
    while (numIter <= maxIter && check_conv > tol) {
        numIter <- numIter + 1
        if (verbose) 
            message("jewel: iteration number ", numIter - 1)
        Theta_old <- Theta
        order1 <- sample(1:(p - 1), size = (p - 1))
        for (j in order1) {
            Act <- (j - 1) + which(Active[j:(p - 1), j], arr.ind = TRUE)
            if (length(Act) == 0) 
                break
            jminus <- setdiff(1:p, j)
            for (a in 1:length(Act)) {
                aminus <- setdiff(1:p, Act[a] + 1)
                za <- c(NA, K)
                zb <- c(NA, K)
                za <- sapply(index, function(k) (1/n_k[k]) * 
                                 Xl[nindex == k, jminus[Act[a]]] %*% R_common[nindex == 
                                                                                  k, j] + Xi[(p - 1) * (k - 1) + Act[a], j])
                zb <- sapply(index, function(k) (1/n_k[k]) * 
                                 Xl[nindex == k, aminus[j]] %*% R_common[nindex == 
                                                                             k, Act[a] + 1] + Xi[(p - 1) * (k - 1) + j, 
                                                                                                 Act[a] + 1])
                z <- c(za, zb)
                W_avr <- mean(sapply(index, function(k) W[(p - 
                                                               1) * (k - 1) + j, Act[a] + 1]))
                thrld <- 1 - lambda1 * sqrt(2 * K) * W_avr/(sqrt(sum(z^2)) + 
                                                                eps)
                if (thrld <= 0) {
                    z <- z * 0
                    Active[j, Act[a] + 1] <- FALSE
                    Active[Act[a], j] <- FALSE
                }
                else {
                    z <- z * thrld
                }
                za <- z[1:K]
                zb <- z[(K + 1):(2 * K)]
                R_common[, j] <- unlist(lapply(index, function(k) R_common[nindex == 
                                                                               k, j] - Xl[nindex == k, jminus[Act[a]]] * (za[k] - 
                                                                                                                              Xi[(p - 1) * (k - 1) + Act[a], j])))
                R_common[, Act[a] + 1] <- unlist(lapply(index, 
                                                        function(k) R_common[nindex == k, Act[a] + 
                                                                                 1] - Xl[nindex == k, aminus[j]] * (zb[k] - 
                                                                                                                        Xi[(p - 1) * (k - 1) + j, Act[a] + 1])))
                for (k in index) {
                    Xi[(p - 1) * (k - 1) + Act[a], j] <- za[k]
                    Xi[(p - 1) * (k - 1) + j, Act[a] + 1] <- zb[k]
                }
            }
        }
        {
            for (i in 1:K) {
                Gamma_list[[i]] <- jewel:::addZeroDiagonal(Gamma[((p - 
                                                                       1) * (i - 1) + 1):((p - 1) * i), ])
                Xi_list[[i]] <- jewel:::addZeroDiagonal(Xi[((p - 1) * 
                                                                (i - 1) + 1):((p - 1) * i), ])
            }
            update <- mapply(function(x, y, z) x - eigenMapMatMult(x, 
                                                                   y) - eigenMapMatMult(x, z), X, Gamma_list, Xi_list, 
                             SIMPLIFY = FALSE)
            R_specific <- do.call(rbind, update)
        }
        order2 <- sample(1:(p - 1), size = (p - 1))
        for (k in 1:K) {
            for (j in order2) {
                Act <- (p - 1) * (k - 1) + (j - 1) + which(Active_K_long[((p - 
                                                                               1) * (k - 1) + j):((p - 1) * k), j], arr.ind = TRUE)
                if (length(Act) == 0) 
                    break
                for (a in 1:length(Act)) {
                    i <- Act[a] + k - p * (k - 1)
                    za <- (1/n_k[k]) * Xl[nindex == k, i] %*% R_specific[nindex == 
                                                                             k, j] + Gamma[Act[a], j]
                    zb <- (1/n_k[k]) * Xl[nindex == k, j] %*% R_specific[nindex == 
                                                                             k, i] + Gamma[(p - 1) * (k - 1) + j, i]
                    z <- c(za, zb)
                    thrld <- 1 - lambda2 * sqrt(2) * W[(p - 1) * 
                                                           (k - 1) + j, i]/(sqrt(sum(z^2)) + eps)
                    if (thrld <= 0) {
                        z <- z * 0
                        Active_K_long[Act[a], j] <- FALSE
                        Active_K_long[(p - 1) * (k - 1) + j, i] <- FALSE
                    }
                    else {
                        z <- z * thrld
                    }
                    za <- z[1]
                    zb <- z[2]
                    R_specific[nindex == k, j] <- R_specific[nindex == 
                                                                 k, j] - Xl[nindex == k, i] * (za - Gamma[Act[a], 
                                                                                                          j])
                    R_specific[nindex == k, i] <- R_specific[nindex == 
                                                                 k, i] - Xl[nindex == k, j] * (zb - Gamma[(p - 
                                                                                                               1) * (k - 1) + j, i])
                    Gamma[Act[a], j] <- za
                    Gamma[(p - 1) * (k - 1) + j, i] <- zb
                }
            }
        }
        {
            for (i in 1:K) {
                Gamma_list[[i]] <- jewel:::addZeroDiagonal(Gamma[((p - 
                                                                       1) * (i - 1) + 1):((p - 1) * i), ])
                Xi_list[[i]] <- jewel:::addZeroDiagonal(Xi[((p - 1) * 
                                                                (i - 1) + 1):((p - 1) * i), ])
            }
            update <- mapply(function(x, y, z) x - eigenMapMatMult(x, 
                                                                   y) - eigenMapMatMult(x, z), X, Gamma_list, Xi_list, 
                             SIMPLIFY = FALSE)
            R_common <- do.call(rbind, update)
        }
        Theta <- Xi + Gamma
        check_conv <- sum(abs(Theta - Theta_old))/(sum(abs(Theta_old)) + 
                                                       eps)
    }
    if (verbose) 
        message("jewel: total number of iterations is ", numIter - 
                    1, " and the error is ", round(check_conv, digits = 5))
    if (verbose) 
        message("2/3 Iterations completed. Assembling the output...")
    BIC <- sum(sapply(index, function(c) n_k[c] * sum(apply(R_common[nindex == 
                                                                         c, ], 2, function(y) log(sum(y^2)))))) + sum(sapply(n_k, 
                                                                                                                             log) * sapply(index, function(c) sum(Active_K_long[nindex == 
                                                                                                                                                                                    c, ])/2))
    Theta_list <- rep(list(NA), K)
    A_list <- rep(list(NA), K)
    names(A_list) <- names(X)
    for (i in 1:K) {
        Theta_list[[i]] <- jewel:::addZeroDiagonal(Theta[((p - 1) * (i - 
                                                                         1) + 1):((p - 1) * i), ])
        A_list[[i]] <- (Theta_list[[i]] != 0)
        colnames(A_list[[i]]) <- rownames(A_list[[i]]) <- colnames(X[[i]])
    }
    A <- Reduce("+", A_list)
    A <- (A == K)
    if (verbose) 
        message("3/3 Completed.")
    return(list(G_list = A_list, CommonG = A, Theta = Theta_list, BIC = BIC))
}

