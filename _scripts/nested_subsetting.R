library(qtl2)
library(doParallel)
library(foreach)
library(tidyverse)

get_sample_list <- function(data) {
  data <- as.data.frame(data)
  samp_350 <- slice_sample(data, n = 350, replace = FALSE)
  samp_300 <- slice_sample(samp_350, n = 300, replace = FALSE)
  samp_250 <- slice_sample(samp_300, n = 250, replace = FALSE)
  samp_200 <- slice_sample(samp_250, n = 200, replace = FALSE)
  samp_150 <- slice_sample(samp_200, n = 150, replace = FALSE)
  samp_100 <- slice_sample(samp_150, n = 100, replace = FALSE)
  samp_50 <- slice_sample(samp_100, n = 50, replace = FALSE)
  
  sample_list <- tibble::lst(samp_350, samp_300, samp_250, samp_200, samp_150, samp_100, samp_50)
  return(sample_list)
}

# need to have addcovar, K, and genoprobs (the overall ones for all samples) already defined

nested_subsetting <- function(sample_list) {
  
  results_list <- list() # this will store the results for all subsets in this nest
  
  for (i in seq_along(sample_list)) {
    
    # get mouse ids from the rows of the subsetted data and filter addcovar to only have these rows
    match_ids <- rownames(sample_list[[i]])
    addcovar_sub <- addcovar[match_ids, ,drop=FALSE]
    
    # now perform genome scan to get LOD scores
    # uses scan1 to calculate lod scores
    sub_lods <- scan1(genoprobs = genoprobs,
                         pheno = sample_list[[i]],
                         kinship = K,
                         addcovar = addcovar_sub,
                         cores = 10) 
    
    # print("lods found")
    
    # perform the permutations with scan1perm
    sub_perms <- scan1perm(genoprobs = genoprobs,
                          pheno = sample_list[[i]],
                          addcovar = addcovar_sub,
                          n_perm = 1000, 
                          cores = 10) 
    # print("perms done")
    
    # get the 95% significance threshold from the permutation analysis
    sub_thr <- summary(object = sub_perms,
                   alpha = 0.05)
    
    # print("Threshold found")
    
    # identify the significant peaks with find_peaks
    sub_peaks <- find_peaks(scan1_output = sub_lods, 
                            map          = map, 
                            threshold = sub_thr, 
                            prob         = 0.95)
    
    # print("peaks found")
      
    #use scan1blup to calculate the blups
    sub_blups <- scan1blup(genoprobs = genoprobs[,10],
                          pheno = sample_list[[i]],
                          kinship = K[[10]],
                          addcovar = addcovar_sub,
                          cores = 10) 
    
    # print("blups found")
    
    # subset the kinship matrix to match 
    sub_kin <- kin_overall[rownames(sample_list[[i]]), rownames(sample_list[[i]])]
    
    # print("kinship matrix subsetted")
    # print(dim(sub_kin))
    
    # calculate heritability to be used for % variance explained
    sub_herit <- est_herit(pheno = sample_list[[i]],
                           kinship = sub_kin, 
                           addcovar = addcovar_sub,
                           cores = 10)
   
    # print("heritability found")
    
    # Here add all of this to a big results list and return it
    
    # store blups in a list (NULL or not) and add sample size
    sub_blups <- list(
      blups = sub_blups,
      sample_size = nrow(sample_list[[i]]))
    
    # add sample size column to thr data frame
    sub_thr_df <- as.data.frame(sub_thr)
    sub_thr_df$sample_size <- nrow(sample_list[[i]])
    
    # add columns to the peaks that indicate if peak was found or not so we can sum them
    if (nrow(sub_peaks) > 0) { # check that there are peaks and if so, add sample size and peak_found and make a list
      sub_peaks$marker <- rownames(sub_lods)[which(sub_lods == max(sub_lods))]
      sub_peaks$max_lod <- sub_lods[which(sub_lods == max(sub_lods))]
      sub_peaks = list(
        peaks = sub_peaks,
        sample_size = nrow(sample_list[[i]]),
        peak_found = TRUE)
    } else {
      placeholder <- data.frame( 
        lodindex = NA,
        lodcolumn = NA,
        chr = NA,
        pos = NA,
        lod = sub_lods[which(rownames(sub_lods) == "UNC18956869")],
        max_lod = sub_lods[which(sub_lods == max(sub_lods))],
        ci_lo = NA,
        ci_hi = NA,
        marker = rownames(sub_lods)[which(sub_lods == max(sub_lods))],
        stringsAsFactors = FALSE
      )
      sub_peaks <- list(
        peaks = placeholder,
        sample_size = nrow(sample_list[[i]]),
        peak_found = FALSE
      )
    }
    
    # add sample size column to herit data frame
    sub_herit_df <- as.data.frame(sub_herit)
    sub_herit_df$sample_size <- nrow(sample_list[[i]])
    
    results_list[[i]] <- list(
      blups = sub_blups,
      thr = sub_thr_df,
      peaks = sub_peaks,
      herit = sub_herit)
    
    # print(paste("Result", i))
  }
  
  names(results_list) <- sapply(sample_list, nrow)
  
  return(results_list)
}

# read in the data
LiverLipids_Phenotypes_V11 <- readRDS("/projects/munger-lab/projects/DO_Attie500/jared_power/_data/LiverLipids_Phenotypes_V11.rds")
load("/projects/munger-lab/projects/DO_Attie500/data/rdata/attie.core.GRCm39.v1.Rdata")
full_data_peaks <- readRDS("/projects/munger-lab/projects/DO_Attie500/jared_power/_data/full_dataset_peaks.rds")

#subset rank Z data column
rz_data <- LiverLipids_Phenotypes_V11[[4]][[3]]
rz_sing_data <- as.data.frame(rz_data[, "UNK_9.698_1251.79456_minus", drop=FALSE])

missing_rows <- setdiff( rownames(rz_sing_data), dimnames(genoprobs$`1`)[[1]])
rz_sing_data <- rz_sing_data[!rownames(rz_sing_data) %in% missing_rows, , drop=FALSE]

#subset the covariate info from the annotation data
samp_annot <- LiverLipids_Phenotypes_V11[[2]]

# make the covariates factors
samp_annot$Sex <- factor(samp_annot$Sex)
samp_annot$Wave <- factor(samp_annot$Wave)
samp_annot$Batch <- factor(samp_annot$Batch)

# build the covariate matrix
addcovar <- model.matrix(~Sex + Wave + Batch, data = samp_annot)[, -1, drop=FALSE]

# addocovar only includes rows which don't have NA values. Filter these out so the rows of 
# addcovar match with the rows in the rz_sing_data and store the mouse ids
mouse_ids <- samp_annot$MouseID[!is.na(samp_annot$Batch)]

# replace the rownames of addcovar with the mouse ids
rownames(addcovar) <- mouse_ids # this addcovar is now the one subsetted in the function

# calculate the overall kinship matrix used in the function for heritability
kin_overall <- calc_kinship(probs <- genoprobs, 
                            type = "overall")

# generate 10 random seeds
num_iter <- 10
seed_list <- sample(1000:10000, num_iter)
samp_lists <- list() # initialize empty sample list

# use get sample list function to generate nested samples from the random seeds and store in big list
for (i in 1:num_iter) {
  set.seed(seed_list[i])
  
  samp_list <- get_sample_list(rz_sing_data)
  
  samp_lists[[i]] <- samp_list
}

# save the seed list in case we need it later
write(seed_list, file = "/projects/munger-lab/projects/DO_Attie500/jared_power/_data/seed_list_rep10.txt", ncolumns = 1)

print(Sys.time())

# do the scans on the nested subsets in parallel
doParallel::registerDoParallel(cores = 70) 
nested_iteration_data <- foreach::foreach(samp_list = samp_lists) %dopar% {
  nested_subsetting(sample_list = samp_list)
}
doParallel::stopImplicitCluster()

print(Sys.time())

# nested_iteration_data will be a list where each element 
# is an iteration that is a list of data from all sample sizes

saveRDS(nested_iteration_data, file = "/projects/munger-lab/projects/DO_Attie500/jared_power/_data/nested_iterations_rep10.rds")
