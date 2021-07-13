cwlVersion: v1.2
class: CommandLineTool
label: Ancestry_MAF_estimate_global_proportion
$namespaces:
  sbg: https://sevenbridges.com

requirements:
- class: ShellCommandRequirement
- class: ResourceRequirement
  coresMin: $(inputs.cpu * inputs.cpu_multiplier)
  ramMin: 8000
- class: DockerRequirement
  dockerPull: |-
    images.sb.biodatacatalyst.nhlbi.nih.gov/einatgranot/ancestral_maf_admixed_population:20210507
- class: InitialWorkDirRequirement
  listing:
  - entryname: freq_ancestry_global.R
    writable: true
    entry: |2+

      #Calculate SOL ancetries MAFs based on Tamar's function using global proportion ancestries (without prior knowledge on Ancestral MAFs)
      # scan if packages are installed, if not install

      library(GWASTools)
      library(gdsfmt)
      library(SeqArray)
      library(SeqVarTools)
      library(dplyr)
      library(stringr)
      library(tictoc)

      source("cwl_input.R")

      plist<-read.csv(plist_path, header=T)

      #get command line args
      args = commandArgs(trailingOnly=TRUE)
      myindex <- args[1]

      #j is chromosomes
      j  <- plist[myindex,1]

      #b is blocks
      b <- plist[myindex,2]


      ########################Define functions for allele frequency caculation in admixed populations#################################
      source('decompose_two_alleles_one_person.R', chdir = TRUE)
      source('generate_smoothing_observations.R', chdir = TRUE)
      source('estimate_frequencies.R', chdir = TRUE)
      source('prep_dat_for_binomial_likelihood.R', chdir = TRUE)
      source('dbinom_approx.R', chdir = TRUE)
      source('estimate_frequencies_search_boundary.R', chdir = TRUE)
      source('estimate_frequencies_w_known_freqs.R', chdir = TRUE)



      #Read SOL global proportion ancestry
      tic("Read SOL global proportion ancestry"); write("Read SOL global proportion ancestry", stderr()); print("Read SOL global proportion ancestry")

      anc <- read.table(global_file_path, header=TRUE)
      annot <- read.csv(annotation_file_path)
      toc(log = TRUE)

      tic("merge 1"); write("merge", stderr())
      anc <- merge(annot[,c("SUBJECT_ID","genetic_consent")], anc[,c("SOL_ID","Africa","Europe","America")], by.x = "SUBJECT_ID", by.y = "SOL_ID")#10547
      anc <- anc[which(anc$genetic_consent == 1),]#10490 
      anc<-anc[,-c(2)]#10490

      toc(log = TRUE)

      tic("Read list of unrelated individuals"); write("Read list of unrelated individuals", stderr())
      sol_unrel<-read.table(sol_unrel_path)[1]#10223
      names(sol_unrel)[1]<-"SUBJECT_ID"
      #Create a list of the final SOL participants available for MAF calculations based on global proportions.
      write("merge 2", stderr())
      anc<-merge(anc, sol_unrel, by="SUBJECT_ID", all=F)#8966

      #Read gds imputed genotypes file for chromosome j
      write("Read gds imputed genotypes file for chromosome", stderr())
      gds <- seqOpen(chromosome_full_path, FALSE)

      #Subset and sort global ancestry proportions to the sample.id in the gds file
      sample.id <- seqGetData(gds, "sample.id")#11928
      sample.id_logical<-sample.id %in% anc$SUBJECT_ID#8933
      sample.id<-as.data.frame(sample.id[sample.id_logical])#8933
      names(sample.id)[1]<-"SUBJECT_ID"
      write("inner join", stderr())

      var<-inner_join(sample.id, anc, by="SUBJECT_ID")#8933
      prop_mat <- as.matrix(var[,c("Africa", "Europe", "America")])
      #Normalize so it will sum up to 1 (since the rounding sometimes causes it to sum to above or below 1)
      write("Normalize", stderr())

      prop_mat <- prop_mat/rowSums(prop_mat)

      snp.rs.id <- seqGetData(gds, "variant.id")#1166746
      nvar <- length(snp.rs.id)
      chr.id <-seqGetData(gds, "chromosome")#1166746
      position.id <-seqGetData(gds, "position")#1166746
      write("alleles", stderr())
      alleles<-seqGetData(gds, "allele")#1166746
      allele_a<-str_split_fixed(alleles, ",",2)[,1]#1166746
      allele_b <-str_split_fixed(alleles, ",",2)[,2]#1166746

      #Partition the gds genotyped data to segments with 3000 variants

      # compute the number of blocks
      b.size <- 3000

      # set filters by blocks and read data
      write("set filters", stderr())

      b.start.inds <- seq(1, nvar, by = b.size)
      b.end.inds <- seq(b.size, nvar, by = b.size)
      if (length(b.start.inds) > length(b.end.inds)) b.end.inds <- c(b.end.inds, nvar)
      b.size_seg<-c(rep(b.size,length(b.start.inds)-1),b.end.inds[length(b.start.inds)]-b.start.inds[length(b.start.inds)]+1)


      # set sample and variant filters
      write("set sample and vairant filters", stderr())

      seqSetFilter(gds, c(b.start.inds[b]:b.end.inds[b]), sample.sel=sample.id_logical, verbose=TRUE)

      #Check compatibility of order of individuals in the gds file and in prop_mat
      #sample.id <- seqGetData(gds, "sample.id")
      #table(sample.id==var$SUBJECT_ID)#8933

      write("imputed Dosage", stderr())
      gen<-as.data.frame(imputedDosage(gds, dosage.field="DS", use.names=TRUE))
      closefn.gds(gds)
      ancestry_mafs <-matrix(NA,nrow=ncol(gen),ncol=12)


      for (i in c(1:ncol(gen))){
        allele_counts<-as.numeric(gen[,i])
        est_freq<-estimate_frequencies_dynamic_boundary (allele_counts, prop_mat, confidence = 0.95, 
                                                        frequency_boundary_grid = c(0.00001, 0.01),
                                                        use_smoothing_data = TRUE,
                                                        chromosome_x = FALSE,
                                                        sex = NULL, 
                                                        male_label = "M", 
                                                        mac_filter = 5)
        est_freq1<-est_freq[["res"]]
        est_freq1<-as.vector(t(est_freq1))
        est_freq1<-est_freq1[-c(1,5,9)]
        nll<-est_freq[["nll"]]
        boundary<-est_freq[["boundary"]]
        flag<-est_freq[["flag"]]
        ancestry_mafs[i,]<-c(est_freq1,nll,boundary,flag)
        }


      ancestry_mafs<- data.frame(ancestry_mafs)

      names(ancestry_mafs)<-c("Africa_est", "Africa_low_CI", "Africa_high_CI","Europe_est", "Europe_low_CI", "Europe_high_CI","America_est", "America_low_CI", "America_high_CI", "nll","boundary", "flag")

      ancestry_mafs<- data.frame(ancestry_mafs)
      ancestry_mafs$CHR<-chr.id[b.start.inds[b]:b.end.inds[b]]
      ancestry_mafs$POS<-position.id[b.start.inds[b]:b.end.inds[b]]
      ancestry_mafs$allele_a<-allele_a[b.start.inds[b]:b.end.inds[b]]
      ancestry_mafs$allele_b<-allele_b[b.start.inds[b]:b.end.inds[b]]


      write.csv(ancestry_mafs, paste0("chr",j,"_",b,"_gds_ancestry_mafs_global_121120.csv"), quote=F, row.names=F)




  - entryname: cwl_input.R
    writable: false
    entry: |-
      plist_path="$(inputs.plist.path)" 
      global_file_path="$(inputs.global_ancestry.path)"
      annotation_file_path="$(inputs.annotation.path)"
      chromosome_full_path="$(inputs.chromosome_gds.path)"
      chromosome_number=$(inputs.chromosome_gds.path.split("chr")[1].split("_")[0])
      sol_unrel_path="$(inputs.sol_unrel.path)"
      cpu=$(inputs.cpu)
  - entryname: decompose_two_alleles_one_person.R
    writable: false
    entry: |-
      
      # Take an allele count (values 0, 1, or 2) and return a vector of allele counts
      # in two chromosomes (values 0 or 1)
      
      decompose_two_alleles_one_person <- function(allele_count){
        
        # check that only valid allele counts are present
        stopifnot(is.element(allele_count, c(0,1,2)))
        
        if (allele_count == 0){
          allele_count_two <- c(0,0)
        }
        if (allele_count == 1){
          allele_count_two <- c(1,0)
        }
        if (allele_count == 2){
          allele_count_two <- c(1,1)
        }
        
        return(c(decomposed_allele_count =  allele_count_two))
      }
  - entryname: generate_smoothing_observations.R
    writable: false
    entry: |
      # Generate artificial observations with alleles from each provided ancestries.
      # These observations are later added to the real data to protect from reaching boundary condition
      # and failure of the estimation function. 

      .generate_smoothing_observations <- function(ancestry_names){
        
        n_ancestry <- length(ancestry_names)
        
        # construct a matrix with two rows per ancestry, corresponding
        # to two observations with 100% global proportions of that ancestry
        prop_mat <- matrix(0, nrow =n_ancestry*2, ncol = n_ancestry)
        colnames(prop_mat) <- ancestry_names
        
        for (i in 1:n_ancestry){
          row_inds <- (i-1)*2+c(1,2)
          prop_mat[row_inds,i] <- 1
        }
         
        # for each ancestry, we set one observation with one variant allele
        # and a second observations with non variant alleles. 
        allele_counts <- rep(c(1,0), n_ancestry)
        
        return(list(simulated_prop_mat = prop_mat, simulated_allele_count = allele_counts))
         
      }
  - entryname: estimate_frequencies.R
    writable: false
    entry: |
      
      # Function to estimate ancestry-specific allele frequencies for a given variant.
      # takes a vector of allele_counts allele counts for a given variant for n people 
      # and an n x d matrix prop_mat giving, for each of n people, d ancestral ancestry proportion
      # at the locus (could be genome-wide) for d ancestries. 
      estimate_frequencies <- function(allele_counts, prop_mat, confidence = 0.95, 
                                       low_freq_bound = 0.001, high_freq_bound = 0.999,
                                       use_smoothing_data = FALSE,
                                       chromosome_x = FALSE,
                                       sex = NULL, 
                                       male_label = "M", 
                                       mac_filter = 5){
        stopifnot(length(allele_counts) == nrow(prop_mat))
        
        prop_mat <- as.matrix(prop_mat)
        
        # check for NAs, if there are observations with missing values, remove them.
        inds_na_alleles <- which(is.na(allele_counts))
        inds_na_prop <- which(apply(prop_mat, 1, function(x) sum(is.na(x))) > 0)
        inds_na <- c(inds_na_alleles, inds_na_prop)
        if (length(inds_na) >0){
          message(paste(length(inds_na), "observations with missing values, removing them..."))
          allele_counts <- allele_counts[-inds_na]
          prop_mat <- prop_mat[-inds_na,]
          sex <- sex[-inds_na]
        }
        
        prep_dat <- .prep_dat_for_binomial_likelihood(allele_counts, prop_mat,
                                                    chromosome_x = chromosome_x,
                                                    sex = sex, 
                                                    male_label = male_label)
        
        prop_mat <- prep_dat$prop_mat
        allele_counts <- prep_dat$allele_counts
        max_counts <- prep_dat$max_counts
        
        # check if the number of minor alleles is higher than the mac_filter,
        # stop if MAC is too low.
        stopifnot(min(sum(allele_counts), sum(max_counts) - sum(allele_counts)) > mac_filter)
        
        # add made-up data to avoid boundaries of the frequency parameter space   
        if (use_smoothing_data){
          smoothing_data <- .generate_smoothing_observations(colnames(prop_mat))
          prop_mat <- rbind(prop_mat, smoothing_data$simulated_prop_mat)
          allele_counts <- c(allele_counts, smoothing_data$simulated_allele_count)
          max_counts <- c(max_counts, rep(1, nrow(smoothing_data$simulated_prop_mat)))
        }
        
        
       
        return(.estimate_frequencies_after_prep(allele_counts, 
                                                prop_mat, 
                                                max_counts,
                                                low_freq_bound,
                                                high_freq_bound,
                                                confidence))
      }
      
      
      ## separate function that is called after checks and preparations were done
      .estimate_frequencies_after_prep <- function(allele_counts,
                                                   prop_mat, 
                                                   max_counts,
                                                   low_freq_bound,
                                                   high_freq_bound,
                                                   confidence){
        ## compute the negative log likelihood function
        nll <- function(freqs){
          allele_probs <- as.numeric(prop_mat %*% freqs)
          nll_by_obs <- log(dbinom_approx(allele_counts, max_counts, allele_probs))
          return(-sum(nll_by_obs))		
        }
        
        
        ## optimize to estimate frequencies
        fit <- optim(par = rep(1/ncol(prop_mat), ncol(prop_mat)), fn = nll, hessian = TRUE, 
                     lower = rep(low_freq_bound ,ncol(prop_mat)), upper = rep(high_freq_bound, ncol(prop_mat)), 
                     method = "L-BFGS-B")
        
        estimated_freqs <- fit$par
        hessian <- fit$hessian
        ses <- sqrt(diag(solve(hessian)))
        res <- data.frame(ancestry = colnames(prop_mat), 
                          estimated_freq = estimated_freqs, 
                          low_CI = estimated_freqs - ses*sqrt(qchisq(confidence, 1)), 
                          high_CI= estimated_freqs + ses*sqrt(qchisq(confidence, 1)))
        
        # return the estimated frequencies, and the negative log likelihood. 
        return(list(res = res, nll = nll(res$estimated_freq)))
        
      }
  - entryname: prep_dat_for_binomial_likelihood.R
    writable: false
    entry: |-
      .prep_dat_for_binomial_likelihood <- function(allele_counts, prop_mat,  
                                                  chromosome_x = FALSE,
                                                  sex = NULL, 
                                                  male_label = "M"){
        max_counts <- rep(NA, length(allele_counts))
        # If alleles are from chromosome X, need to be separately handled for males and females.
        if (chromosome_x){
          stopifnot(length(allele_counts) == length(sex))
          male_counts <- allele_counts[which(sex == male_label)]
          male_props <- prop_mat[which(sex == male_label),]
          
          # check that male counts are not higher than one (imputed data may have fractions)
          # if some dosages are higher than 1, divided all dosages/counts by 2. 
          max_male_dosage <- max(male_counts)
          if (max_male_dosage > 1){
            message("some male allele counts/dosages are larger than 1, dividing all counts by two...")
            male_counts <- male_counts/2
          }
          
          female_counts <- allele_counts[which(sex != male_label)]
          female_props <- prop_mat[which(sex != male_label),]
          
          # merge male and females:
          prop_mat <- rbind(male_props, female_props)
          allele_counts <- c(male_counts, female_counts)
          max_counts <- c(rep(1, length(male_counts)), rep(2, length(female_counts)))
        } else{  # not chromosome X
          
          # just add max_count
          max_counts <- rep(2, nrow(prop_mat))
          
        }
        
        return(list(allele_counts = allele_counts, 
                    prop_mat = prop_mat, 
                    max_counts = max_counts))
        
      }
  - entryname: dbinom_approx.R
    writable: false
    entry: |
      # function that returns binomial probability mass function, 
      # with linear interpolation to the case where we have dosage data.
      
      dbinom_approx <- function(x, size, prob){
        
        # add checks to account for numerical error
        if (sum(prob > 1 + 1e-10 | prob < -1e-10) > 0)  
          stop("probability higher than 1 or lower than 0")
        prob <- ifelse(prob > 1, 1, prob)
        prob <- ifelse(prob <0 , 0, prob)
        
        arg2 <- prob^(x)*(1-prob)^(size-x)
        
        # now compute arg1
        floor_nums <- floor(x) 
        ceiling_nums <- ceiling(x)
        
        floor_vals <- ifelse(floor_nums == 1, 2, 1)
        ceiling_vals <- ifelse(ceiling_nums == 1, 2, 1)
        
        arg1 <- floor_vals + (ceiling_vals - floor_vals)*(x-floor_nums)
        
        return(arg1*arg2)
      }
  - entryname: estimate_frequencies_w_known_freqs.R
    writable: false
    entry: |
      ##
      # Function to estimate ancestry-specific allele frequencies for a given variant,
      # when some (one or more) ancestral allele frequencies are known and provided.
      # takes a vector of allele_counts allele counts for a given variant for n people 
      # and an n x d matrix matrix prop_mat giving, for each of n people, d ancestral ancestry proportion
      # at the locus (could be genome-wide) for d ancestries. 
      # known_freq need to be a named vector, specifying known frequencies. The names of elements 
      # in known_freq are some of the column names of prop_mat. 
      estimate_frequencies_w_known_freqs <- function(allele_counts, prop_mat, known_freqs, 
                                                     confidence = 0.95, 
                                                  low_freq_bound = 0.001, high_freq_bound = 0.999,
                                                  use_smoothing_data = FALSE,
                                                  chromosome_x = FALSE,
                                                  sex = NULL, 
                                                  male_label = "M",
                                                  mac_filter = 5){
        
        stopifnot(length(allele_counts) == nrow(prop_mat))
        stopifnot(all(is.element(names(known_freqs), colnames(prop_mat))))
        stopifnot(length(known_freqs) > 0 & length(known_freqs) < ncol(prop_mat))
        
        prop_mat <- as.matrix(prop_mat)
        
        # check for NAs, if there are observations with missing values, remove them.
        inds_na_alleles <- which(is.na(allele_counts))
        inds_na_prop <- which(apply(prop_mat, 1, function(x) sum(is.na(x))) > 0)
        inds_na <- c(inds_na_alleles, inds_na_prop)
        if (length(inds_na) >0){
          message(paste(length(inds_na), "observations with missing values, removing them..."))
          allele_counts <- allele_counts[-inds_na]
          prop_mat <- prop_mat[-inds_na,]
          sex <- sex[-inds_na]
        }
        
        # re-order the column so that the ones with known frequencies are at the beginning:
        n_known <- length(known_freqs)
        n_unknown <- ncol(prop_mat) - n_known
        inds_known <- which(is.element(colnames(prop_mat), names(known_freqs)))
        inds_unknown <- setdiff(1:ncol(prop_mat), inds_known)
        prop_mat <- prop_mat[,c(inds_known, inds_unknown)]
        
        prep_dat <- .prep_dat_for_binomial_likelihood(allele_counts, prop_mat,
                                                    chromosome_x = chromosome_x,
                                                    sex = sex, 
                                                    male_label = male_label)
        
        prop_mat <- prep_dat$prop_mat
        allele_counts <- prep_dat$allele_counts
        max_counts <- prep_dat$max_counts
        
        # check if the number of minor alleles is higher than the mac_filter,
        # stop if MAC is too low.
        stopifnot(min(sum(allele_counts), sum(max_counts) - sum(allele_counts)) > mac_filter)

        # add made-up data to avoid boundaries of the frequency parameter space 
        if (use_smoothing_data){
          smoothing_data <- .generate_smoothing_observations(colnames(prop_mat))
          # add only the rows corresponding to the ancestries with unknown frequencies
          inds_keep <- (n_known*2 + 1):((n_known + n_unknown)*2)
          prop_mat <- rbind(prop_mat, smoothing_data$simulated_prop_mat[inds_keep,])
          allele_counts <- c(allele_counts, smoothing_data$simulated_allele_count[inds_keep])
          max_counts <- c(max_counts, rep(1, length(inds_keep)))
        }
        
        
        return(.estimate_frequencies_w_known_freq_after_prep(allele_counts, 
                                                             prop_mat, 
                                                             max_counts,
                                                             known_freqs, 
                                                             n_unknown,
                                                             n_known,
                                                             low_freq_bound, 
                                                             high_freq_bound,
                                                             confidence))
        
      }



      .estimate_frequencies_w_known_freq_after_prep <- function(allele_counts, 
                                                                prop_mat, 
                                                                max_counts,
                                                                known_freqs, 
                                                                n_unknown,
                                                                n_known,
                                                                low_freq_bound, 
                                                                high_freq_bound,
                                                                confidence){
        
        
        known_sums <- prop_mat[,names(known_freqs), drop = FALSE] %*% known_freqs
        

        ## compute the negative log likelihood function
        nll <- function(unknown_freqs){
          allele_probs <- as.numeric(known_sums + prop_mat[,c(n_known+1):ncol(prop_mat), drop = FALSE] %*% unknown_freqs)
          nll_by_obs <- log(dbinom_approx(allele_counts, max_counts, allele_probs))
          return(-sum(nll_by_obs))		
        }
        
        ## optimize to obtain allele frequency estimates
        
        if (n_unknown == 1){
          fit <- optim(par = rep(0.5, n_unknown), fn = nll, hessian = TRUE,
                       method = "Brent", lower = low_freq_bound, upper = high_freq_bound)
        } else{
          fit <- optim(par = rep(0.5, n_unknown), fn = nll, hessian = TRUE,
                       lower = rep(low_freq_bound, n_unknown), 
                       upper = rep(high_freq_bound, n_unknown), 
                       method = "L-BFGS-B")
        }
        
        
        estimated_freqs <- fit$par
        hessian <- fit$hessian
        ses <- sqrt(diag(solve(hessian)))
        res <- data.frame(ancestry = colnames(prop_mat)[-c(1:n_known)], 
                          estimated_freq = estimated_freqs, 
                          low_CI = estimated_freqs - ses*sqrt(qchisq(confidence, 1)), 
                          high_CI= estimated_freqs + ses*sqrt(qchisq(confidence, 1)))
        
        # return the estimated frequencies, and the negative log likelihood. 
        return(list(res = res, nll = nll(res$estimated_freq)))
        
      }
  - entryname: estimate_frequencies_search_boundary.R
    writable: false
    entry: |
      # Function to estimate ancestry-specific allele frequencies for a given variant.
      # takes a vector of allele_counts allele counts for a given variant for n people 
      # and an n x d matrix prop_mat giving, for each of n people, d ancestral ancestry proportion
      # at the locus (could be genome-wide) for d ancestries. 
      # this is a wrapper function that changes the frequency boundaries and 
      # if frequencies for some ancestries are estimated at the set boundary conditions,
      # compares the likelihood to the likelihood when the frequencies are set at the 
      # boundaries of the parameter space (0 or 1, as needed).
      estimate_frequencies_dynamic_boundary <- function(allele_counts, prop_mat, confidence = 0.95, 
                                       frequency_boundary_grid = c(0.001, 0.01, 0.02, 0.05),
                                       use_smoothing_data = TRUE,
                                       chromosome_x = FALSE,
                                       sex = NULL, 
                                       male_label = "M", 
                                       mac_filter = 5){
       
        stopifnot(length(allele_counts) == nrow(prop_mat))
        
        prop_mat <- as.matrix(prop_mat)
        
        # check for NAs, if there are observations with missing values, remove them.
        inds_na_alleles <- which(is.na(allele_counts))
        inds_na_prop <- which(apply(prop_mat, 1, function(x) sum(is.na(x))) > 0)
        inds_na <- c(inds_na_alleles, inds_na_prop)
        if (length(inds_na) >0){
          message(paste(length(inds_na), "observations with missing values, removing them..."))
          allele_counts <- allele_counts[-inds_na]
          prop_mat <- prop_mat[-inds_na,]
          sex <- sex[-inds_na]
        }
        
        # prepare data.frame for returning null results if the function doesn't run
        # or does not converge
        return_val_not_run <- data.frame(ancestry = colnames(prop_mat),
                                         estimated_freq = NA,
                                         low_CI = NA,
                                         high_CI = NA)
        
         # check that MAC is higher than mac_filter
        prep_dat <- .prep_dat_for_binomial_likelihood(allele_counts, prop_mat,
                                                      chromosome_x = chromosome_x,
                                                      sex = sex, 
                                                      male_label = male_label)
        
        prop_mat <- prep_dat$prop_mat
        allele_counts <- prep_dat$allele_counts
        max_counts <- prep_dat$max_counts
        
      
        if (min(sum(allele_counts), sum(max_counts) - sum(allele_counts)) <= mac_filter){
          return(list(res = return_val_not_run, 
                      nll = NA,
                      boundary = NA, 
                      flag = "MAC filter"))
        }
        
        # add made-up data to avoid boundaries of the frequency parameter space   
        if (use_smoothing_data){
          smoothing_data <- .generate_smoothing_observations(colnames(prop_mat))
          prop_mat <- rbind(prop_mat, smoothing_data$simulated_prop_mat)
          allele_counts <- c(allele_counts, smoothing_data$simulated_allele_count)
          max_counts <- c(max_counts, rep(1, nrow(smoothing_data$simulated_prop_mat)))
        }
        
        
        # start loop to find boundaries in which the analysis converges:
        converged <- FALSE
        ind_grid <- 1
        while (!converged & ind_grid <= length(frequency_boundary_grid)){
          low_freq_bound <- frequency_boundary_grid[ind_grid]
          high_freq_bound <- 1 - frequency_boundary_grid[ind_grid]
          ## optimize to estimate frequencies
          
          cur_res <- tryCatch({.estimate_frequencies_after_prep(allele_counts = allele_counts, 
                                                                prop_mat = prop_mat, 
                                                                max_counts = max_counts,
                                                                low_freq_bound = low_freq_bound,
                                                                high_freq_bound = high_freq_bound,
                                                                confidence = confidence)},
                              error = function(error){
                                "not converged"
                              }
                            )
          if (!identical(cur_res, "not converged")){
            converged <- TRUE
          } else{
            ind_grid <- ind_grid + 1
          }
        }
        
        
        if (!converged) {
          return(list(res = return_val_not_run, 
                      nll = NA,
                      boundary = NA, 
                      flag = "Not converged"))
        } else{ # if converged, prepare potential return value
          potential_return_val_converged <- list(res = cur_res$res,
                                                 nll = cur_res$nll,
                                                 boundary = low_freq_bound,
                                                 flag = "converged")
        }
        
        # check if convergence is right at the boundary condition for some ancestries
        low_bound_inds <- which(cur_res$res$estimated_freq == low_freq_bound)
        high_bound_inds <- which(cur_res$res$estimated_freq == high_freq_bound)
        
        # if nothing was estimated at the boundary, return the computed frequencies.
        if (length(low_bound_inds) + length(high_bound_inds) == 0){
          return(potential_return_val_converged)
        }
        
        # otherwise, we continue to identify which frequencies were 
        # estimated at the boundaries and...
        # prepare a vector called set_freqs with in which the ancestries
        # where frequencies were estimated in the boundaries are set to the 
        # boundaries of the parameter space. Start with low and then high boundary.
          
        # if frequencies for some ancestries were estimated right at the lower bound:
        if (length(low_bound_inds) >0){
            set_freqs <- rep(0, length = length(low_bound_inds))
            names(set_freqs) <- cur_res$res$ancestry[low_bound_inds]
          }
        
        # if frequencies for some ancestries were estimated right at the high bound:
        if (length(high_bound_inds) >0 ){
            set_freqs_high <- rep(1, length = length(high_bound_inds))
            names(set_freqs_high) <- cur_res$res$ancestry[high_bound_inds]
            
            # should we augment a previous vector, and have a new set_freqs vector?
            if (length(low_bound_inds) > 0){
              set_freqs <- c(set_freqs, set_freqs_high)
            } else{
              set_freqs <- set_freqs_high
            }
          } # finished preparing a vector of set frequencies
       
        # first: if all frequencies were estimated at the boundary:
        if (length(set_freqs) == ncol(prop_mat)){
          allele_probs <- prop_mat[,names(set_freqs), drop = FALSE] %*% set_freqs
          like_by_obs <- dbinom_approx(allele_counts, max_counts, allele_probs)
          
          # check if we have an impossible value (probability zero)
          if (sum(like_by_obs == 0) > 0) return(potential_return_val_converged)
          
          nll <- -sum(log(like_by_obs))
          
          if (nll < potential_return_val_converged$nll){
            return_val <- return_val_not_run
            return_val[match(names(set_freqs), return_val$ancestry),"estimated_freq"] <- set_freqs
            return(list(res = return_val,
                        nll = nll,
                        boundary = potential_return_val_converged$boundary,
                        flag = "converged and evaluated monomorphic"))
          } else{
            return(potential_return_val_converged)
          }
        }
        
        # at this point: some frequencies but not all were estimated at the boundary
        boundary_res <- 
          tryCatch({.estimate_frequencies_w_known_freq_after_prep(allele_counts = allele_counts, 
                                                        prop_mat = prop_mat, 
                                                        max_counts = max_counts,
                                                        known_freqs= set_freqs, 
                                                        n_unknown = ncol(prop_mat) - length(set_freqs),
                                                        n_known = length(set_freqs),
                                                        low_freq_bound = low_freq_bound, 
                                                        high_freq_bound = high_freq_bound,
                                                        confidence = confidence)},
                   error = function(error){
                     "not converged"
                   }
          )
          # if this did not converge...
          if (identical(boundary_res, "not converged")){
           return(potential_return_val_converged)
          }
        
          # if it did converge, continue...
        
          if (boundary_res$nll < cur_res$nll) {
            # we need to return the new result, with the estimated 
            # frequency for the ancestry with boundary value set to that value
            return_val <- return_val_not_run
            return_val[match(boundary_res$res$ancestry, return_val$ancestry),] <- boundary_res$res
            return_val[match(names(set_freqs), return_val$ancestry),"estimated_freq"] <- set_freqs
            return(list(res = return_val,
                        nll = boundary_res$nll,
                        boundary = low_freq_bound,
                        flag = "converged some freqs evaluated monomorphic")
                        )
          } else{ 
            return(potential_return_val_converged)
          } 
      
      }
      
      
- class: InlineJavascriptRequirement
  expressionLib:
  - |2-

    var setMetadata = function(file, metadata) {
        if (!('metadata' in file)) {
            file['metadata'] = {}
        }
        for (var key in metadata) {
            file['metadata'][key] = metadata[key];
        }
        return file
    };
    var inheritMetadata = function(o1, o2) {
        var commonMetadata = {};
        if (!o2) {
            return o1;
        };
        if (!Array.isArray(o2)) {
            o2 = [o2]
        }
        for (var i = 0; i < o2.length; i++) {
            var example = o2[i]['metadata'];
            for (var key in example) {
                if (i == 0)
                    commonMetadata[key] = example[key];
                else {
                    if (!(commonMetadata[key] == example[key])) {
                        delete commonMetadata[key]
                    }
                }
            }
            for (var key in commonMetadata) {
                if (!(key in example)) {
                    delete commonMetadata[key]
                }
            }
        }
        if (!Array.isArray(o1)) {
            o1 = setMetadata(o1, commonMetadata)
            if (o1.secondaryFiles) {
                o1.secondaryFiles = inheritMetadata(o1.secondaryFiles, o2)
            }
        } else {
            for (var i = 0; i < o1.length; i++) {
                o1[i] = setMetadata(o1[i], commonMetadata)
                if (o1[i].secondaryFiles) {
                    o1[i].secondaryFiles = inheritMetadata(o1[i].secondaryFiles, o2)
                }
            }
        }
        return o1;
    };

inputs:
- id: global_ancestry
  label: Global  ancestry file
  doc: Proportion of global ancestries SOL
  type: File?
- id: chromosome_gds
  doc: 1-22 chromosomes
  type: File?
- id: annotation
  type: File?
- id: plist
  type: File?
  sbg:fileTypes: CSV
- id: myIndex
  type: int?
  inputBinding:
    position: 1
    shellQuote: false
- id: sol_unrel
  doc: List of third_degree unrelated SOL individuals
  type: File?
  sbg:fileTypes: .txt
- id: cpu
  type: int
- id: cpu_multiplier
  type: int

outputs:
- id: ancestry_csv
  type: File?
  outputBinding:
    glob: '*.csv'
    outputEval: $(inheritMetadata(self, inputs.chromosome_raw))
- id: workspace
  type: File?
  outputBinding:
    glob: '*.RData'
stdout: standard.out

baseCommand:
- Rscript freq_ancestry_global.R

hints:
- class: sbg:SaveLogs
  value: '*.R'
- class: sbg:SaveLogs
  value: standard.out
id: |-
  dave/build-ancestry-maf-admixed-population/ancestry-maf-estimate-global-proportion/20
sbg:appVersion:
- v1.2
sbg:content_hash: a3301f27024722caf8af8ccf8e04c9bec4c4a7866ef58191a82b47b5b7779648c
sbg:contributors:
- dave
- einatgranot
sbg:createdBy: dave
sbg:createdOn: 1620395842
sbg:id: |-
  dave/build-ancestry-maf-admixed-population/ancestry-maf-estimate-global-proportion/20
sbg:image_url:
sbg:latestRevision: 20
sbg:modifiedBy: einatgranot
sbg:modifiedOn: 1625667733
sbg:project: dave/build-ancestry-maf-admixed-population
sbg:projectName: 'BUILD: Ancestry MAF Admixed Population'
sbg:publisher: sbg
sbg:revision: 20
sbg:revisionNotes: Back to 'for' lopp due to minor discrepancies
sbg:revisionsInfo:
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1620395842
  sbg:revision: 0
  sbg:revisionNotes: |-
    Copy of einatgranot/2020-freq-estimation-latino-ancestries/ancestry-maf-estimate-global-proportion/100
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1620396434
  sbg:revision: 1
  sbg:revisionNotes: added tic toc
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1620408646
  sbg:revision: 2
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1620428633
  sbg:revision: 3
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1620428797
  sbg:revision: 4
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1620429115
  sbg:revision: 5
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1620441527
  sbg:revision: 6
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1620451495
  sbg:revision: 7
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1620452076
  sbg:revision: 8
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1620452593
  sbg:revision: 9
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1620452641
  sbg:revision: 10
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1620525916
  sbg:revision: 11
  sbg:revisionNotes: fixing do loop
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1620527988
  sbg:revision: 12
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1620551303
  sbg:revision: 13
  sbg:revisionNotes: time 2 cpu
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1620551955
  sbg:revision: 14
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1620555692
  sbg:revision: 15
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1620556058
  sbg:revision: 16
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1620559179
  sbg:revision: 17
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1620561152
  sbg:revision: 18
  sbg:revisionNotes: ''
- sbg:modifiedBy: einatgranot
  sbg:modifiedOn: 1620663061
  sbg:revision: 19
  sbg:revisionNotes: added annotation columns to variants
- sbg:modifiedBy: einatgranot
  sbg:modifiedOn: 1625667733
  sbg:revision: 20
  sbg:revisionNotes: Back to 'for' lopp due to minor discrepancies
sbg:sbgMaintained: false
sbg:validationErrors: []
