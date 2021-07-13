cwlVersion: v1.2
class: CommandLineTool
label: Ancestry MAF estimate local proportion
$namespaces:
  sbg: https://sevenbridges.com

requirements:
- class: ShellCommandRequirement
- class: ResourceRequirement
  coresMin: $(inputs.cpu * inputs.cpu_multiplyer)
  ramMin: 8000
- class: DockerRequirement
  dockerPull: |-
    images.sb.biodatacatalyst.nhlbi.nih.gov/einatgranot/ancestral_maf_admixed_population:20210507
- class: InitialWorkDirRequirement
  listing:
  - entryname: freq_ancestry_local.R
    writable: false
    entry: |2+

      #Calculate SOL ancetries MAFs based on Tamar's function using global proportion ancestries (without prior knowledge on Ancestral MAFs)

      library(GWASTools)
      library(gdsfmt)
      library(SeqArray)
      library(SeqVarTools)
      library(dplyr)
      library(stringr)
      library(rtracklayer) #need this

      source("cwl_input.R")

      plist<-read.csv(plist_path, header=T)

      #get command line args
      args = commandArgs(trailingOnly=TRUE)
      myindex <- args[1]

      #j is chromosomes
      j  <- plist[myindex,1]
      j
      #b is blocks
      b <- plist[myindex,2]
      b

      ########################Define functions for allele frequency caculation in admixed populations#################################
      source('decompose_two_alleles_one_person.R', chdir = TRUE)
      source('generate_smoothing_observations.R', chdir = TRUE)
      source('estimate_frequencies.R', chdir = TRUE)
      source('prep_dat_for_binomial_likelihood.R', chdir = TRUE)
      source('dbinom_approx.R', chdir = TRUE)
      source('estimate_frequencies_search_boundary.R', chdir = TRUE)
      source('estimate_frequencies_w_known_freqs.R', chdir = TRUE)

      #Read SOL global proportion ancestry
      write("Read SOL local proportion ancestry", stderr())

      LAI_amer_annot <- read.csv(local_file_path_amer)
      LAI_afr_annot <- read.csv(local_file_path_afr)
      LAI_eur_annot <- read.csv(local_file_path_eur)

      write("Read list of unrelated individuals", stderr())
      #Read list of unrelated individuals
      sol_unrel<-read.table(sol_unrel_path)[1] #10223
      names(sol_unrel)[1] <-"SUBJECT_ID"

      #Create a list of the final SOL participants available for MAF calculations 
      LAI_amer_annot<-merge(sol_unrel, LAI_amer_annot, by="SUBJECT_ID", all=F)#10157
      LAI_afr_annot<-merge(sol_unrel, LAI_afr_annot, by="SUBJECT_ID", all=F)#10157
      LAI_eur_annot<-merge(sol_unrel, LAI_eur_annot, by="SUBJECT_ID", all=F)#10157

      #Read gds imputed genotypes file for chromosome j
      gds <- seqOpen(chromosome_full_path, FALSE)

      #Subset and sort local ancestry proportions to the sample.id in the gds file
      sample.id <- seqGetData(gds, "sample.id")#11928
      sample.id_logical<-sample.id %in% LAI_eur_annot$SUBJECT_ID#9512
      sample.id<-as.data.frame(sample.id[sample.id_logical])#9512
      names(sample.id)[1]<-"SUBJECT_ID"

      LAI_amer_annot<-inner_join(sample.id, LAI_amer_annot, by="SUBJECT_ID")#9512
      LAI_afr_annot<-inner_join(sample.id, LAI_afr_annot, by="SUBJECT_ID")#9512
      LAI_eur_annot<-inner_join(sample.id, LAI_eur_annot, by="SUBJECT_ID")#9512

      snp.rs.id <- seqGetData(gds, "variant.id")
      nvar <- length(snp.rs.id)
      chr.id <-seqGetData(gds, "chromosome")
      position.id <-seqGetData(gds, "position")
      alleles<-seqGetData(gds, "allele")
      allele_a<-str_split_fixed(alleles, ",",2)[,1]
      allele_b <-str_split_fixed(alleles, ",",2)[,2]

      #Partition the gds genotyped data to segments with 3000 variants

      # compute the number of blocks
      b.size <- 3000

      # set filters by blocks and read data
      b.start.inds <- seq(1, nvar, by = b.size)
      b.start.inds
      b.end.inds <- seq(b.size, nvar, by = b.size)
      b.end.inds
      if (length(b.start.inds) > length(b.end.inds)) b.end.inds <- c(b.end.inds, nvar)
      b.size_seg<-c(rep(b.size,length(b.start.inds)-1),b.end.inds[length(b.start.inds)]-b.start.inds[length(b.start.inds)]+1)

      # set sample and variant filters
      seqSetFilter(gds, c(b.start.inds[b]:b.end.inds[b]), sample.sel=sample.id_logical, verbose=TRUE)
      gen<-as.data.frame(imputedDosage(gds, dosage.field="DS", use.names=TRUE))
      closefn.gds(gds)

      #Check compatibility of order of individuals in the gds file and in prop_mat
      #table(rownames(gen)==LAI_eur_annot$SUBJECT_ID)#9512

      snp.chromosome_gen <-chr.id[b.start.inds[b]:b.end.inds[b]]
      snp.position_gen <-position.id[b.start.inds[b]:b.end.inds[b]]


      #Liftover the hg38 positions to the hg19- so it could be compared to the LAI
      #Read liftover file 

      chain <- import.chain(chain_file_path)

      gr <- GRanges(seqnames=paste0("chr",snp.chromosome_gen), ranges=IRanges(start=snp.position_gen, end=snp.position_gen))
      lift <- liftOver(gr, chain)

      # report out number of successful liftovers
      # one.rslt <- sum(elementNROWS(lift) == 1)
      # no.rslt <- sum(elementNROWS(lift) == 0)
      # message(one.rslt, " input positions successfully converted; ", no.rslt, " positions failed conversion (will be set to NA)")

      # extract converted chrom positions in hg19
      lift.pos <- as.integer(start((lift)))
      lift.pos[is.na(lift.pos)] <- 0
      snp.position_gen_19 <- lift.pos

      #Read LAI annotation file
      admixAnnot <- getobj(lai_annot_file_path)
      head(pData(admixAnnot))
      dim(admixAnnot)
      LAI_pos<-pData(admixAnnot)


      ancestry_mafs <-matrix(NA,nrow=ncol(gen),ncol=12)

      for (i in c(1:ncol(gen))){
        
        allele_counts<-as.numeric(gen[,i])
        var_pos<-snp.position_gen_19[i]
        
        LAI_position<-as.numeric(LAI_pos[j==LAI_pos$chromosome & var_pos>LAI_pos$pos.start & var_pos<LAI_pos$pos.end,]$snpID)
        if (length(LAI_position)==1){
          prop_mat<-as.matrix(cbind(LAI_afr_annot[,LAI_position+1],LAI_eur_annot[,LAI_position+1],LAI_amer_annot[,LAI_position+1]))
          colnames(prop_mat) <-c("afr","eur","amer")
          #Change the LAI counts to proportions
          prop_mat<-as.matrix(prop_mat/2)
          
          
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
      }

      ancestry_mafs<- data.frame(ancestry_mafs)

      names(ancestry_mafs)<-c("Africa_est", "Africa_low_CI", "Africa_high_CI","Europe_est", "Europe_low_CI", "Europe_high_CI","America_est", "America_low_CI", "America_high_CI", "nll","boundary", "flag")

      ancestry_mafs$POS_hg38<-snp.position_gen
      ancestry_mafs$CHR<-snp.chromosome_gen
      ancestry_mafs$POS_hg19<-snp.position_gen_19
      ancestry_mafs$allele_a<-allele_a[b.start.inds[b]:b.end.inds[b]]
      ancestry_mafs$allele_b<-allele_b[b.start.inds[b]:b.end.inds[b]]
      write.csv(ancestry_mafs, paste0("chr",j,"_",b,"_gds_ancestry_mafs_local_121120.csv"), quote=F, row.names=F)

  - entryname: cwl_input.R
    writable: false
    entry: |-
      plist_path="$(inputs.plist.path)" 
      local_file_path_amer="$(inputs.local_ancestry_amerindian.path)"
      local_file_path_afr="$(inputs.local_ancestry_african.path)"
      local_file_path_eur="$(inputs.local_ancestry_european.path)"
      chromosome_full_path="$(inputs.chromosome_gds.path)"
      chromosome_number=$(inputs.chromosome_gds.path.split("chr")[1].split("_")[0])
      sol_unrel_path="$(inputs.sol_unrel.path)"
      chain_file_path="$(inputs.chain_file.path)"
      lai_annot_file_path="$(inputs.lai_annotation.path)"
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
    entry: |+
      # Generate artificial observations with alleles from each provided ancestries.
      # These observations are later added to the real data to protect from reaching boundary condition
      # and failure of the estimation function. 

      generate_smoothing_observations <- function(ancestry_names){
        
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
    entry: |2+

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
        
        prep_dat <- prep_dat_for_binomial_likelihood(allele_counts, prop_mat,
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
          smoothing_data <- generate_smoothing_observations(colnames(prop_mat))
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
      estimate_frequencies_after_prep <- function(allele_counts,
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
      prep_dat_for_binomial_likelihood <- function(allele_counts, prop_mat,  
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
        
        prep_dat <- prep_dat_for_binomial_likelihood(allele_counts, prop_mat,
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
          smoothing_data <- generate_smoothing_observations(colnames(prop_mat))
          # add only the rows corresponding to the ancestries with unknown frequencies
          inds_keep <- (n_known*2 + 1):((n_known + n_unknown)*2)
          prop_mat <- rbind(prop_mat, smoothing_data$simulated_prop_mat[inds_keep,])
          allele_counts <- c(allele_counts, smoothing_data$simulated_allele_count[inds_keep])
          max_counts <- c(max_counts, rep(1, length(inds_keep)))
        }
        
        
        return(estimate_frequencies_w_known_freq_after_prep(allele_counts, 
                                                             prop_mat, 
                                                             max_counts,
                                                             known_freqs, 
                                                             n_unknown,
                                                             n_known,
                                                             low_freq_bound, 
                                                             high_freq_bound,
                                                             confidence))
        
      }



      estimate_frequencies_w_known_freq_after_prep <- function(allele_counts, 
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
    entry: |+
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
        prep_dat <- prep_dat_for_binomial_likelihood(allele_counts, prop_mat,
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
          smoothing_data <- generate_smoothing_observations(colnames(prop_mat))
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
          
          cur_res <- tryCatch({estimate_frequencies_after_prep(allele_counts = allele_counts, 
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
          tryCatch({estimate_frequencies_w_known_freq_after_prep(allele_counts = allele_counts, 
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
- id: local_ancestry_amerindian
  label: Amerindian local ancestry file
  doc: Amerindian local ancestry file SOL
  type: File?
- id: chromosome_gds
  doc: 1-22 chromosomes
  type: File?
- id: lai_annotation
  label: LAI annotation file
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
- id: local_ancestry_african
  label: African local ancestry file
  doc: African local ancestry file SOL
  type: File?
- id: local_ancestry_european
  label: European local ancestry file
  doc: European local ancestry file SOL
  type: File?
- id: chain_file
  label: chain file
  doc: chain file liftover from hg38 to hg19
  type: File?
- id: cpu_multiplyer
  type: int

outputs:
- id: ancestry_csv
  type: File?
  outputBinding:
    glob: '*.csv'
    outputEval: $(inheritMetadata(self, inputs.chromosome_gds))
- id: workspace
  type: File?
  outputBinding:
    glob: '*.RData'
    outputEval: $(inheritMetadata(self, inputs.chromosome_gds))
stdout: standard.out

baseCommand:
- Rscript
- freq_ancestry_local.R

hints:
- class: sbg:SaveLogs
  value: '*.R'
- class: sbg:SaveLogs
  value: standard.out
- class: sbg:SaveLogs
  value: '*.RData'
id: |-
  dave/build-ancestry-maf-admixed-population/ancestry-maf-estimate-local-proportion/37
sbg:appVersion:
- v1.2
sbg:content_hash: a5d415dad685af7e0b95b5dce75224c1343ba42aad184c9d5b1e6b4ac3d412be8
sbg:contributors:
- dave
- einatgranot
sbg:createdBy: dave
sbg:createdOn: 1620395860
sbg:id: |-
  dave/build-ancestry-maf-admixed-population/ancestry-maf-estimate-local-proportion/37
sbg:image_url:
sbg:latestRevision: 37
sbg:modifiedBy: einatgranot
sbg:modifiedOn: 1625666832
sbg:project: dave/build-ancestry-maf-admixed-population
sbg:projectName: 'BUILD: Ancestry MAF Admixed Population'
sbg:publisher: sbg
sbg:revision: 37
sbg:revisionNotes: Back to 'for' loop- due to minor discrepancies
sbg:revisionsInfo:
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1620395860
  sbg:revision: 0
  sbg:revisionNotes: |-
    Copy of einatgranot/2020-freq-estimation-latino-ancestries/ancestry-maf-estimate-local-proportion/18
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1620678038
  sbg:revision: 1
  sbg:revisionNotes: added some part of the do loop
- sbg:modifiedBy: einatgranot
  sbg:modifiedOn: 1620745764
  sbg:revision: 2
  sbg:revisionNotes: change for loop with something faster
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1620825280
  sbg:revision: 3
  sbg:revisionNotes: just saving workspace
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1620825769
  sbg:revision: 4
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1620826303
  sbg:revision: 5
  sbg:revisionNotes: |-
    images.sb.biodatacatalyst.nhlbi.nih.gov/einatgranot/ancestral_maf_admixed_population:20210507
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1620826774
  sbg:revision: 6
  sbg:revisionNotes: rm tic toc
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1620827266
  sbg:revision: 7
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1620852542
  sbg:revision: 8
  sbg:revisionNotes: back to older script - just truncated before for loop and saving
    workspace
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1621045569
  sbg:revision: 9
  sbg:revisionNotes: change function names
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1621108905
  sbg:revision: 10
  sbg:revisionNotes: fixing function names
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1621110769
  sbg:revision: 11
  sbg:revisionNotes: fixing function names
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1621127292
  sbg:revision: 12
  sbg:revisionNotes: added in full do loop
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1621127357
  sbg:revision: 13
  sbg:revisionNotes: save csv
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1621133014
  sbg:revision: 14
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1621133652
  sbg:revision: 15
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1621134584
  sbg:revision: 16
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1621135212
  sbg:revision: 17
  sbg:revisionNotes: dplyr::rename
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1621135920
  sbg:revision: 18
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1621137488
  sbg:revision: 19
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1621175982
  sbg:revision: 20
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1621176902
  sbg:revision: 21
  sbg:revisionNotes: unbound cpu multiplier
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1621190139
  sbg:revision: 22
  sbg:revisionNotes: ''
- sbg:modifiedBy: einatgranot
  sbg:modifiedOn: 1621196960
  sbg:revision: 23
  sbg:revisionNotes: removed levels of frequency boundary
- sbg:modifiedBy: einatgranot
  sbg:modifiedOn: 1621197124
  sbg:revision: 24
  sbg:revisionNotes: corrected order of boundaries
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1621214694
  sbg:revision: 25
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1621215387
  sbg:revision: 26
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1621216936
  sbg:revision: 27
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1621217342
  sbg:revision: 28
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1621219080
  sbg:revision: 29
  sbg:revisionNotes: removing nas
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1621220801
  sbg:revision: 30
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1621222378
  sbg:revision: 31
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1621223664
  sbg:revision: 32
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1621242312
  sbg:revision: 33
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1621455498
  sbg:revision: 34
  sbg:revisionNotes: reformat na filter
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1621475195
  sbg:revision: 35
  sbg:revisionNotes: ''
- sbg:modifiedBy: einatgranot
  sbg:modifiedOn: 1621610436
  sbg:revision: 36
  sbg:revisionNotes: added pos19 pos38 columns
- sbg:modifiedBy: einatgranot
  sbg:modifiedOn: 1625666832
  sbg:revision: 37
  sbg:revisionNotes: Back to 'for' loop- due to minor discrepancies
sbg:sbgMaintained: false
sbg:validationErrors: []
