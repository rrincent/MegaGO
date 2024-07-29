#######################################################################################
################# The Script to Apply MegaLMM, MegaGO, or MegaGAO #####################
#######################################################################################
# adapted from https://github.com/deruncie/MegaLMM/blob/master/vignettes/MultiEnvironmentTrial.Rmd
#### loading packages ####
library(MegaLMM)
library(rrBLUP)
library(caret)
#--------------------------------------------------------------------------------------

#### loading input data ####

# phenotypes

phenotypes <- read.csv("DROPS.csv") # grain yield or platform traits. The first column has genotype IDs. Rest of them are traits
genotypes <- data.frame(phenotypes[,1]);colnames(genotypes) <- "GENOTYPE" # the first column contains IDs

# genomic kinship matrix

kin_G <- readRDS("K.RDS") # genomic kinship matrix
kin_G <- kin_G[genotypes$GENOTYPE,genotypes$GENOTYPE]

# omics data

tww <- readRDS("tww04.RDS");tww <- tww[genotypes$GENOTYPE,]  # WW transcript expressions (additive or standard)
twd <- readRDS("twd04.RDS");twd <- twd[genotypes$GENOTYPE,]  # WD transcript expressions (additive or standard)
pww <- readRDS("pww04.RDS");pww <- pww[genotypes$GENOTYPE,]  # WW protein abundances (additive or standard)
pwd <- readRDS("pwd04.RDS");pwd <- pwd[genotypes$GENOTYPE,]  # WD protein abundances (additive or standard)

#--------------------------------------------------------------------------------------------------
#### creating Y matrix for MegaLMM, MegaGO/MegaGAO ####

pheno_mat_Y_megalmm <- phenotypes[,-1] #id column is removed

pheno_mat_Y_megaGO <- cbind(phenotypes[,-1],tww,twd,pww,pwd) #data for megaGO/megaGAO

NphenoCols <- ncol(phenotypes)-1

#--------------------------------------------------------------------------------------------------

#### creating repeats and folds ####

nreps <- 1
nfolds <- 5
set.seed(1)
folds <- replicate(nreps,createFolds(rnorm(genotypes$GENOTYPE),k=nfolds,list = F,returnTrain = F))

#--------------------------------------------------------------------------------------------------

#### lists to store the results ####

gblup_cors_pearson <- list()
Eta_mean_cors_pearson <- list()
U_hat_cors_pearson <- list()

#--------------------------------------------------------------------------------------------------


for (i in 1:nfolds) {
  
  pheno <- pheno_mat_Y_megalmm #change it to pheno_mat_Y_megaGO to run megaGO or megaGAO
  pheno[which(folds==i),1:NphenoCols]=NA # it masks the data for test set.
  
  #GBLUP predictions
  GBLUP_predictions = matrix(NA,nrow(pheno),NphenoCols,dimnames = dimnames(pheno[,1:NphenoCols]))
  for(j in 1:NphenoCols) {
    GBLUP_predictions[,j] = mixed.solve(y = pheno[,j],
                                        X = NULL,
                                        K = kin_G)$u
  }
  
  gblup_cors_pearson[[i]] <- diag(cor(pheno_mat_Y_megalmm[which(folds==i),1:NphenoCols],
                                      GBLUP_predictions[which(folds==i),1:NphenoCols],use = "p"))
  
  
  
  #### running megalmm or megaGO/megaGAO ####
  run_parameters = MegaLMM_control(
    h2_divisions = 20, # Each variance component is allowed to explain between 0% and 100% of the
    # total variation. How many segments should the range [0,100) be divided 
    # into for each random effect?
    burn = 0,  
    thin = 2,
    K = 10,#150 # number of latent factors. increase when running megaGO or MegaGAO (we used 150)
    scale_Y = T #megalmm scales and unscales automatically when making predictions. 
    #if we don't scale here, it will mess with the factor loadings due to different scales, e.g., some  traits will overload on a single factor
  )
  
  # setup model for megalmm
  
  MegaLMM_state = setup_model_MegaLMM(
    Y = pheno,  
    formula = ~ 1 + (1|GENOTYPE),  #y=mu+g+e 
    # This is syntax like lme4 for mixed effect models. 
    data = genotypes, #this file should have the same order as pheno file. We can add different fixed effects here and then specify them in the model.        
    # the data.frame with information for constructing the model matrices
    relmat = list(GENOTYPE = kin_G), #can have multiple kins
    run_parameters=run_parameters,
    # This list of control parameters created above
    run_ID = sprintf('MegaLMM_fold_%02d',i)
    # A run identifier. The function will create a folder with this name 
    # and store lots of useful data inside it
  )
  
  
  # setting priors sampler
  
  Lambda_prior = list(      # prior for factor loadings. 
    sampler = sample_Lambda_prec_ARD,
    Lambda_df = 3,
    delta_1   = list(shape = 2,  rate = 1),
    delta_2   = list(shape = 3, rate = 1)
  )
  
  
  # setup megalmm priors
  
  priors = MegaLMM_priors(
    tot_Y_var = list(V = 0.5,   nu = 5),      
    # Prior variance of trait residuals after accounting for fixed effects and factors
    # See MCMCglmm for meaning of V and nu
    tot_F_var = list(V = 1, nu = 1e6), 
    # Prior variance of factor traits. This is included to improve MCMC mixing, 
    # but can be turned off by setting nu very large
    h2_priors_resids_fun = function(h2s,n)  1,  
    # Function that returns the prior density for any value of the h2s vector 
    # (ie the vector of random effect proportional variances across all random effects. 
    # 1 means constant prior. 
    # n is the number of h2 divisions above (here=20)
    # 1-n*sum(h2s)/n linearly interpolates between 1 and 0, 
    # giving more weight to lower values
    h2_priors_factors_fun = function(h2s,n) 1, 
    # See above. 
    # sum(h2s) linearly interpolates between 0 and 1,
    # giving more weight to higher values
    # Another choice is one that gives 50% weight to h2==0: ifelse(h2s == 0,n,n/(n-1))
    Lambda_prior = Lambda_prior
    # from above
  )
  
  ###assign priors to the state object
  MegaLMM_state = set_priors_MegaLMM(MegaLMM_state,priors)
  
  ##map: forms different mappings of the data to identify missing values
  maps = make_Missing_data_map(MegaLMM_state,max_NA_groups = ncol(pheno)+1,verbose=F) #ncol(pheno)+1 if pheno doesnt have id column
  maps$map_results
  MegaLMM_state = set_Missing_data_map(MegaLMM_state,maps$Missing_data_map_list[[1]])
  
  
  ###random starting values
  MegaLMM_state = initialize_variables_MegaLMM(MegaLMM_state)
  
  ####memory estimation
  estimate_memory_initialization_MegaLMM(MegaLMM_state)
  
  ###calculate useful matrices (precalculations)
  MegaLMM_state = initialize_MegaLMM(MegaLMM_state,verbose = T)
  
  ####
  MegaLMM_state$Posterior$posteriorMean_params
  
  #####which values to keep
  MegaLMM_state$Posterior$posteriorSample_params = c('F_h2','resid_h2','tot_Eta_prec','F','Lambda') #lambda can take a lot space when number of traits increase 
  MegaLMM_state$Posterior$posteriorMean_params = 'Eta_mean'
  
  ####calculating predicted values
  tt=1:NphenoCols #target traits
  ntt = (1+length(tt)):ncol(pheno) #not target traits. we need this line if we are running megago or megagao
  MegaLMM_state$Posterior$posteriorFunctions = list(
    U = 'U_F %*% Lambda[,tt] + U_R[,tt]'#,
    ##  G_sub = 't(Lambda[,tt]) %*% diag(F_h2[1,]) %*% Lambda[,ntt]', #this is when we want to get correlations between traits and omics. remove all ntt and even tt when running megaLMM
    ##  G_diag = 'colSums(F_h2[1,]*Lambda^2)+resid_h2[1,]/tot_Eta_prec[1,]',#this is when we want to get correlations between traits and omics. remove all ntt and even tt when running megaLMM
    
    ##  P_diag = 'colSums(Lambda^2)+1/tot_Eta_prec[1,]'
    #U = 'U_F %*% Lambda + U_R'
    #G = 't(Lambda) %*% diag(F_h2[1,]) %*% Lambda + diag(resid_h2[1,]/tot_Eta_prec[1,])',
    #R = 't(Lambda) %*% diag(1-F_h2[1,]) %*% Lambda + diag((1-resid_h2[1,])/tot_Eta_prec[1,])',
    #h2 = '(colSums(F_h2[1,]*Lambda^2)+resid_h2[1,]/tot_Eta_prec[1,])/(colSums(Lambda^2)+1/tot_Eta_prec[1,])'
  )
  ####Now that we've decided which values to save, we initialize the posterior database:
  MegaLMM_state = clear_Posterior(MegaLMM_state) 
  #we can check its memory estimate
  estimate_memory_posterior(MegaLMM_state,100)
  
  ################
  ##model fitting or GIBBS sampling
  n_iter = 500
  for(k in 1:5) {
    print(sprintf('Burnin run %d',k))
    # Factor order doesn't "mix" well in the MCMC.
    # We can help it by manually re-ordering from biggest to smallest
    MegaLMM_state = reorder_factors(MegaLMM_state,drop_cor_threshold = 0.6)
    # clear any previous collected samples because we've re-started the chain 
    MegaLMM_state = clear_Posterior(MegaLMM_state)
    # Draw n_iter new samples, storing the chain
    MegaLMM_state = sample_MegaLMM(MegaLMM_state,n_iter)
    # make diagnostic plots
    #traceplot_array(MegaLMM_state$Posterior$Lambda,name = paste0("Lambda_",i,".pdf"))
    ##fold_ID_matrix <- matrix(pheno[,-1])
    #traceplot_array(MegaLMM_state$Posterior$U,name = paste0("U_",i,".pdf"),
    #                facet_dim = 3 )#,mask = fold_ID_matrix != 1)
    #print(sprintf('Completed %d burnin samples', MegaLMM_state$current_state$nrun))
  }
  MegaLMM_state = clear_Posterior(MegaLMM_state)
  
  #######################
  n_iter = 100
  for(k in 1:25) {
    print(sprintf('Sampling run %d',k))
    MegaLMM_state = sample_MegaLMM(MegaLMM_state,n_iter) 
    MegaLMM_state = save_posterior_chunk(MegaLMM_state)
    print(MegaLMM_state)
  }
  
  
  ###loading posteriors
  # Lambda_samples = load_posterior_param(MegaLMM_state,'Lambda')
  U_samples = load_posterior_param(MegaLMM_state,'U')
  ##  G_sub = load_posterior_param(MegaLMM_state,'G_sub')
  ##  G_diag = load_posterior_param(MegaLMM_state,'G_diag')
  
  ##  MegaLMM_state$Posterior$G_sub = G_sub
  ##  MegaLMM_state$Posterior$G_diag = G_diag
  
  ##  G_cor_sub = get_posterior_FUN(MegaLMM_state, sweep(sweep(G_sub,1,sqrt(G_diag[tt,1]),'/'),2,sqrt(G_diag[ntt,1]),'/'))
  ##  MegaLMM_state$Posterior$posteriorSample_params = c(MegaLMM_state$Posterior$posteriorSample_params,'G_sub','G_diag')
  ##  G_cor_sub = get_posterior_mean(MegaLMM_state, {
  #recover()
  ##    G_sub_std_tt = sweep(G_sub,1,sqrt(G_diag[tt,1]),'/')
  ##    G_sub_std_tt_ntt = sweep(G_sub_std_tt,2,sqrt(G_diag[ntt,1]),'/')
  ##    G_sub_std_tt_ntt
  ##  })
  
  #Image(G_cor_sub) # genetic correlations
  
  #plot(G_cor_sub[1,],G_cor_sub[80,])
  #cor(G_cor_sub[1,],t(G_cor_sub[70:92,]))
  
  ##  contributions[[i]] <- G_cor_sub
  #Image(G_cor_sub)
  
  
  dim(U_samples)
  U_hat = get_posterior_mean(U_samples)
  Eta_mean = load_posterior_param(MegaLMM_state,'Eta_mean')
  
  
  ####accuracies
  U_hat_sel <- U_hat#[,colnames(U_hat) %in% sel_envs$Experiment]
  rownames(U_hat_sel) <- gsub("::GENOTYPE","",rownames(U_hat_sel))
  
  U_hat_cors_pearson[[i]] = diag(cor(pheno_mat_Y_megalmm[which(folds==i),1:NphenoCols]
                                     ,U_hat[which(folds==i),1:NphenoCols],use='p'))
  rownames(Eta_mean) <- rownames(U_hat)
  Eta_mean_sel <- Eta_mean#[,colnames(Eta_mean) %in% sel_envs$Experiment]
  rownames(Eta_mean_sel) <- gsub("::GENOTYPE","",rownames(Eta_mean_sel))
  
  
  Eta_mean_cors_pearson[[i]] = diag(cor(pheno_mat_Y_megalmm[which(folds==i),1:NphenoCols]
                                        ,Eta_mean[which(folds==i),1:NphenoCols],use='p'))
  
  
  
  
}

save(Eta_mean_cors_pearson,gblup_cors_pearson,
     U_hat_cors_pearson,
     file="resultMegaLMM.Rdata")
