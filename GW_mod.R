##################################################################################
# The following script will implement GxW-BLUP, GOxW-BLUP, and GAOxW-BLUP models #
#
##################################################################################
setwd(" ") #directory as defined by the user

#### 1/ load necessary packages ####
library(caret)
library(BGLR)
library(tidyverse)
KinshipTransform <- function(Matrix) {
  nn <- ncol(Matrix)
  return((sum(diag(Matrix))-as.numeric(matrix(1,1,nn) %*% Matrix %*% matrix(1,nn,1))/nn)/(nn-1))
}


#### 2/ Here we will load the data ####
# list of DROPS environments
env_ids <- read.table("drops_ids.txt",header = T) #contained one column named "ids"


#Put your data in the right format :
grain_yield <- read.csv("~/work/megalmm/Scenario 1/grain.csv")
colnames(grain_yield) <- gsub(".grain_yield_15","",colnames(grain_yield))
grain_yield <- data.frame(grain_yield[,1],scale(grain_yield[,-1],center = T,scale = F))
grain_yield <- data.frame(grain_yield[,1],grain_yield[,match(env_covs$Experiment,colnames(grain_yield))])

colnames(grain_yield)[1] <- "GENOTYPE"
rownames(grain_yield) <- grain_yield$GENOTYPE

ngeno <- nrow(grain_yield)
ntrials <- ncol(grain_yield)[-1]

#here we will load DROPS Environmental Covariate Data
env_covs <- read.csv("~/work/GP_ECs/3b-Indices_Env_level.csv")

env_covs <- env_covs[match(colnames(grain_yield)[-1],env_covs$Experiment),-c(2,5,6,13)] # here removed unnecessary columns


# load multi-omics kinship matrices and scale them
K_G <- readRDS("genomic_kin.rds")
K_G <- K_G[match(grain_yield$GENOTYPE,rownames(K_G)),match(grain_yield$GENOTYPE,colnames(K_G))]
K_G <- K_G/KinshipTransform(K_G) #scaling the matrix

K_tr <- readRDS("transcriptomic_kin.rds")
K_tr <- K_tr[match(grain_yield$GENOTYPE,rownames(K_tr)),match(grain_yield$GENOTYPE,colnames(K_tr))]
K_tr <- K_tr/KinshipTransform(K_tr) #scaling the matrix

K_tw <- readRDS("twanscriptomic_kin.rds")
K_tw <- K_tw[match(grain_yield$GENOTYPE,rownames(K_tw)),match(grain_yield$GENOTYPE,colnames(K_tw))]
K_tw <- K_tw/Kinshiptwansform(K_tw) #scaling the matwix


K_pr <- readRDS("proteomic_kin.rds")
K_pr <- K_P[match(grain_yield$GENOTYPE,rownames(K_pr)),match(grain_yield$GENOTYPE,colnames(K_pr))]
K_pr <- K_pr/KinshipTransform(K_pr) #scaling the matrix


K_pw <- readRDS("pwoteomic_kin.rds")
K_pw <- K_P[match(grain_yield$GENOTYPE,rownames(K_pw)),match(grain_yield$GENOTYPE,colnames(K_pw))]
K_pw <- K_pw/KinshipTransform(K_pw) #scaling the matrix
###################################################################################################


# now we will change the data to long format to be used for next steps
data_format <- pivot_longer(data = grain_yield,cols = 2:ncol(grain_yield),names_to = "Env",values_to = "pheno")

# here we will add integer based IDs for environments and genotypes
data_format$IndE <- as.integer(factor(data_format$Env,labels = c(1:ntrials)))
data_format$IndG <- as.integer(factor(data_format$GENOTYPE,labels = c(1:ngeno)))

#data_format is a dataframe, with a column pheno corresponding to all the traits one after the other,
# a column of integer IndE indicating the environment


#-------------------------------------------------------------------------------
# 3/ Compute environmental similarity matrix W

# ECs are env covariates, each column corresponds to an environment (in the same order as in data_format)
EC <- env_covs

rownames(EC)  <- EC$Experiment

EC <- EC[,-1]

EC_cs = scale(EC,center=T,scale=T)



distce=dist(x=EC_cs, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)

W=1-as.matrix(distce)/(max(distce))



covEC4 = W[data_format$IndE,data_format$IndE]

#-------------------------------------------------------------------------------
# 4/ Compute the covariance matrix for the random GxW effect

# design matrix Z

Z=matrix(0,length(data_format$pheno),ngeno)   # ngeno is the number of genotypes
for (i in 1:nrow(Z)) {
  Z[i,data_format$IndG[i]]=1
}

cov_GW = Z %*% K_G %*% t(Z) * covEC4   # K is the kinship matrix in the correct order (same as in data_format)
cov_GW = cov_GW/KinshipTransform(cov_GW)
EVD_cov_GW=eigen(cov_GW)   # eigen value decomposition
rm(cov_GW)

#Twd
cov_trW = Z %*% K_tr %*% t(Z) * covEC4   # K is the kinship matrix in the correct order (same as in data_format)
EVD_cov_trW=eigen(cov_trW)   # eigen value decomposition
cov_trW = cov_trW/KinshipTransform(cov_trW)
rm(cov_TrW)

#Tww
cov_twW = Z %*% K_tw %*% t(Z) * covEC4   # K is the kinship matrix in the correct order (same as in data_format)
cov_twW = cov_twW/KinshipTransform(cov_twW)
EVD_cov_twW=eigen(cov_twW)   # eigen value decomposition
rm(cov_twW)

#Pwd
cov_prW = Z %*% K_pr %*% t(Z) * covEC4   # K is the kinship matrix in the correct order (same as in data_format)
cov_prW = cov_prW/KinshipTransform(cov_prW)
EVD_cov_prW=eigen(cov_prW)   # eigen value decomposition
rm(cov_prW)

#Pww
cov_pwW = Z %*% K_pw %*% t(Z) * covEC4   # K is the kinship matrix in the correct order (same as in data_format)
cov_pwW = cov_pwW/KinshipTransform(cov_pwW)
EVD_cov_pwW=eigen(cov_pwW)   # eigen value decomposition
rm(cov_pwW)



#-------------------------------------------------------------------------------
# 5/ Compute the covariance matrix for the random main G effect

kin2=kronecker(X=matrix(1,ntrials,ntrials), Y=K_G, FUN = "*", make.dimnames = FALSE)  # ntrials is the total nb of envts
EVD_kin2=eigen(kin2)   # eigen value decomposition
rm(kin2)


kin2=kronecker(X=matrix(1,ntrials,ntrials), Y=K_tr, FUN = "*", make.dimnames = FALSE)  # ntrials is the total nb of envts
EVD_kin2_tr=eigen(kin2)   # eigen value decomposition
rm(kin2)

kin2=kronecker(X=matrix(1,ntrials,ntrials), Y=K_tw, FUN = "*", make.dimnames = FALSE)  # ntrials is the total nb of envts
EVD_kin2_tw=eigen(kin2)   # eigen value decomposition
rm(kin2)


kin2=kronecker(X=matrix(1,ntrials,ntrials), Y=K_pr, FUN = "*", make.dimnames = FALSE)  # ntrials is the total nb of envts
EVD_kin2_pr=eigen(kin2)   # eigen value decomposition
rm(kin2)

kin2=kronecker(X=matrix(1,ntrials,ntrials), Y=K_pw, FUN = "*", make.dimnames = FALSE)  # ntrials is the total nb of envts
EVD_kin2_pw=eigen(kin2)   # eigen value decomposition
rm(kin2)

#-------------------------------------------------------------------------------
# 6/ Declare the model to the BGLR format

####for GxW Model
#ETA_EG_GXW=list(list(~factor(data_format$IndE), model='BRR'),list(V=EVD_kin2$vectors,d=EVD_kin2$values,model='RKHS'),
#               list(V=EVD_cov_GW$vectors,d=EVD_cov_GW$values,model='RKHS'))

#### for GOxW and GAOxW Models
ETA_EG_GOXW=list(list(~factor(data_format$IndE), model='BRR'),list(V=EVD_kin2$vectors,d=EVD_kin2$values,model='RKHS'),
                 list(V=EVD_kin2_tr$vectors,d=EVD_kin2_tr$values,model='RKHS'),list(V=EVD_kin2_tw$vectors,d=EVD_kin2_tw$values,model='RKHS'),
                 list(V=EVD_kin2_pr$vectors,d=EVD_kin2_pr$values,model='RKHS'),list(V=EVD_kin2_pw$vectors,d=EVD_kin2_pw$values,model='RKHS'),
                 list(V=EVD_cov_GW$vectors,d=EVD_cov_GW$values,model='RKHS'),list(V=EVD_cov_trW$vectors,d=EVD_cov_TrW$values,model='RKHS'),
                 list(V=EVD_cov_TwW$vectors,d=EVD_cov_TwW$values,model='RKHS'),list(V=EVD_cov_PrW$vectors,d=EVD_cov_PrW$values,model='RKHS'),
                 list(V=EVD_cov_PwW$vectors,d=EVD_cov_PwW$values,model='RKHS'))
rm(data_format,EVD_kin2,kin2,EVD_kin2_tr,EVD_kin2_tw,EVD_kin2_pr,EVD_kin2_pw,EVD_cov_trW,EVD_cov_twW,EVD_cov_prW,EVD_cov_pwW)
#we no longer need the above objects as we have everything that we need in ETA_EG_GXW or ETA_EG_GOXW

#-------------------------------------------------------------------------------
# 7/ Apply the model in reps and folds

###make folds for CV###
nreps <- 5
nfolds <- 5
set.seed(1)
folds <- replicate(nreps,createFolds(rnorm(grain_yield$GENOTYPE),k=nfolds,list = F,returnTrain = F))

######
cors_rf <- list() #to save the final results of all reps and folds

for (reps in 1:nreps) { # to scroll through different reps
  
  
  cors_folds <- list() #to save the results of different folds
  
  for (i in 1:nfolds) {
    
    
    pheno <- grain_yield # this one will be used to mask the test data
    pheno[which(folds==i),2:ncol(pheno)]=NA # this introduces NAs for the test set
    
    pheno_backup <- grain_yield # this one have the phenotypes for test set, and will be used to get predictive abilities
    
    # now we will change the data to long format
    data_format <- pivot_longer(data = pheno,cols = 2:ncol(pheno),names_to = "Env",values_to = "pheno")
    data_format_backup <- pivot_longer(data = pheno_backup,cols = 2:ncol(pheno_backup),names_to = "Env",values_to = "pheno")
    
    # here we will add integer based IDs for environments and genotypes
    data_format$IndE <- as.integer(factor(data_format$Env,labels = c(1:ntrials)))
    data_format$IndG <- as.integer(factor(data_format$GENOTYPE,labels = c(1:ngeno)))
    
    #data_format is a dataframe, with a column pheno corresponding to all the traits one after the other,
    # a column of integer IndE indicating the environment
    
    pheno_backup <- pheno_backup[,-1] #now we no longer need the geno id column
    
    
    #-------------------------------------------------------------------------------
    # 5/ Run BGLR
    
    N_iter=20000
    N_burnin=5000
    N_thin=2
    
    a=Sys.time()
    print(a)
    res=BGLR(y=data_format$pheno, response_type = "gaussian",ETA = ETA_EG_GXW, nIter = N_iter,
             burnIn = N_burnin, thin = N_thin, S0 = NULL,
             df0 =5, R2 = 0.5, weights = NULL,
             verbose = FALSE, rmExistingFiles = TRUE, groups=NULL)
    b=Sys.time()
    
    
    res_out <- data_format_backup
    
    res_out$pred <- res$yHat # These are the predictions
    
    
    pred <- data.frame(pivot_wider(data = res_out[,c(1,2,4)],names_from =  "Env",values_from = "pred"))
    
    rownames(pred) <- pred$GENOTYPE
    
    pred <- pred[,-1] #removing the geno id column
    
    
    cors_folds[[i]] = diag(cor(pheno_backup[which(folds==i),]
                               ,pred[which(folds==i),],use='p'))
    
    print(paste0("status:"," ",i))
    
    
  }
  
  cors_rf[[reps]] <- cors_folds
}
saveRDS(cors_rf,file ="res_all.RDS")
