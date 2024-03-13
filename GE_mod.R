##################################################################################################
#           The following script implements GxE-BLUP GOxE-BLUP and GAO-E-BLUP models   		     #
#   By default it runs CV1, but we can modify line 57 according to the given comment to run CV2  #
##################################################################################################

setwd("") #Working directory defined by the user.
#### 1/ loading necessary packages####
library(sommer)
library(caret)
library(tidyverse)
KinshipTransform <- function(Matrix) { # function for scaling kinship
  nn <- ncol(Matrix)
  return((sum(diag(Matrix))-as.numeric(matrix(1,1,nn) %*% Matrix %*% matrix(1,nn,1))/nn)/(nn-1))
}

###################################

#### 2/ Loading Data ####
# Environment names
env_ids <- read.table("env_ids.txt",header=T)

# GY Data. The GY data is in wide format, with the first column for the genotype ids, and the remaining for Environments/Trials. (NxT)
grain_yield <- read.csv("grain_yield.csv")
grain_yield <- data.frame(grain_yield[,1],scale(grain_yield[,-1])) #if scaling is required
grain_yield <- data.frame(grain_yield[,1],grain_yield[,match(env_ids$ids,colnames(grain_yield))])

colnames(grain_yield)[1] <- "GENOTYPE"
rownames(grain_yield) <- grain_yield$GENOTYPE

n.trials <- ncol(grain_yield)-1 #specifies the number of GY trials. Can be set manually to run a customized example
n.lines <- nrow(grain_yield) # specifies the number of lines/hybrids. Can be set manually to run customized examples

##### kinship matrices####
# genomic
geno <- readRDS("kin.RDS")
K_G <- geno$A #additive kinship matrix
K_G <- K_G[match(grain_yield$GENOTYPE,rownames(K_G)),match(grain_yield$GENOTYPE,colnames(K_G))]
K_G <- K_G/KinshipTransform(K_G)
# Transcriptomic
K_tr <- readRDS("kin_tr.RDS")
K_tr <- K_tr[match(grain_yield$GENOTYPE,rownames(K_tr)),match(grain_yield$GENOTYPE,colnames(K_tr))]
K_tr <- K_tr/KinshipTransform(K_tr)

K_tw <- readRDS("kin_tw.RDS")
K_tw <- K_tw[match(grain_yield$GENOTYPE,rownames(K_tw)),match(grain_yield$GENOTYPE,colnames(K_tw))]
K_tw <- K_tw/KinshipTransform(K_tw)
# Proteomic
K_pr <- readRDS("kin_pr.RDS")
K_pr <- K_pr[match(grain_yield$GENOTYPE,rownames(K_pr)),match(grain_yield$GENOTYPE,colnames(K_pr))]
K_pr <- K_pr/KinshipTransform(K_pr)

K_pw <- readRDS("kin_pw.RDS")
K_pw <- K_pw[match(grain_yield$GENOTYPE,rownames(K_pw)),match(grain_yield$GENOTYPE,colnames(K_pw))]
K_pw <- K_pw/KinshipTransform(K_pw)

#### 3/ Setting CV #####
nreps <- 5
nfolds <- 5
set.seed(1) # we use the same seed for other models to ensure same fold-repeat combinations
folds <- replicate(nreps,createFolds(rnorm(grain_yield$GENOTYPE),k=nfolds,list = F,returnTrain = F))

#### 4/ Creating a list to store accuracies ####
ansCS_cors_rf <- list()

for(reps in 1:nreps){ #for different reps
selected_fold <- folds[,reps]
ansCS_cors <- list()

for(i in 1:nfolds){  #for different folds
pheno <- grain_yield #pheno is created for replacing test phenotypes with NAs

pheno[which(selected_fold==i),c(2:n.trials)]=NA # here c(2:n.trials) vector can be changed to selected trials.

# we will change the data to long format as the software needs it.
gy <- pivot_longer(pheno[,1:(n.trials+1)],cols = 2:ncol(pheno[,1:(n.trials+1)]),
                   names_to = "Env",values_to = "Yield")

# genotype and environment ids need to be changed to factors.
gy$GENOTYPE_G <- as.factor(gy$GENOTYPE) #for K_G
gy$GENOTYPE_tr <- as.factor(gy$GENOTYPE) #for K_tr
gy$GENOTYPE_tw <- as.factor(gy$GENOTYPE) #for K_tw
gy$GENOTYPE_pr <- as.factor(gy$GENOTYPE) #for K_pr
gy$GENOTYPE_pw <- as.factor(gy$GENOTYPE) #for K_pw
gy$Env <- as.factor(gy$Env)

# here we will create the E matrix and form its Kronecker product with Kinship matrix
E <- diag(length(unique(gy$Env)))
rownames(E) <- colnames(E) <- unique(gy$Env)
EK_G <- kronecker(E,K_G, make.dimnames = TRUE)
EK_tr <- kronecker(E,K_tr, make.dimnames = TRUE)
EK_tw <- kronecker(E,K_tw, make.dimnames = TRUE)
EK_pr <- kronecker(E,K_pr, make.dimnames = TRUE)
EK_pw <- kronecker(E,K_pw, make.dimnames = TRUE)


#model fitting with mmer for GxE-BLUP
# NOTE: GxE-blup can be run by uncommenting the lines here

#ansCS <- mmer(Yield~Env,
#              random= ~ vsr(GENOTYPE_G, Gu=K_G)+
#                vsr(Env:GENOTYPE, Gu=EK_G),
#              rcov= ~ units,data=gy, verbose = FALSE)


#model fitting with mmer for GOxE and GAOxE-BLUP.
ansCS <- mmer(Yield~Env,
              random= ~ vsr(GENOTYPE_G, Gu=K_G) +vsr(GENOTYPE_tr, Gu=K_tr)+vsr(GENOTYPE_tw, Gu=K_tw)+vsr(GENOTYPE_pr, Gu=K_pr)+vsr(GENOTYPE_pw, Gu=K_pw)+
                vsr(Env:GENOTYPE_G, Gu=EK_G)+ vsr(Env:GENOTYPE_tr, Gu=EK_tr)+ vsr(Env:GENOTYPE_tw, Gu=EK_tw)+ vsr(Env:GENOTYPE_pr, Gu=EK_pr)+ vsr(Env:GENOTYPE_pw, Gu=EK_pw),
              rcov= ~ units,data=gy, verbose = FALSE)




# extracting blups for different interaction random effects. only need to run line 115 for GxE BLUP
blup <- matrix(ansCS$U$`u:Env:GENOTYPE_G`$Yield, nrow = length(unique(gy$GENOTYPE_G)),ncol = length(unique(gy$Env)))
blup2 <- matrix(ansCS$U$`u:Env:GENOTYPE_tr`$Yield, nrow = length(unique(gy$GENOTYPE_tr)),ncol = length(unique(gy$Env)))
blup3 <- matrix(ansCS$U$`u:Env:GENOTYPE_tw`$Yield, nrow = length(unique(gy$GENOTYPE_tw)),ncol = length(unique(gy$Env)))
blup4 <- matrix(ansCS$U$`u:Env:GENOTYPE_pr`$Yield, nrow = length(unique(gy$GENOTYPE_pr)),ncol = length(unique(gy$Env)))
blup5 <- matrix(ansCS$U$`u:Env:GENOTYPE_pw`$Yield, nrow = length(unique(gy$GENOTYPE_pw)),ncol = length(unique(gy$Env)))


# here we set the dimnames of different blup matrices. Only need to run lines 123 and 124 for GxE BLUP
colnames(blup) <- unique(gsub(":.*","",names(ansCS$U$`u:Env:GENOTYPE_G`$Yield)))
rownames(blup) <- unique(gsub(".*:","",names(ansCS$U$`u:Env:GENOTYPE_G`$Yield)))
dimnames(blup2) <- dimnames(blup);dimnames(blup3) <- dimnames(blup)
dimnames(blup4) <- dimnames(blup);dimnames(blup5) <- dimnames(blup)

# now we will add different interaction effects together. Don't run it for GxE BLUP
blup <- blup + blup2 + blup3 + blup4 + blup5 

## here we will add different main effects
#only for GxE-BLUP
#for(j in 1:ncol(blup)){blup[,j] <- blup[,j]+as.numeric(ansCS$U$`u:GENOTYPE_G`$Yield)}

# for GOxE and GAOxE-BLUP
for(j in 1:ncol(blup)){blup[,j] <- blup[,j]+as.numeric(ansCS$U$`u:GENOTYPE_G`$Yield)+as.numeric(ansCS$U$`u:GENOTYPE_tr`$Yield)+as.numeric(ansCS$U$`u:GENOTYPE_tw`$Yield)+as.numeric(ansCS$U$`u:GENOTYPE_pr`$Yield)+as.numeric(ansCS$U$`u:GENOTYPE_pw`$Yield)}


# here we will store the predictive abilities in the list ansCS_cors
ansCS_cors[[i]] <- diag(cor(grain_yield[which(selected_fold==i),match(colnames(blup),colnames(grain_yield))],
                            blup[which(selected_fold==i),],use = "p"))

#rm(ansCS)

#####################################
}

ansCS_cors_rf[[reps]] <- ansCS_cors
rm(ansCS_cors)
}

saveRDS(ansCS_cors_rf,file = "res.rds")


