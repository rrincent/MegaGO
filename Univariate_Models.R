#################################################################################
##### In this script we will run GBLUP, GOBLUP, OBLUP, GAOBLUP,AOBLUP ###########
#################################################################################

# 1/ Here we will load necessary functions and packaages
KinshipTransform <- function(Matrix) {
  
  nn <- ncol(Matrix)
  
  return((sum(diag(Matrix))-as.numeric(matrix(1,1,nn) %*% Matrix %*% matrix(1,nn,1))/nn)/(nn-1))
  
}

library(sommer)
library(caret)
#_______________________________________________________________________________

# 2/ Here we will load different data
#Phenotypes
grain_yield <- read.csv("~/work/megalmm/Scenario 1/grain.csv")
colnames(grain_yield) <- gsub(".grain_yield_15","",colnames(grain_yield))
grain_yield <- data.frame(grain_yield[,1],scale(grain_yield[,-1]))
colnames(grain_yield)[1] <- "GENOTYPE"

#Kinship matrices
K_G <- readRDS("kin.RDS")
K_G <- K_G[match(grain_yield$GENOTYPE,rownames(K_G)),match(grain_yield$GENOTYPE,colnames(K_G))]
K_G <- K_G/KinshipTransform(K_G) #Function for scaling the matrix

K_tr <- readRDS("kin.RDS") #standard or additive omics based
K_tr <- K_tr[match(grain_yield$GENOTYPE,rownames(K_tr)),match(grain_yield$GENOTYPE,colnames(K_tr))]
K_tr <- K_tr/KinshipTransform(K_tr)

K_tw <- readRDS("kin.RDS") #standard or additive omics based
K_tw <- K_tw[match(grain_yield$GENOTYPE,rownames(K_tw)),match(grain_yield$GENOTYPE,colnames(K_tw))]
K_tw <- K_tw/KinshipTransform(K_tw)

K_pr <- readRDS("kin_pr.RDS") #standard or additive omics based
K_pr <- K_pr[match(grain_yield$GENOTYPE,rownames(K_pr)),match(grain_yield$GENOTYPE,colnames(K_pr))]
K_pr <- K_pr/KinshipTransform(K_pr)

K_pw <- readRDS("kin_pw.RDS") #standard or additive omics based
K_pw <- K_pw[match(grain_yield$GENOTYPE,rownames(K_pw)),match(grain_yield$GENOTYPE,colnames(K_pw))]
K_pw <- K_pw/KinshipTransform(K_pw)

#____________________________________________________________________________________

# 3/ In this part, we will set up the cross validation

nreps <- 5
nfolds <- 5
set.seed(1)
folds_mat <- replicate(nreps,createFolds(rnorm(grain_yield$GENOTYPE),k=nfolds,list = F,returnTrain = F))

# 4/ select the phenotype here. you can also use cbind to bind multiple phenotypes
pheno <- grain_yield # the first column corresponds to the genotype, and the remaining are phenotypes
rownames(pheno) <- pheno$GENOTYPE
pheno <- pheno[,-1]


# 5/ here we will create empty objects to store predictive abilities 

gblup_cors <- data.frame(matrix(NA,nrow = 1,ncol = 3))
colnames(gblup_cors) <- c("Env","Pcor","Fold")

oblup_cors <- data.frame(matrix(NA,nrow = 1,ncol = 3))
colnames(oblup_cors) <- c("Env","Pcor","Fold")

goblup_cors <- data.frame(matrix(NA,nrow = 1,ncol = 3))
colnames(goblup_cors) <- c("Env","Pcor","Fold")


# 6/ Start of the loop


for (sel_rep in 1:nreps) { #to go through different reps
  
 
  folds <- folds_mat[,sel_rep] #selects the folds column to go into the nested folds loop
  
  for (i in 1:nfolds) { # to go thorugh different folds
    
    pheno2 <- pheno # here we create pheno2 to mask test phenotypes
    
    pheno2[which(folds==i),] <- NA
    
    for(j in 1:ncol(pheno2)){
      
      sel_data <- data.frame(rownames(pheno2),rownames(pheno2),rownames(pheno2),rownames(pheno2),rownames(pheno2),pheno2[,j])
      
      colnames(sel_data) <- c("GID1","GID2","GID3","GID4","GID5","Y") # multiple ID columns to be fed to different random terms
      
      ### GBLUP
      mod_gblup <- mmer(Y~1, random = ~vsr(GID1,Gu=K_G), rcov = ~units, data = sel_data,verbose = F)
      
      mod_gblup$U$`u:GID1`$Y <- mod_gblup$U$`u:GID1`$Y[match(sel_data$GID1,names(mod_gblup$U$`u:GID1`$Y))] #sommer tends to reorder predicted phenotypes
      
      
      preds <- data.frame(sel_data$GID1,pheno[,j],mod_gblup$U$`u:GID1`$Y)
      
      colnames(preds) <- c("ID","obs","pred")
      
      gblup <- data.frame(matrix(NA,nrow = 1,ncol = 3))
      colnames(gblup) <- c("Env","Pcor","Fold")
      
     
      gblup$Env <- colnames(pheno)[j]
 
      
      gblup$Pcor <- cor(preds$obs,preds$pred,use="p")
      
      gblup$Fold <- paste0("Rep_",sel_rep,"Fold_",i)
      
      gblup_cors <- rbind(gblup_cors,gblup)
      
      rm(mod_gblup,gblup,preds)
      ###########################
      
      #### goblup
      
      mod_goblup <- mmer(Y~1, random = ~vsr(GID1,Gu=K_G) + vsr(GID2,Gu=K_tr) + vsr(GID3,Gu=K_tw) + vsr(GID4,Gu=K_pr) + vsr(GID5,K_pw), rcov = ~units, data = sel_data,verbose = F)
      
      mod_goblup_res <- mod_goblup$U$`u:GID1`$Y +  mod_goblup$U$`u:GID2`$Y +  mod_goblup$U$`u:GID3`$Y +  mod_goblup$U$`u:GID4`$Y +  mod_goblup$U$`u:GID5`$Y
      
      mod_goblup_res <- mod_goblup_res[match(sel_data$GID1,names(mod_goblup_res))] #sommer tends to reorder the predicted phenotypes
      
      preds <- data.frame(sel_data$GID1,pheno[,j],as.numeric(mod_goblup_res))
      
      colnames(preds) <- c("ID","obs","pred") # here are the observed and the predicted phenotypes
      
      goblup <- data.frame(matrix(NA,nrow = 1,ncol = 3))
      colnames(goblup) <- c("Env","Pcor","Fold")
      
      
      goblup$Env <- colnames(pheno)[j]
      
      
      goblup$Pcor <- cor(preds$obs,preds$pred,use="p")
      
      goblup$Fold <- paste0("Rep_",sel_rep,"Fold_",i)
      
      goblup_cors <- rbind(goblup_cors,goblup)
      
      rm(mod_goblup,goblup,preds)
      ################################################
      #### oblup
      mod_oblup <- mmer(Y~1, random = ~vsr(GID2,Gu=K_tr) + vsr(GID3,Gu=K_tw) + vsr(GID4,Gu=K_pr) + vsr(GID5,K_pw), rcov = ~units, data = sel_data,verbose = F)
      
      mod_oblup_res <- mod_oblup$U$`u:GID2`$Y +  mod_oblup$U$`u:GID3`$Y +  mod_oblup$U$`u:GID4`$Y +  mod_oblup$U$`u:GID5`$Y
      
      mod_oblup_res <- mod_oblup_res[match(sel_data$GID2,names(mod_oblup_res))]
      
      preds <- data.frame(sel_data$GID2,pheno[,j],as.numeric(mod_oblup_res))
      
      colnames(preds) <- c("ID","obs","pred")
      
      oblup <- data.frame(matrix(NA,nrow = 1,ncol = 3))
      colnames(oblup) <- c("Env","Pcor","Fold")
      
      
      oblup$Env <- colnames(pheno)[j]
      
      
      oblup$Pcor <- cor(preds$obs,preds$pred,use="p")
      
      oblup$Fold <- paste0("Rep_",sel_rep,"Fold_",i)
      
      oblup_cors <- rbind(oblup_cors,oblup)
      
      rm(mod_oblup,oblup,preds)
      
      print(paste0("Fold:",i," Env:",j))
      
    }
    
    
  }
 

}

gblup_cors <- gblup_cors[-1,];goblup_cors <- gbolup_cors[-1,];oblup_cors <- oblup_cors[-1,] # we will remove the NA rows we initially created

saveRDS(gblup_cors,"gblup_cors.rds");saveRDS(goblup_cors,"goblup_cors.rds");saveRDS(oblup_cors,"oblup_cors.rds")


#end
