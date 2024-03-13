################################################################################
#### In this script we will extract the additive part of omics expressions ####
################################################################################

# 1/ load the package
library(rrBLUP)

# 2/ load the data

transc_wd <- readRDS("twd.rds");transc_ww <- readRDS("tww.rds") # Transcript expression data. geno ids in rows and omics in columns
spwd <- readRDS("spwd.rds");spww <- readRDS("spww.rds") # protein expression data. geno ids in rows and omics in columns
K <- readRDS("kin.rds") # genomic kinship matrix. It has the same number of lines as in T and P data
h2_prot <- readRDS("h2_prot.rds") # contains heritabilities of adjusted means of proteins expressions for WD and WW
herit_wd <- readRDS("herit_wd.rds") # contains h2 for T WD
herit_wd <- readRDS("herit_ww.rds") # contains h2 for T WW

# 3/ extracting additive part from omics
# ww transcripts
transc_ww_preds <- matrix(NA,nrow(Tww_com),ncol(Tww_com),dimnames = dimnames(Tww_com)) #an object to store the results

for (j in 1:ncol(Tww_com)) {
  transc_ww_preds[,j] <- mixed.solve(y=Tww_com[,j],
                                     X=NULL, K=K)$u
}

#wd transcripts
transc_wd_preds <- matrix(NA,nrow(Twd_com),ncol(Twd_com),dimnames = dimnames(Twd_com))


for (j in 1:ncol(Twd_com)) {
  transc_wd_preds[,j] <- mixed.solve(y=Twd_com[,j],
                                     X=NULL, K=K)$u
}

# ww proteins
spww_preds <- matrix(NA,nrow(spww),ncol(spww),dimnames = dimnames(spww))
dim(spww_preds)

for (j in 1:ncol(spww)) {
  spww_preds[,j] <- mixed.solve(y=spww[,j],
                                X=NULL, K=K)$u
}

# wd proteins
spwd_preds <- matrix(NA,nrow(spwd),ncol(spwd),dimnames = dimnames(spwd))
dim(spwd_preds)

for (j in 1:ncol(spwd)) {
  spwd_preds[,j] <- mixed.solve(y=spwd[,j],
                                X=NULL, K=K)$u
}

# 4/ here we will estimate kinship matrices from the additive part of the omics data
# first let's filter the data based on heritabilities

twd_04 <- transc_wd_preds[,colnames(transc_wd_preds) %in% herit_wd$transcript[herit_wd$value >= 0.4]]
tww_04 <- transc_ww_preds[,colnames(transc_ww_preds) %in% herit_ww$transcript[herit_ww$value >= 0.4]]

spwd_04 <- spwd_preds[,colnames(spwd_preds) %in% h2_prot$group[h2_prot$h2_WD >=0.4]]
spww_04 <- spww_preds[,colnames(spww_preds) %in% h2_prot$group[h2_prot$h2_WW >=0.4]]

spwd_04 <- scale(spwd_04);spww_04 <- scale(spww_04);tww_04 <- scale(tww_04);twd_04 <- scale(twd_04)


# estimate kinship for T 
K_tr_ad <- (twd_04 %*% t(twd_04))/ncol(twd_04)
diag(K_Tr_ad) <- diag(K_tr_ad) + 0.001 #incase the matrix is not positive semi.indefinite 
K_tr_ad <- K_tr_ad/KinshipTransform(K_tr_ad)

K_tw_ad <- (tww_04 %*% t(tww_04))/ncol(tww_04)
diag(K_tw_ad) <- diag(K_tw_ad) + 0.001
K_tw_ad <- K_tw_ad/KinshipTransform(K_tw_ad)

# estimate kinship for P

K_pr_ad <- (spwd_04 %*% t(spwd_04))/ncol(spwd_04)
diag(K_pr_ad) <- diag(K_pr_ad) + 0.001 #incase the matrix is not positive semi.indefinite 
K_pr_ad <- K_pr_ad/KinshipTransform(K_pr_ad)

K_pw_ad <- (spww_04 %*% t(spww_04))/ncol(spww_04)
diag(K_pw_ad) <- diag(K_pw_ad) + 0.001
K_pw_ad <- K_pw_ad/KinshipTransform(K_pw_ad)


# 5/ let's save the data

saveRDS(K_tr_ad,file="K_tr_ad.rds");saveRDS(K_tw_ad,file="K_tw_ad.rds")
saveRDS(K_pr_ad,file="K_pr_ad.rds");saveRDS(K_pw_ad,file="K_pw_ad.rds")


############end####################################


