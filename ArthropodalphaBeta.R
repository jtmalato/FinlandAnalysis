# Name: Rui Miguel Carvalho
# Date of creation: 10/13/2012
# Date of last update: 


# Index ---------------------------------------------------

# Datasets utilizados:
#   1.  GLMM_Variables.csv       - data.variables
#   2.  Control250_Fin.csv       - data.alpha.controls
#   3.  Traits_All_Fin.csv       - data.traits.all.fin
#   4.  Traits_Nat_Fin.csv       - data.traits.nat.fin
#   5.  Traits_NInd_Fin.csv      - data.traits.nind.fin
#   6.  All_Fin.csv              - data.all.fin
#   7.  Nat_Fin.csv              - data.nat.fin
#   8.  NInd_Fin.csv             - data.nind.fin
#   9.  MDS_vectors.csv          - data.mds
#   10. Weight_Ratios_Traits.csv - data.wr

# CALC COMPLETENESS - Faz as curvas de acumulação para os controlos 250
# SETTTING FILES
# CALC_COMPLETENESS
# CALC_TD_ALPHA
# CALC_TD_BETA
# CALC_FD_ALPHA
# CALC_FD_BETA
# Summing up Alpha values - 162
# Summing up Beta values - 336
# RESULTS - 427
# Data Exploration - 451
# GLMM - 478


# Libraries -----------------------------------------------
library(BAT)
library(readr)
library(FD) 
library(alphahull)
library(hypervolume)
library(car)
library(MASS)
library(lme4)
# source(file = "HighstatLibV10.R") #cool tools to support
# library(factoextra) # Useful for PCA analysis

# jtm added these packages
library(here)       # to locate files
library(data.table) # to work with data
library(dplyr)      # to manage data
library(magrittr)   # to use the pipe operator %>% 


# Load files ----------------------------------------------

# 1. data.variables ----------------------
# Ficheiro com as variável de distâncias com edge, trilhos e etc - para o GLMM

data.variables <- fread(here("data", "GLMM_Variables.csv"), stringsAsFactors = T)
colnames(data.variables)
str(data.variables)
summary(data.variables)
# first look at the data.variables dataset
barplot(table(data.variables$ForestID), las=2, main="Forest ID", ylab="Frequency")
barplot(table(data.variables$Treatment), las=2, main="Treatment", ylab="Frequency")

# 2. data.alpha.controls -----------------
# Ficheiro com as abundâncias por amostra para os controlos 250/Max
data.alpha.controls <- fread(here("data", "Control250_Fin.csv"), header = T)
# change column name and transform in in a factor
setnames(data.alpha.controls, old="Trail_Sampling Area", new="TSA")
data.alpha.controls[, TSA:=factor(TSA)]
summary(data.alpha.controls)


# 3. data.traits.all.fin -----------------
# Matrix with species (rows)x traits (cols) for all species in all plots
data.traits.all.fin  <- read.csv2(here("data", "Traits_All_Fin.csv"), header=TRUE, dec=".")
setDT(data.traits.all.fin)
# edit data.traits.all.fin colnames
setnames(data.traits.all.fin,
         old = c("MF","SPECIES","FAMILY","FAMILY.1","N.E.I","IND.NONNAIND","HABITATNALeaves..branches..flowers..seeds..fruits..surface.","HABITATNAinside.stems..roots..fruits..pods..fungi","HABITATNAGround","HABITATNAUnder.stones..bark..twigs","HABITATNADecaying.matter","HABITATNASubterranean","Stenophagous...Not.Euryphagous"),
         new = c("Mf","Sp","Family.name","Family","Nei","Ind.nonind","Hab.leaves","Hab.inside","Hab.ground","Hab.under","Hab.matter","Hab.sub","Stenophagous"))
data.traits.all.fin
# data.traits.all.fin <- data.traits.all.fin[,c(1:3)]
str(data.traits.all.fin) # as colunas Hab. poderao ser factores 0/1
pairs(data.traits.all.fin[,c(13:18)]) # all same data, might justify the positive pattern
pairs(data.traits.all.fin[,c(19:21)], pch=16)


# 4. data.traits.nat.fin -----------------
# Matrix with species (rows)x traits (cols) for all natives species in all plots
data.traits.nat.fin <- read.csv2(here("data", "Traits_Nat_Fin.csv"), header=TRUE, dec=".")
setDT(data.traits.nat.fin)
colnames(data.traits.nat.fin)
setnames(data.traits.nat.fin,
         old = c("MF","SPECIES","FAMILY","FAMILY.1","N.E.I","IND.NONNAIND","HABITATNALeaves..branches..flowers..seeds..fruits..surface.","HABITATNAinside.stems..roots..fruits..pods..fungi","HABITATNAGround","HABITATNAUnder.stones..bark..twigs","HABITATNADecaying.matter","HABITATNASubterranean","Stenophagous...Not.Euryphagous"),
         new = c("Mf","Sp","Family.name","Family","Nei","Ind.nonind","Hab.leaves","Hab.inside","Hab.ground","Hab.under","Hab.matter","Hab.sub","Stenophagous"))
# mudar os valores de "E " para "E" (sem espaco)
levels(data.traits.nat.fin$Nei) <- c("E", "N")
# ver a estrutura da coluna
as.character(data.traits.nat.fin$Nei)
# estrutura dos dados
str(data.traits.nat.fin)


# 5. data.traits.nind.fin ----------------
# Matrix with species (rows)x traits (cols) for all Non-Indigenous species in all plots
data.traits.nind.fin  <- read.csv2(here("data", "Traits_NInd_Fin.csv"), header=TRUE, dec = ".")
setDT(data.traits.nind.fin)
colnames(data.traits.nind.fin)
setnames(data.traits.nind.fin,
         old = c("MF","SPECIES","FAMILY","FAMILY.1","N.E.I","IND.NONNAIND","HABITATNALeaves..branches..flowers..seeds..fruits..surface.","HABITATNAinside.stems..roots..fruits..pods..fungi","HABITATNAGround","HABITATNAUnder.stones..bark..twigs","HABITATNADecaying.matter","HABITATNASubterranean","Stenophagous...Not.Euryphagous"),
         new = c("Mf","Sp","Family.name","Family","Nei","Ind.nonind","Hab.leaves","Hab.inside","Hab.ground","Hab.under","Hab.matter","Hab.sub","Stenophagous"))
str(data.traits.nind.fin)
summary(data.traits.nind.fin) # há duas linhas NA


# 6. data.all.fin ------------------------
# Ficheiro com as abundâncias por Área de amostragem, para todas as amostras
SAAll <- fread(here("data", "All_Fin.csv"), header = T) %>% 
  setnames(old="Trail_Sampling Area", new="TSA")
sites <- as.vector(SAAll$TSA)
# SAAll <- SAAll[,-1]
# rownames(SAAll) <- sites
# colnames(SAAll) <- species

# HEllinger transformation
SAAll.standard <- decostand(SAAll[,-1], "hellinger")
SAAll.standard


# 7. data.nat.fin ------------------------
# Ficheiro com as abundâncias por área de amostragem, para as espécies indígenas
SANat <- fread(here("data", "Nat_Fin.csv"), header = T) %>% 
  setnames(old="Trail_Sampling Area", new="TSA")
# sites are the same as the previous object
# SANat <- SANat[,-1] 
# rownames(SANat) <- sites
# colnames(SANat) <- species_nat
# HEllinger transformation
SANat.standard <- decostand(SANat[,-1], "hellinger") 
SANat.standard


# 8. data.nind.fin -----------------------
# Ficheiro com as abundâncias por área de amostragem, para as espécies não-indígenas
SANInd <- fread(here("data", "NInd_Fin.csv"), header = T) %>% 
  setnames(old="Trail_Sampling Area", new="TSA")
# SANInd <- SANInd[,-1] 
# rownames(SANInd) <- sites
# colnames(SANInd) <- species_nind
# HEllinger transformation
SANInd.standard <- decostand(SANInd[,-1], "hellinger") 
SANInd.standard


# 9. data.mds ----------------------------
# Vector with weight trail and treatment names for MDS
data.mds <- fread(here("data", "MDS_vectors.csv"))
trail <- unique(as.vector(data.mds$Trail_name))
str(data.mds)


# 10. data.wr -----------------------------
# Vector with weight ratios for trails - used in cluster::daisy to compensate decomposition of one variable by
# multiple columns
data.wr <- read.csv2(here("data", "Weight_Ratios_Traits.csv"), header = TRUE, dec = ".")
setDT(data.wr)
colnames(data.wr)
setnames(data.wr,
         old=c("MF","SPECIES","FAMILY","FAMILY.1","N.E.I","IND.NON.IND","HABITAT.Leaves..branches..flowers..seeds..fruits..surface.","HABITAT.inside.stems..roots..fruits..pods..fungi","HABITAT.Ground","HABITAT.Under.stones..bark..twigs","HABITAT.Decaying.matter","HABITAT.Subterranean","Stenophagous...Not.Euryphagous"),
         new=c("Mf","Sp","Family.name","Family","Nei","Ind.nonind","Hab.leaves","Hab.inside","Hab.ground","Hab.under","Hab.matter","Hab.sub","Stenophagous"))
# colnames(data.wr)
# weights <- data.wr[,-c(1:4)]
# weights <- as.vector(data.wr)


# Looking into the datasets -------------------------------

# CALCULATING COMPLETENESS FOR THE CONTROLS FOR BETA ANALYSYS -------
alphaaccumlist <- list()

tsa <- levels(data.alpha.controls$TSA)
# pdf(here("results", "accumCurvesControls.pdf"), width = 10, height = 8)
par(mfrow = c(2, 2))
for(i in seq_along(tsa)) {
  
  # definir a matriz a analisar
  data.alpha.controls.matrix <- data.alpha.controls[which(data.alpha.controls$TSA == tsa[i]), -1]
  alphaaccumlist[[i]] <- alpha.accum(data.alpha.controls.matrix, prog = F)

  # plot(alphaaccumlist[[i]][,1], alphaaccumlist[[i]][,2], xlab="Samples", ylab="Individuals", )
  plot(alphaaccumlist[[i]][, 1], alphaaccumlist[[i]][, 3], xlab = "Samples", ylab = "Individuals",
       ylim = c(0,40), xlim = c(0,25), main = tsa[i], las =1, type="n")
  lines(alphaaccumlist[[i]][, 1], alphaaccumlist[[i]][, 4], col="gray30", lwd=2)
  lines(alphaaccumlist[[i]][, 1], alphaaccumlist[[i]][, 5], col = "#C8631B", lwd=2)
  lines(alphaaccumlist[[i]][, 1], alphaaccumlist[[i]][, 8], col = "#BF020A", lwd=2)
  lines(alphaaccumlist[[i]][, 1], alphaaccumlist[[i]][, 17], col = "#0ABF02", lwd=2)
  points(alphaaccumlist[[i]][, 1], alphaaccumlist[[i]][, 3], pch=21, bg = "#1B80C8", cex=1.2)
  points(alphaaccumlist[[i]][, 1], alphaaccumlist[[i]][, 13], pch=21, bg = "#BF020A", cex=1.2)
  points(alphaaccumlist[[i]][, 1], alphaaccumlist[[i]][, 19], pch=21, bg = "#0ABF02", cex=1.2)
}
# dev.off()


# CALC_COMPLETENESS -----

# Quero usar o fichero dos contrlos para calcular os estimadores, e daí obter completeness.
# rownames(data.alpha.controls) <-tsa
# colnames(data.alpha.controls) <- tsa
# data.alpha.controls <- data.alpha.controls[,-1]
alpha.estimate(data.alpha.controls[,-1])


# Calculating functional elements -----

# All species tree and FD Alpha
# Calculating distances between species
dissimAll <- cluster::daisy(data.traits.all.fin, metric = "gower", weights = data.wr)
# https://www.rdocumentation.org/packages/cluster/versions/2.0.7-1/topics/daisy

# o warning diz que estas colunas binarias vao ser tratadas como um scaled interval
# data.traits.all.fin[, c(7, 9, 10, 11, 12,22, 24, 26, 27, 28, 29, 31, 31)]


treeAll <- hclust(dissimAll, "average")
# pdf(here("results", "dissimAll.dendrogram.pdf"), width = 10, height = 8)
par(mfrow=c(1,1), mar=c(1, 4, 3, 1))
plot(treeAll, hang = -1, las=1, xlab="", main="Dendrogram of \nAll_Fin Species")
# dev.off()
# plot(as.dendrogram(tree))


# All natives tree and FD Alpha ----
# Calculating distances between species
dissimNat <- cluster::daisy(data.traits.nat.fin, "gower", weights = data.wr)
treeNat <- hclust(dissimNat, "average") # building the dendrogrma
cor(dissimNat, cophenetic(treeNat))

# pdf(here("results", "dissimNat.dendrogram.pdf"), width = 10, height = 8)
par(mfrow=c(1,1), mar=c(1, 4, 5, 1))
plot(treeNat, hang = -1, las=1, xlab="", main="Dendrogram of \nNat_Fin Species")
# dev.off()
# plot(as.dendrogram(treeNat))


# All Non-Ind tree and FD Alpha ----
# Calculating distances between species
dissimNInd <- cluster::daisy(data.traits.nind.fin, "gower", weights = data.wr)
treeNInd <- hclust(dissimNInd, "average")
# pdf(here("results", "dissimNInd.dendrogram.pdf"), width = 10, height = 8)
par(mfrow=c(1,1), mar=c(1, 4, 5, 1))
plot(treeNInd, hang = -1, las=1, xlab="", main="Dendrogram of \nNInd_Fin Species")
# dev.off()
# plot(as.dendrogram(treeNInd))


# ALPHAs for All species, Natives and Non-Indigenous ----

# Taxonomical Alpha
alpha.All  <- alpha(SAAll[,-1])
alpha.Nat  <- alpha(SANat[, -1])
alpha.NInd <- alpha(SANInd[, -1])

# Functional Alpha
# Functional Alpha All devolve um erro!!
FDalpha.All  <- alpha(SAAll[,-1], treeAll)
FDalpha.Nat  <- alpha(SANat[,-1], treeNat)
FDalpha.NInd <- alpha(SANInd[,-1], treeNInd)

alphas <- data.table(alpha.All=alpha.All, alpha.Nat=alpha.Nat, alpha.NInd=alpha.NInd, FDalpha.All=NA, FDalpha.Nat=FDalpha.Nat, FDalpha.NInd=FDalpha.NInd)
colnames(Alphas) <- cbind("TAlphaAll", "TAlphaNat", "TAlphaNInd","FAlphaAll", "FAlphaNat", "FAlphaNInd" )
Alphas # Will be printed along with BETAS in the RESULTS file

# Calculating Betas for All species, Natives and Non-Indigenous ----
# Beta results come as a list (class(BetaAll) = list), so we transform them into data frame to export them into csv

# Taxonomical Beta
beta.All <- beta(SAAll[,-1])
beta.Nat <- beta(SANat[,-1])
beta.Nat <- beta(SANInd[,-1])

# Functional Beta
# Functional beta All devolve um erro 
betaFunc.All <- beta(SAAll[,-1], treeAll, abund =  T)
betaFunc.Nat <- beta(SANat[,-1], treeNat, abund = T) 
betaFunc.NInd <- beta(SANInd[,-1], treeNInd, abund= T)


## Separating Total, Richness and Replacement  TAXONOMICAL betas into different data frames
BetaAllTotal  <- data.table(as.matrix(beta.All[["Btotal"]]))
BetaAllRich   <- data.table(as.matrix(beta.All[["Brich"]]))
BetaAllRepl   <- data.table(as.matrix(beta.All[["Brepl"]]))
BetaNatTotal  <- data.table(as.matrix(beta.Nat[["Btotal"]]))
BetaNatRich   <- data.table(as.matrix(beta.Nat[["Brich"]]))
BetaNatRepl   <- data.table(as.matrix(beta.Nat[["Brepl"]]))
BetaNIndTotal <- data.table(as.matrix(beta.Nat[["Btotal"]]))
BetaNIndRich  <- data.table(as.matrix(beta.Nat[["Brich"]]))
BetaNIndRepl  <- data.table(as.matrix(beta.Nat[["Brepl"]]))

## Separating Total, Richness and Replacement  FUNCTIONAL betas into different data frames
BetaFuncAllTotal  <- data.table(as.matrix(beta.FuncAll[["Btotal"]]))
BetaFuncAllRich   <- data.table(as.matrix(beta.FuncAll[["Brich"]]))
BetaFuncAllRepl   <- data.table(as.matrix(beta.FuncAll[["Brepl"]]))
BetaFuncNatTotal  <- data.table(as.matrix(betaFunc.Nat[["Btotal"]]))
BetaFuncNatRich   <- data.table(as.matrix(betaFunc.Nat[["Brich"]]))
BetaFuncNatRepl   <- data.table(as.matrix(betaFunc.Nat[["Brepl"]]))
BetaFuncNIndTotal <- data.table(as.matrix(betaFunc.Nat[["Btotal"]]))
BetaFuncNIndRich  <- data.table(as.matrix(betaFunc.Nat[["Brich"]]))
BetaFuncNIndRepl  <- data.table(as.matrix(betaFunc.Nat[["Brepl"]]))

# Separating the TAXONOMICAL beta values between the Control 250/Max and the other sampling areas from each trail ----

A01 <- BetaAllTotal[1,1]
AA1 <- BetaAllTotal[6,c(2:6)]
BB1 <- BetaAllTotal[12,c(7:12)]
CC1 <- BetaAllTotal[18,c(13:18)]
DD1 <- BetaAllTotal[22,c(19:22)]
BetaAllTotalVector <- c(A01,AA1, BB1,  CC1, DD1)

A02 <- BetaAllRich[1,1]
AA2 <- BetaAllRich[6,c(2:6)]
BB2 <- BetaAllRich[12,c(7:12)]
CC2 <- BetaAllRich[18,c(13:18)]
DD2<- BetaAllRich[22,c(19:22)]
BetaAllRichVector <- c(A02,AA2, BB2,  CC2, DD2)

A03 <- BetaAllRepl[1,1]
AA3 <- BetaAllRepl[6,c(2:6)]
BB3 <- BetaAllRepl[12,c(7:12)]
CC3 <- BetaAllRepl[18,c(13:18)]
DD3 <- BetaAllRepl[22,c(19:22)]
BetaAllReplVector <- c(A03,AA3, BB3, CC3, DD3)

A01 <- BetaNatTotal[1,1]
AA1 <- BetaNatTotal[6,c(2:6)]
BB1 <- BetaNatTotal[12,c(7:12)]
CC1 <- BetaNatTotal[18,c(13:18)]
DD1 <- BetaNatTotal[22,c(19:22)]
BetaNatTotalVector <- c(A01,AA1, BB1,  CC1, DD1)

A02 <- BetaNatRich[1,1]
AA2 <- BetaNatRich[6,c(2:6)]
BB2 <- BetaNatRich[12,c(7:12)]
CC2 <- BetaNatRich[18,c(13:18)]
DD2<- BetaNatRich[22,c(19:22)]
BetaNatRichVector <- c(A02,AA2, BB2,  CC2, DD2)

A03 <- BetaNatRepl[1,1]
AA3 <- BetaNatRepl[6,c(2:6)]
BB3 <- BetaNatRepl[12,c(7:12)]
CC3 <- BetaNatRepl[18,c(13:18)]
DD3 <- BetaNatRepl[22,c(19:22)]
BetaNatReplVector <- c(A03,AA3, BB3, CC3, DD3)

A01 <- BetaNIndTotal[1,1]
AA1 <- BetaNIndTotal[6,c(2:6)]
BB1 <- BetaNIndTotal[12,c(7:12)]
CC1 <- BetaNIndTotal[18,c(13:18)]
DD1 <- BetaNIndTotal[22,c(19:22)]
BetaNIndTotalVector <- c(A01,AA1, BB1,  CC1, DD1)

A02 <- BetaNIndRich[1,1]
AA2 <- BetaNIndRich[6,c(2:6)]
BB2 <- BetaNIndRich[12,c(7:12)]
CC2 <- BetaNIndRich[18,c(13:18)]
DD2<- BetaNIndRich[22,c(19:22)]
BetaNIndRichVector <- c(A02,AA2, BB2,  CC2, DD2)

A03 <- BetaNIndRepl[1,1]
AA3 <- BetaNIndRepl[6,c(2:6)]
BB3 <- BetaNIndRepl[12,c(7:12)]
CC3 <- BetaNIndRepl[18,c(13:18)]
DD3 <- BetaNIndRepl[22,c(19:22)]
BetaNIndReplVector <- c(A03,AA3, BB3, CC3, DD3)

## Separating the FUNCTIONAL beta values between the Control 250/Max and the other sampling areas from each trail


A01 <- BetaFuncAllTotal[1,1]
AA1 <- BetaFuncAllTotal[6,c(2:6)]
BB1 <- BetaFuncAllTotal[12,c(7:12)]
CC1 <- BetaFuncAllTotal[18,c(13:18)]
DD1 <- BetaFuncAllTotal[22,c(19:22)]
BetaFuncAllTotalVector <- c(A01,AA1, BB1,  CC1, DD1)

A02 <- BetaFuncAllRich[1,1]
AA2 <- BetaFuncAllRich[6,c(2:6)]
BB2 <- BetaFuncAllRich[12,c(7:12)]
CC2 <- BetaFuncAllRich[18,c(13:18)]
DD2<- BetaFuncAllRich[22,c(19:22)]
BetaFuncAllRichVector <- c(A02,AA2, BB2,  CC2, DD2)

A03 <- BetaFuncAllRepl[1,1]
AA3 <- BetaFuncAllRepl[6,c(2:6)]
BB3 <- BetaFuncAllRepl[12,c(7:12)]
CC3 <- BetaFuncAllRepl[18,c(13:18)]
DD3 <- BetaFuncAllRepl[22,c(19:22)]
BetaFuncAllReplVector <- c(A03,AA3, BB3, CC3, DD3)

A01 <- BetaFuncNatTotal[1,1]
AA1 <- BetaFuncNatTotal[6,c(2:6)]
BB1 <- BetaFuncNatTotal[12,c(7:12)]
CC1 <- BetaFuncNatTotal[18,c(13:18)]
DD1 <- BetaFuncNatTotal[22,c(19:22)]
BetaFuncNatTotalVector <- c(A01,AA1, BB1,  CC1, DD1)

A02 <- BetaFuncNatRich[1,1]
AA2 <- BetaFuncNatRich[6,c(2:6)]
BB2 <- BetaFuncNatRich[12,c(7:12)]
CC2 <- BetaFuncNatRich[18,c(13:18)]
DD2<- BetaFuncNatRich[22,c(19:22)]
BetaFuncNatRichVector <- c(A02,AA2, BB2,  CC2, DD2)

A03 <- BetaFuncNatRepl[1,1]
AA3 <- BetaFuncNatRepl[6,c(2:6)]
BB3 <- BetaFuncNatRepl[12,c(7:12)]
CC3 <- BetaFuncNatRepl[18,c(13:18)]
DD3 <- BetaFuncNatRepl[22,c(19:22)]
BetaFuncNatReplVector <- c(A03,AA3, BB3, CC3, DD3)

A01 <- BetaFuncNIndTotal[1,1]
AA1 <- BetaFuncNIndTotal[6,c(2:6)]
BB1 <- BetaFuncNIndTotal[12,c(7:12)]
CC1 <- BetaFuncNIndTotal[18,c(13:18)]
DD1 <- BetaFuncNIndTotal[22,c(19:22)]
BetaFuncNIndTotalVector <- c(A01,AA1, BB1,  CC1, DD1)

A02 <- BetaFuncNIndRich[1,1]
AA2 <- BetaFuncNIndRich[6,c(2:6)]
BB2 <- BetaFuncNIndRich[12,c(7:12)]
CC2 <- BetaFuncNIndRich[18,c(13:18)]
DD2<- BetaFuncNIndRich[22,c(19:22)]
BetaFuncNIndRichVector <- c(A02,AA2, BB2,  CC2, DD2)

A03 <- BetaFuncNIndRepl[1,1]
AA3 <- BetaFuncNIndRepl[6,c(2:6)]
BB3 <- BetaFuncNIndRepl[12,c(7:12)]
CC3 <- BetaFuncNIndRepl[18,c(13:18)]
DD3 <- BetaFuncNIndRepl[22,c(19:22)]
BetaFuncNIndReplVector <- c(A03,AA3, BB3, CC3, DD3)


#### COMPILING ALL BETA INFORMATION INTO A TABLE, AND EXPORTING IT TO A FILE

Betas <- as.data.frame( t(rbind(
  BetaAllTotalVector, BetaAllRichVector, BetaAllReplVector, 
  BetaFuncAllTotalVector, BetaFuncAllRichVector, BetaFuncAllReplVector,
  BetaNatTotalVector, BetaNatRichVector, BetaNatReplVector, 
  BetaFuncNatTotalVector, BetaFuncNatRichVector, BetaFuncNatReplVector,
  BetaNIndTotalVector, BetaNIndRichVector, BetaNIndReplVector,
  BetaFuncNIndTotalVector, BetaFuncNIndRichVector, BetaFuncNIndReplVector)))

str(Betas)


## Exporting TAXONOMICAL beta results to .csv (for the record - it won't be usedin analysis, as the RESULTS 
##      file compiles all the  useable results )



write.csv(BetaAllTotal,file = "C:/ArthropodsArticle/2.0 _Finland_Analysis/BetaResults/BetaAllTotal.csv")
write.csv(BetaAllRich,file = "C:/ArthropodsArticle/2.0 _Finland_Analysis/BetaResults/BetaAllRich.csv")
write.csv(BetaAllRepl,file = "C:/ArthropodsArticle/2.0 _Finland_Analysis/BetaResults/BetaAllRepl.csv")
write.csv(BetaNatTotal,file = "C:/ArthropodsArticle/2.0 _Finland_Analysis/BetaResults/BetaNatTotal.csv")
write.csv(BetaNatRich,file = "C:/ArthropodsArticle/2.0 _Finland_Analysis/BetaResults/BetaNatRich.csv")
write.csv(BetaNatRepl,file = "C:/ArthropodsArticle/2.0 _Finland_Analysis/BetaResults/BetaNatRepl.csv")
write.csv(BetaNIndTotal,file = "C:/ArthropodsArticle/2.0 _Finland_Analysis/BetaResults/BetaNIndTotal.csv")
write.csv(BetaNIndRich,file = "C:/ArthropodsArticle/2.0 _Finland_Analysis/BetaResults/BetaNIndRich.csv")
write.csv(BetaNIndRepl,file = "C:/ArthropodsArticle/2.0 _Finland_Analysis/BetaResults/BetaNIndRepl.csv")

## Exporting FUNCTIONAL beta results to .csv (for the record - it won't be usedin analysis, as the RESULTS 
#       file compiles all the  useable results )

write.csv(BetaFuncAllTotal,file = "C:/ArthropodsArticle/2.0 _Finland_Analysis/BetaFuncResults/BetaFuncAllTotal.csv")
write.csv(BetaFuncAllRich,file = "C:/ArthropodsArticle/2.0 _Finland_Analysis/BetaFuncResults/BetaFuncAllRich.csv")
write.csv(BetaFuncAllRepl,file = "C:/ArthropodsArticle/2.0 _Finland_Analysis/BetaFuncResults/BetaFuncAllRepl.csv")
write.csv(BetaFuncNatTotal,file = "C:/ArthropodsArticle/2.0 _Finland_Analysis/BetaFuncResults/BetaFuncNatTotal.csv")
write.csv(BetaFuncNatRich,file = "C:/ArthropodsArticle/2.0 _Finland_Analysis/BetaFuncResults/BetaFuncNatRich.csv")
write.csv(BetaFuncNatRepl,file = "C:/ArthropodsArticle/2.0 _Finland_Analysis/BetaFuncResults/BetaFuncNatRepl.csv")
write.csv(BetaFuncNIndTotal,file = "C:/ArthropodsArticle/2.0 _Finland_Analysis/BetaFuncResults/BetaFuncNIndTotal.csv")
write.csv(BetaFuncNIndRich,file = "C:/ArthropodsArticle/2.0 _Finland_Analysis/BetaFuncResults/BetaFuncNIndRich.csv")
write.csv(BetaFuncNIndRepl,file = "C:/ArthropodsArticle/2.0 _Finland_Analysis/BetaFuncResults/BetaFuncNIndRepl.csv")

###########################################################################
#####                            RESULTS                              #####
###########################################################################


Results <- cbind.data.frame(Variables, Alphas, Betas)
str(Results)

#MAKING THE RESULTS EXPORTABLE INTO CSV
Results <- apply(Results, 2 , as.character)

#NAMING THE TRAIL SEGMENTS
rownames(Results) <- rownames(Alphas)

#PASSING RESULTS TO FILE
write.csv(Results, file = ("C:/ArthropodsArticle/FinlandAnalysis/RESULTS.csv"), row.names = TRUE)

Results2 <- read.csv2(here("RESULTS.CSV"), header=TRUE, row.names = 1,  stringsAsFactors = T,sep = ",", dec = ".")
Results2

Results3 <- read.csv("RESULTS.csv", header = TRUE, row.names = 1, sep = ",", dec = ".")
ResultsWithoutControls <- Results3[-c(5,6,11,12,17,18),-c(1:2)]