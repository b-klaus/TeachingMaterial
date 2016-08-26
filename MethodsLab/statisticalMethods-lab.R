## ----style, echo=FALSE, results='asis'-----------------------------------
BiocStyle::latex()

## ----options, include=FALSE----------------------------------------------
options(digits=3, width=80)
opts_chunk$set(echo=TRUE,tidy=FALSE,include=TRUE,
               dev='pdf', fig.width = 6, fig.height = 3.5, comment = '  ', dpi = 300,
		cache = T, lazy.load = FALSE, background="grey93" )

## ----required packages and data, echo = TRUE, message = FALSE------------
set.seed(777)
library(TeachingDemos)
library(geneplotter)
library(ggplot2)
library(plyr)
library(LSD)
library(boot)
library(ade4)
library(DESeq2)
library(vsn)
library(gplots)
library(RColorBrewer)
library(psych)
library(car)
library(matrixStats)
library(MASS)
library(vegan)
library(locfit)
library(stringr)
library(sva)
library(limma)
library(corpcor)


ggplotRegression <- function (fit) {

ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
  geom_point() +
  stat_smooth(method = "lm", col = "coral2") +
  labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                     "Intercept =",signif(fit$coef[[1]],5 ),
                     " Slope =",signif(fit$coef[[2]], 5),
                     " P =",signif(summary(fit)$coef[2,4], 5)))
}


## ----simulateRNASeq------------------------------------------------------
dds <- makeExampleDESeqDataSet(m=8, betaSD = 2)

## add sex as a batch effect, there are two males and two females in each group
colData(dds)$sex <- factor(rep(c("m", "f"), 4))

## modify counts to add a batch effect, we add normally distributed random noise
## with mean 2 to randomly selected genes of male samples and then round the result
cts <- counts(dds)
ranGenes <- floor(runif(300)*1000)

for(i in ranGenes){
  cts[i, colData(dds)$sex == "m"] <- as.integer(cts[i, colData(dds)$sex == "m"]
                                                + round(rnorm(1,4, sd = 1)))
  }
counts(dds) <- cts

counts(dds)[5:10,]

## ----meanSd_Sim_Data, dependson="simulateRNASeq", echo = TRUE, fig.show='hide'----
pl <- meanSdPlot(counts(dds))

## ----meanSd_Sim_Data_pl--------------------------------------------------
pl$gg + ylim(0,100)

## ----meanSdrlog,  dependson="meanSd_Sim_Data",  echo = TRUE--------------
rldSim <- assay(rlogTransformation(dds, blind=TRUE))
meanSdPlot(rldSim)

## ----PCA_simData, dependson="meanSdrlog"---------------------------------

ntop = 500

pvars <- rowVars(rldSim)
select <- order(pvars, decreasing = TRUE)[seq_len(min(ntop,
                                                      length(pvars)))]

PCA <- prcomp(t(rldSim)[, select], scale = F)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)


dataGG = data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                    PC3 = PCA$x[,3], PC4 = PCA$x[,4],
                    sex = colData(dds)$sex,
                    condition = colData(dds)$condition)

(qplot(PC1, PC2, data = dataGG, color =  condition, shape = sex,
       main = "PC1 vs PC2, top variable genes", size = I(6))
+ labs(x = paste0("PC1, VarExp:", round(percentVar[1],4)),
       y = paste0("PC2, VarExp:", round(percentVar[2],4)))
+ scale_colour_brewer(type="qual", palette=2)
)



## ----SimData_heatmap_and_clustering, fig.height = 5, dependson="meanSdrlog"----
dists <- as.matrix(dist(t(rldSim)))

rownames(dists) <-  paste0(colData(dds)$sex,
                                  "_", str_sub(rownames(colData(dds)), 7))
colnames(dists) <- paste0(colData(dds)$condition,
                                  "_", str_sub(rownames(colData(dds)), 7))

hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(dists, trace="none", col = rev(hmcol))

## ----r scalingSingleCell-------------------------------------------------
scData <- t(read.csv("http://www-huber.embl.de/users/klaus/nbt.3102-S7.csv",
            row.names  = 1))


scData[1:10,1:10]

distEucSC <- dist(t(scale(scData)))

scalingEuSC <- as.data.frame(cmdscale(distEucSC, k = 2))
names(scalingEuSC) <- c("MDS_Dimension_1", "MDS_Dimension_2")
head(scalingEuSC)

## ----plotScalingSingleCell, dependson="scalingSingleCell"----------------

rownames(scData)[727]

gata3 <- cut(scale(scData[727,]), 5,
             labels = c("very low", "low", "medium", "high", 'very high'))

# reverse label ordering
gata3 <- factor(gata3, rev(levels(gata3)))

scPlot <- qplot(MDS_Dimension_1, MDS_Dimension_2, data = scalingEuSC,
      main = "Metric MDS",
      color = gata3, size = I(3))  + scale_color_brewer(palette = "RdBu")


scPlot
# + scale_color_gradient2(mid = "#ffffb3")
# + scale_colour_manual(values = rev(brewer.pal(5, "RdBu"))))
gata3Colours <- invisible(print(scPlot))$data[[1]]$colour



## ----checkScaling, dependson="scalingSingleCell"-------------------------

distMDS <- dist(scalingEuSC)
## get stress
sum((distEucSC - distMDS)^2)/sum(distEucSC^2)


ord <- order(as.vector(distEucSC))
dataGG <- data.frame(dOrg = as.vector(distEucSC)[ord],
                     dMDS = as.vector(distMDS)[ord])

(qplot(dOrg, dMDS, data=dataGG,
       main = "Shepard plot: original vs MDS distances") + geom_smooth())

## ----Kruskal, dependson="scalingSingleCell"------------------------------

kruskalMDS <- isoMDS(distEucSC, y = as.matrix(scalingEuSC))

#tressplot(kruskalMDS, distEucSC)

kruskalMDS$stress


ord <- order(as.vector(distEucSC))
dataGGK <- data.frame(dOrg = as.vector(distEucSC)[ord],
                     dMDS = as.vector(dist(scores(kruskalMDS)))[ord])

(qplot(dOrg, dMDS, data=dataGGK,
       main = "Shepard plot for Kruskal: original vs MDS distances") + geom_smooth())


## ----kruskalMDSPlot, dependson=c("Kruskal", "plotScalingSingleCell")-----

scalingK <- as.data.frame(scores(kruskalMDS))
names(scalingK) <- c("MDS_Dimension_1", "MDS_Dimension_2")

qplot(MDS_Dimension_1, MDS_Dimension_2, data = scalingK,
      main = "Non--metric MDS",
      color = gata3, size = I(3)) + scale_color_brewer(palette = "PuOr")


## ----loadBodyfat,   echo = TRUE------------------------------------------
load(url("http://www-huber.embl.de/users/klaus/BasicR/bodyfat.rda"))
dim   (bodyfat)    # how many rows and columns in the dataset?
names (bodyfat)    # names of the columns)

## ----age-summary,   echo = TRUE------------------------------------------
## compute descriptive statistics for "age"
summary(bodyfat$age)
sd(bodyfat$age)
mean(bodyfat$age)
IQR(bodyfat$age)/1.349
## mean value of every variable in the bodyfat data set
sapply(bodyfat,  FUN = mean)

## ----InspectionPlots,   echo = TRUE, eval = TRUE, fig.width = 15, fig.height = 12----
if("bodyfat" %in% search()) detach(bodyfat)
bodyfat <- bodyfat[-39,]



PCA <- prcomp(bodyfat, center = TRUE, scale = TRUE)
biplot(PCA)

## plot original PCs
dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], no = seq_along(PCA$x[,1]))
qplot(PC1, PC2, data = dataGG, label = no, geom = "text")


heatpairs(as.matrix(bodyfat[,3:5]))


## ----frenchPCA,   echo = TRUE, eval = TRUE, fig.width = 8----------------
frenchPCA <- dudi.pca(bodyfat, scannf = FALSE)
frenchPCA
# biplot
scatter(frenchPCA, type ="lines")
# PCA plot
s.label(frenchPCA$l1)
# correlation circle
s.corcircle(frenchPCA$co)

## ----select_variable, echo = TRUE----------------------------------------
### look at correlations of bodyfat with other variables
select.vars <- abs(cor(bodyfat))["percent.fat" ,]
select.vars
select.vars[select.vars > 0.6]

## ----fit_model,   echo = TRUE--------------------------------------------
lm.fat <- lm(bodyfat$percent.fat ~ bodyfat$abdomen.circum )


ggplotRegression(lm.fat)

## ----plot_selected_model,   echo = TRUE----------------------------------
qplot(.fitted, .resid, data = fortify(lm.fat)) +
  geom_hline(yintercept = 0) + geom_smooth(se = FALSE)

## ----loessExampleLinFit, dependson="fit_model"---------------------------

# create_data
y <- seq(from=1, to=10, length.out=100)
a <- y^3 +y^2  + rnorm(100,mean=0, sd=30)
dataL <- data.frame(a=a, y=y)
qplot(y, a, data = dataL)

# linear fit
linreg <- lm(y~a, data = dataL)

ggplotRegression(linreg)

dataL$LinReg <- predict(linreg)

## ----loessExampleFit, dependson="loessExampleLinFit"---------------------

dataL$locFit  <- predict(locfit(y~lp(a, nn=0.5, deg=1), data=dataL),
                         newdata = dataL$a)


 (qplot(a, y, data = dataL, main = "Linear vs. local regression")
 +  geom_line(aes(x = a, y = locFit), color = "dodgerblue3")
 +  geom_line(aes(x = a, y = LinReg), color = "coral3"))


## ----getBottomly Data, echo = TRUE---------------------------------------
load(url("http://www-huber.embl.de/users/klaus/bottomly_eset.RData"))
bottomly.eset
pData(bottomly.eset)
idx.nn <- apply(exprs(bottomly.eset), 1, function(x) { all(x > 5)})
bottomly.eset <- subset(bottomly.eset, idx.nn)

## ----normalize_Bottomly_Data, echo = TRUE--------------------------------
bottomly.sf <- estimateSizeFactorsForMatrix(exprs(bottomly.eset))


bottomly.norm <- bottomly.eset

 for(i in seq_along(bottomly.sf)){
   exprs(bottomly.norm)[,i] <-  exprs(bottomly.eset)[,i] /  bottomly.sf[i]
 }


multidensity( exprs(bottomly.eset), xlab="mean counts", xlim=c(0, 500),
           main = "Bottomly raw")
multidensity( exprs(bottomly.norm), xlab="mean counts", xlim=c(0, 500),
           main = "Bottomly normalized")




## ----log2 Trans Bottomly Data, echo = TRUE, eval=TRUE--------------------
bottomly.log2 <- log2(exprs(bottomly.norm))
bottomly.log2.raw <- log2(exprs(bottomly.eset))

## ----MA plots  Bottomly Data, echo = TRUE, eval=FALSE--------------------
## pdf("pairwiseMAsBottomly.pdf")
##  MA.idx = t(combn(1:21, 2))
## 
## 	for( i in  1:dim(MA.idx)[1]){
## 	 MDPlot(bottomly.log2,
## 		c(MA.idx[i,1],MA.idx[i,2]),
## 	main = paste( sampleNames(bottomly.norm)[MA.idx[i,1]], " vs ",
## 	 sampleNames(bottomly.norm)[MA.idx[i,2]] ))
## 	}
## dev.off()
## 

## ----meanSd_Bottomly_Data, echo = TRUE, fig.show='hide'------------------
pl <- meanSdPlot(exprs(bottomly.eset))

## ----meanSd_Bottomly_Data_pl---------------------------------------------
pl$gg + ylim(0,100)

## ----meanSdlog2_Bottomly_Data, echo = TRUE-------------------------------
meanSdPlot(bottomly.log2 )

## ----correlation log vs raw, fig.width=3, fig.height=3-------------------
cor(exprs(bottomly.eset)[,6],exprs(bottomly.eset)[,11])
cor(bottomly.log2[,6], bottomly.log2[,11])

## ----meanSdlog2_pseudocount_Data, echo = TRUE----------------------------
meanSdPlot(log2(exprs(bottomly.norm) + 5) )

#library(LSD)
#heatscatter(rowMeans(bottomly.log2),rowSds(bottomly.log2),
#xlab= "mean", ylab = "sd")

## ----means_and_variances_bottomly----------------------------------------
meansAndVars <- DESeq:::getBaseMeansAndPooledVariances(
  counts=exprs(bottomly.norm),sizeFactors=bottomly.sf,
   conditions=pData(bottomly.eset)$strain)
heatscatter(log(meansAndVars$baseMean),
            log(meansAndVars$baseVar), xlab= "mean of normalized counts",
            ylab = "Variance")


## ----dispersions_bottomly------------------------------------------------
dispersions <-  DESeq:::estimateAndFitDispersionsFromBaseMeansAndVariances(
  means=meansAndVars$baseMean, variances=meansAndVars$baseVar,
  sizeFactors=bottomly.sf, fitType="parametric")

dispersionsLocal <-  DESeq:::estimateAndFitDispersionsFromBaseMeansAndVariances(
  means=meansAndVars$baseMean, variances=meansAndVars$baseVar,
  sizeFactors=bottomly.sf, fitType="local")

    px = meansAndVars$baseMean
    py = dispersions$disps
    xg = 10^seq(0.9, 5, length.out = length(dispersions$disps))
    fitg = dispersions$dispFun(xg)
    fitLocal = dispersionsLocal$dispFun(xg)
    dataGG = data.frame(px, py, xg, fitg)

  (qplot(px, py, data = dataGG, ylab = "dispersion",
         xlab= "mean of normalized counts",
         log = "xy", alpha = I(1/10))
     + geom_line(aes(x = xg, y = fitg), color = "coral3")
      + geom_line(aes(x = xg, y = fitLocal), color = "dodgerblue3")
   )



## ----bottomly_PCA--------------------------------------------------------

b.PCA = prcomp(t(bottomly.log2), scale = T)
dataGG = data.frame(PC1 = b.PCA$x[,1], PC2 = b.PCA$x[,2], strain
                    = pData(bottomly.eset)$strain,
                    exp = as.factor(pData(bottomly.eset)$experiment.number))
(qplot(PC1, PC2, data = dataGG, color = exp, shape = strain) +
  scale_color_brewer( type = "qual", palette = 6))

## ----bottomly_heatmap_and_clustering, fig.height= 5----------------------
dists <- as.matrix(dist(t(bottomly.log2)))

hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(dists, trace="none", col = rev(hmcol))


## ----remove_batch_sim----------------------------------------------------
cleaned_data <- removeBatchEffect(rldSim,
                                  batch = as.character(colData(dds)$sex),
                           design =  model.matrix( ~ colData(dds)$condition))

pc_vars<- rowVars(cleaned_data)
selected_vars <- order(pc_vars, decreasing = TRUE)[seq_len(min(ntop,
                                                            length(pc_vars)))]
PCA <- prcomp(t(cleaned_data[selected_vars, ]), scale. = TRUE)
perc_var <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)


dataGG = data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                    PC3 = PCA$x[,3], PC4 = PCA$x[,4],
                    sex = colData(dds)$sex,
                    condition = colData(dds)$condition)

(qplot(PC1, PC2, data = dataGG, color =  condition, shape = sex,
       main = "PC1 vs PC2, top variable genes, after removal of sex effect", size = I(6))
+ labs(x = paste0("PC1, VarExp:", round(percentVar[1],4)),
       y = paste0("PC2, VarExp:", round(percentVar[2],4)))
+ scale_colour_brewer(palette = "Set1")
)


## ----inspect_stockori_data-----------------------------------------------
stockori <- read.csv("http://www-huber.embl.de/users/klaus/BasicR/stockori.csv")
stockori[sample(dim(stockori)[1],5), 1:10]

stockori_gtypes <- scale(as.matrix(stockori[, -seq(1,3)]))
stockori_anno <- stockori[, seq(1,3)]
stockori_anno$flow = factor(stockori_anno$flowering_binary,
                                           labels=c("short", "long"))

## ----stockori_PCA--------------------------------------------------------
stock_PCA = prcomp(stockori_gtypes, center = F, scale = F)
stock_SVD <- fast.svd(stockori_gtypes)

## ----plot_PCA_stockpori--------------------------------------------------

#table(stockori_anno$country)

  pal_func <- function(n){
  rep(brewer.pal(12, name ="Set3")[-9],(n %/% 11 + 1))[seq_len(n)]
  }
  cls = pal_func(length(levels(stockori_anno$country)))


country =  stockori_anno$country
country = ifelse(country %in% c("Cape Verdi Islands","India",
                   "Kazakhstan", "Lybia", "New Zealand",
                   "Tadjikistan", "Unknown", "USA"), "RestOfWorld", "Europe" )

country[stockori_anno$country == "Germany" ] <- "Germany"
country[stockori_anno$country == "UK" ] <- "UK"
country[stockori_anno$country == "USA" ] <- "USA"
country[stockori_anno$country == "Russia" ] <- "Russia"


ggBar <- qplot(stockori_anno$country, fill = stockori_anno$country)


dataGG <- data.frame(PC1 = stock_PCA$x[,1], PC2 = stock_PCA$x[,2],
                    PC3 = stock_PCA$x[,3], PC4 = stock_PCA$x[,4],
                    PC5 = stock_PCA$x[,5], PC6 = stock_PCA$x[,6],
                    country = stockori_anno$country,
                    flowering = stockori_anno$flow)

qplot(PC1, PC2, data = dataGG, label = country, geom = "text", color = country,
      main = "PCA for the Stockori data")


## ----PCs_from_SVD--------------------------------------------------------
PCA_from_SVD <- stockori_gtypes %*% stock_SVD$v

dataGG = data.frame(PC_SVD1 = PCA_from_SVD[,1], PC_SVD2 = PCA_from_SVD[,2],
                    PC_SVD3 = PCA_from_SVD[,3], PC_SVD4 = PCA_from_SVD[,4],
                    PC_SVD5 = PCA_from_SVD[,5], PC_SVD6 = PCA_from_SVD[,6],
                    country = country, flowering = stockori_anno$flow)

qplot(PC_SVD1, PC_SVD2, data = dataGG, label = country, geom = "text",
      color = country, main ="PCA from SVD")

## ----EIGENSTRAT_calculations---------------------------------------------
no_comp <- 10

mod <- model.matrix(~flow, data=stockori_anno)
modEig <- cbind(mod, stock_SVD$u[,seq_len(no_comp)])
stock_EIG_rotation <- lmFit(t(stockori_gtypes),
                                  modEig)$coefficients[,-c(1,2)]

stock_EIG_adj  <- stockori_gtypes - tcrossprod(stock_SVD$u[,seq_len(no_comp)],
                                               stock_EIG_rotation)
## alternative way of computation
stock_EIG_adj <- stockori_gtypes -  Reduce("+", lapply(seq_len(no_comp),
         function(k){tcrossprod(stock_SVD$u[,k], stock_EIG_rotation[,k])}))

## ----SVA_calculations----------------------------------------------------
mod0 <- model.matrix(~1, data = stockori_anno)
n_sv <- num.sv(t(stockori_gtypes), mod, method="leek")
n_sv # does not find anything
svaobj <- sva(t(stockori_gtypes), n.sv = no_comp, mod, mod0, method="irw")

modSVA <- cbind(mod,svaobj$sv)
stock_sva_rotation <- lmFit(t(stockori_gtypes),
                                  modSVA)$coefficients[,-c(1,2)]
stock_sva_adj <- stockori_gtypes - tcrossprod(svaobj$sv, stock_sva_rotation)

## ----EIGENSTRAT_and_SVA--------------------------------------------------

stock_sva_adj_PCA <-prcomp(stock_sva_adj)
stock_EIG_adj_PCA <-prcomp(stock_EIG_adj)


dataGG = data.frame(PC_adj_SVA1 = stock_sva_adj_PCA$x[,1],
                    PC_adj_SVA2 = stock_sva_adj_PCA$x[,2],
                    PC_adj_EIG1 = stock_EIG_adj_PCA$x[,1],
                    PC_adj_EIG2 = stock_EIG_adj_PCA$x[,2],
        country = country, flowering = stockori_anno$flow)
qplot(PC_adj_SVA1 , PC_adj_SVA2, data = dataGG, label = country, geom = "text",
      color = country, main ="SVA corrected data")

qplot(-PC_adj_EIG1 , -PC_adj_EIG2, data = dataGG, label = country, geom = "text",
      color = country, main ="EIGENSTRAT corrected data")

## ----naive_correction_using_only_SVD-------------------------------------
stock_naive_adj <- stockori_gtypes -
  Reduce("+", lapply(seq_len(no_comp),
         function(k){stock_SVD$d[k]*tcrossprod(stock_SVD$u[,k],
                                               stock_SVD$v[,k])}))
stock_naive_adj_PCA <-prcomp(stock_naive_adj)


dataGG = data.frame(PC_adj_naive1 = stock_naive_adj_PCA$x[,1],
                    PC_adj_naive2 = stock_naive_adj_PCA$x[,2],
        country = country, flowering = stockori_anno$flow)
qplot(PC_adj_naive1 ,PC_adj_naive2, data = dataGG, label = country, geom = "text",
      color = country, main ="naively  corrected data using only the SVD")


## ----solution:  Normalization strategies, echo = TRUE, eval=FALSE--------
## #a
## ###############################################################################
## # Pairwise MA plots for raw data
## 
## pdf("pairwiseMAsBottomlyRaw.pdf")
##  MA.idx = t(combn(1:21, 2))
## 
##   for( i in  1:dim(MA.idx)[1]){
## 	 MDPlot(log2(exprs(bottomly.eset)),
## 		c(MA.idx[i,1],MA.idx[i,2]),
## 	main = paste( sampleNames(bottomly.norm)[MA.idx[i,1]], " vs ",
## 	 sampleNames(bottomly.norm)[MA.idx[i,2]] ))
## 	}
## dev.off()
## 
## #b
## ###############################################################################
## bottomly.cS <- colSums(exprs(bottomly.eset)) * 1e-6
## multidensity( t(t(exprs(bottomly.eset)) /bottomly.cS) , xlab="mean counts", xlim=c(0, 500),
##            main = "Bottomly total sum normalized")

## ----solution:  PCA based quality control, results ='hide', fig.show='hide', fig.height=5----

b.PCA = prcomp(t(bottomly.log2), scale = T)
dataGG = data.frame(PC1 = b.PCA$x[,1], PC2 = b.PCA$x[,2], strain
                    = pData(bottomly.eset)$strain,
                    exp = as.factor(pData(bottomly.eset)$experiment.number)
                    ,lane = as.factor(pData(bottomly.eset)$lane.number))
qplot(PC1, PC2, data = dataGG, color = exp, shape = strain)


