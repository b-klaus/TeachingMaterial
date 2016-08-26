## ----style, echo=FALSE, results='asis'-----------------------------------
BiocStyle::latex()

## ----options, include=FALSE----------------------------------------------
options(digits=3, width=80)
opts_chunk$set(echo=TRUE,tidy=FALSE,include=TRUE,
               dev='pdf', fig.width = 6, fig.height = 3.5, comment = '  ', dpi = 300,
		cache = T, lazy.load = FALSE, background="grey93" )

## ----required packages and data, echo = TRUE, message = FALSE------------

set.seed(999)

library(gplots)
library(RColorBrewer)
library(ggplot2)
library(plyr)
library(dplyr)
library(magrittr)
library(tidyr)
library(mutoss)
library(qvalue)
library(st)
library(ALL)
library(hgu95av2.db)
library(genefilter)
library(IHW)
library("DESeq")
library("ggplot2")
library("methods")
library("airway")
library("DESeq2")
data("airway")


## ----plantData, cache=FALSE----------------------------------------------

fBar <- function(x){
    tmp <- mean(x, na.rm = T)
    data.frame(ymax = tmp, y = tmp, ymin = tmp)
}

variety1 <- rnorm(10, mean = 5)
variety2 <- rnorm(10, mean = 7, sd = 1/sqrt(2))

plantData <- data.frame(height  = c(variety1,variety2) ,
                     variety = sort(rep(c("variety1", "variety2"),10)) )

head(plantData)



(ggplot(aes(x = variety, y = height, color = variety, shape = variety, size = 3),
        data = plantData)

    + geom_jitter(width = 0.3, height = 0)


    + scale_color_brewer(palette = "Set1")

    + ylab("height in cm")

    + geom_errorbar(stat = "summary", width = 0.4,
                    fun.data = fBar, size = 0.4))




## ----permTest, dependson="plantsExample"---------------------------------

## helper function to compute permutation p-value

mDiff <- function(data, group){
  data$group <- NULL
  data$group <- group
  tmp <-  data %>%
    group_by(group)  %>%
    summarize(m = mean(data, na.rm=TRUE))
  as.numeric(tmp[2,"m"] - tmp[1,"m"])
}


## function to compute the permutation test
permTestTwoGroups <- function(group1, group2,
                              twoSided = TRUE, permutations = 1e4){

  stopifnot(is.numeric(group1),  is.numeric(group2),
             length(group1) > 0,  length(group2) > 0,
             is.vector(group1), is.vector(group2),
            is.logical(twoSided))

  inputData <- data.frame(data  = c(group1, group2),
                     group = rep(c("group1", "group2"),
                                   c(length(group1), length(group2))))

  # compute the observed difference between the groups
  obsDiff  <-  mDiff(inputData, inputData$group)

  # compute sampling distribution and p--value
  samplingDist <- c(replicate(as.integer(permutations),
                              mDiff(inputData, sample(inputData$group))),
                              obsDiff)

  # compute two-sided p-value
  pvalP <- 2*min(1 - ecdf(samplingDist)(abs(obsDiff)),
                 ecdf(samplingDist)(abs(obsDiff)))

  # create plot of the sampling distribution
  samplingDistPlot <- (qplot(samplingDist, fill = I("orange4"),
                             main = "Sampling distribution mean difference",
                              binwidth = 0.1)
  +  geom_vline(xintercept = obsDiff, size = 2))

  return(list(
    samplingDist = samplingDist,
    samplingDistPlot = samplingDistPlot,
    obsDiff = obsDiff,
    pval = pvalP
  ))

}


## ----getPvalue, dependson="permTest", cache = TRUE-----------------------
## draw 1e4 times a permutation of the  sample labels  and compute the
## two sided p--value


testResult <- permTestTwoGroups(group1 = variety1,
                  group2 = variety2)

testResult$pval

testResult$samplingDistPlot


## ----setupALL------------------------------------------------------------

data("ALL")
bALL <- ALL[, substr(ALL$BT,1,1) == "B"]
fusALL <- bALL[, bALL$mol.biol %in% c("BCR/ABL", "NEG")]
fusALL$mol.biol <- factor(fusALL$mol.biol)
fusALL
sample_n(pData(fusALL), 10)

groupsALL <- fusALL$mol.biol
expALL <- exprs(fusALL)

anno_fusALL <- plyr::ddply(AnnotationDbi::select(hgu95av2.db,
                                  keys=rownames(expALL),
                                  columns = c("SYMBOL", "GENENAME", "ENSEMBL"),
                                  keytype="PROBEID"), "PROBEID", function(X){X[1,]})




## ----showALL, dependson="setupALL"---------------------------------------

head(groupsALL)
head(expALL[, 1:5])
head(anno_fusALL)



## ----getBCL2-------------------------------------------------------------
anno_fusALL[1152,]

## ----BCL2-density-cdf, echo = TRUE---------------------------------------

f <-function(x){dnorm(x, 8.6, 0.5)}
F <-function(x){pnorm(x, 8.6, 0.5)}

x <- seq(6,11,0.01)
dataGG <- data.frame(x = x, y = f(x))
dataGG <- mutate(dataGG,  area = ifelse(x < 8, "in", "out"  ))

p<-qplot(data = dataGG, x = x, y = y, geom="line")
p<-p + geom_area(aes(ymax = y, fill = area)) +  guides(fill=FALSE)
p<-p + xlab("Gene-Expression") + ylab("density f") + annotate("text", x = 8, y = 0, label = "8")
p + labs(title = "P(X<=8)= 0.115") + scale_fill_brewer(palette = "Dark2")


dataGG = data.frame(x = x, y = F(x))
dataGG <- mutate(dataGG,  area = ifelse(x < 8, "in", "out"  ))

p<-qplot(data = dataGG, x = x, y = y, geom="line")
p<- p + annotate("text", x = 8, y = 0, label = "8")
p<- p	+ annotate("text",y = .16, x = 6, label = "0.115")

p <- p + geom_segment(y=0, yend = F(8), x=8, xend = 8, size = I(1))
p <- p + geom_segment(y=F(8), yend = F(8), x=6, xend = 8, size = I(1))
p + labs(title = "P(X<=8)= 0.115") + xlab("Gene-Expression") + ylab("distribution F")



## ----BCL2-1, echo = TRUE-------------------------------------------------
1 - pnorm(9, 8.6, 0.5)

## ----BCL2-2, echo = TRUE-------------------------------------------------
pnorm(9, 8.6, 0.5) - pnorm(8, 8.6, 0.5)

## ----BCL2-3, echo = TRUE-------------------------------------------------
qnorm(0.025,8.6,0.5)

## ----BCL2-4, echo = TRUE-------------------------------------------------
x <- rnorm(1000,8.6,0.5)

## ----BCL2-5, echo = TRUE-------------------------------------------------
mean(x)
#and
sd(x)

## ----tTestBCL2,  echo = TRUE---------------------------------------------
t.test(expALL[1152,] ~ groupsALL)
### alternative call
t.test(expALL[1152, groupsALL == "BCR/ABL"],
       expALL[1152, groupsALL == "NEG"] )

## ----wilcoxTestBCL2,  echo = TRUE----------------------------------------
wilcox.test(expALL[1152,] ~ groupsALL)

## ----permTestBCL2--------------------------------------------------------
pTest <- permTestTwoGroups(expALL[1152, groupsALL == "BCR/ABL"],
       expALL[1152, groupsALL == "NEG"], permutations = 1e4 )

pTest$pval

pTest$samplingDistPlot

## ----WilcoxSim-----------------------------------------------------------


x <- rnorm(10)
qqnorm(x)
qqline(x)


y <- rnorm(10)
wilcox.test(x,y)



wc <- function(){
  x <- rnorm(10)
  y <- rnorm(10, sd = 15)
  tt <- wilcox.test(x,y)
  tt$p.value
}

pVals_WC <- replicate(1000, expr = wc())

hist(pVals_WC, col = "tan3", main = "Wilcoxon P--values")
prop.table(table(pVals_WC < 0.05))


set.seed(999)

ttest <- function(){
  x <- rnorm(10)
  y <- rnorm(10, sd = 15)
  tt <- t.test(x,y)
  tt$p.value
}

pValsT <- replicate(1000, expr = ttest())

hist(pValsT, col = "coral3", main = "t-test P--values")
prop.table(table(pValsT < 0.05))


## ----Gentetic association test-------------------------------------------
disease=c(rep("no",180),rep("yes",20),rep("no",40),rep("yes",10))
genotype=c(rep("AA",200),rep("aa",50))
tab=table(genotype,disease)
tab

## ----Gentetic association test 2-----------------------------------------
p <- mean(disease == "yes")
p

## ----Gentetic association test 3-----------------------------------------
rbind(c(1-p,p)*sum(genotype=="aa"),c(1-p,p)*sum(genotype=="AA"))

## ----Gentetic association test 4-----------------------------------------

chisq.test(tab)$p.value

## ----Gentetic association test 5-----------------------------------------
tab=tab*10
chisq.test(tab)$p.value

## ----chi2-enrich---------------------------------------------------------
dat1 <- matrix(c(100,200,3000,6000),2,byrow=TRUE)
## Chi2 test
chisq.test(dat1)

## ----chi2-enrich-fisher--------------------------------------------------
dat1 <- matrix(c(100,200,3000,6000),2,byrow=TRUE)

## Fisher test
fisher.test(dat1)


## ----chi2-enrich-prop----------------------------------------------------
dat1 <- matrix(c(100,200,3000,6000),2,byrow=TRUE)
## Comparison of proportions of oncogenes in the two subsets
dat1.prop <- matrix(c(dat1[1,1], dat1[2,1],dat1[1,1] + dat1[1,2],
               dat1[2,1] + dat1[2,2]), 2,2,byrow=TRUE)
prop.test(dat1.prop[1,] ,dat1.prop[2,])


## ----chi2-enrich-2-------------------------------------------------------
dat2 <- matrix(c(50,250,3000,6000),2,byrow=TRUE)
## Chi2 test
chisq.test(dat2)
## Fisher test
fisher.test(dat2)
## Comparison of proportions of oncogenes in the two subsets
dat2.prop <- matrix(c(dat2[1,1], dat2[2,1],dat2[1,1] + dat2[1,2],
               dat2[2,1] + dat2[2,2]), 2,2,byrow=TRUE)
prop.test(dat2.prop[1,] ,dat2.prop[2,])

## ----simulate z scores---------------------------------------------------

sd.true =1
eta0.true = 0.75

get.random.zscore = function(m=200)
{

  m0 = m*eta0.true
  m1 = m-m0
  z = c(  rnorm(m0, mean=0, sd=sd.true),
       rnorm(m1, mean=2, sd=1))
  #z = sign(rnorm(length(z)))*z

  return(z)
}

set.seed(555)
z <- get.random.zscore(200)
pv = 1- pnorm(z, sd = sd.true)


## ----SchwSpot_plot-------------------------------------------------------
(ggplot2::qplot(sort(1-pv), 1:200, xlab = "1 - p-values", ylab = "No. of hypotheses",
       main = "Schweder and Spotvoll plot")
 + geom_abline(intercept = 0, slope = 200*eta0.true, aes(color = "coral3")))

## ----pvalue_histogram----------------------------------------------------
ggplot2::qplot(x = pv, xlab = "p-values", main = "Histogram of p-values, correct null
        distribution",
       fill = I("navyblue"))

## ----pvalue_histogram_wrong_null-----------------------------------------
ggplot2::qplot(x = pnorm(z, sd = 2) , xlab = "p-values", main = "Histogram of p-values,
      variance of null distribution too high",
       fill = I("coral3"))

## ----pvalue_histogram_wrong_null_2---------------------------------------
ggplot2::qplot(x = pnorm(z, sd = 0.5) , xlab = "p-values", main = "Histogram of p-values,
      variance of null distribution too low",
       fill = I("chartreuse4"))

## ----p-value_adjustments-------------------------------------------------
alpha = 0.05
pv.BH <- p.adjust(pv, method = "BH")
table(pv.BH < 0.05)
(ggplot2::qplot(rank(pv), pv, xlab = "p-value rank", ylab="p-values"
       ,main = "Visualization of the BH procedure")
+ geom_abline(intercept = 0, slope = alpha/200, aes(color = "coral3"))
+ ylim(c(0, 0.2) ))
pv.FWER <- p.adjust(pv, method = "bonferroni")
table(pv.FWER < 0.05)

## ----mutoss oracleBH-----------------------------------------------------
 ABH_pi0_est(pv)
pv.OracleBH <- oracleBH(pValue=pv, alpha=alpha, pi0=0.75)

## ----mutoss Qvalue-------------------------------------------------------
pv.Qvals <- Qvalue(pv)
table(pv.Qvals$qValues < 0.05)

## ----modTestBCL2,  echo = TRUE-------------------------------------------

modt.stat(t(expALL), groupsALL)[1152]
### p-value using the normal distribution
2 - 2*pnorm(modt.stat(t(expALL), groupsALL)[1152])
### slightly higher t-value than the ordinary t-test
t.test(expALL[1152,] ~ groupsALL, var.equal=FALSE)$statistic

## ----shrinkage_t_ALL-----------------------------------------------------
sts <- shrinkt.stat(t(expALL), groupsALL)
### p-value using the normal distribution
pval.st <- 2 - 2*pnorm(abs(sts))
hist(pval.st, col = "lavender", main = "Histogram of p-values for
     shrinkage t statistics of the ALL data")

## ----shrinkage t multiple testing----------------------------------------
pv.st.BH <- p.adjust(pval.st, method = "BH")
table(pv.st.BH< 0.05)


## ----libraries,results='hide'--------------------------------------------
data("pasillaGenes")

## ----DESeq2,cache=TRUE,results='hide', warning = FALSE-------------------
cds  <- estimateSizeFactors( pasillaGenes )
cds  <- estimateDispersions( cds )
fit1 <- fitNbinomGLMs( cds, count ~ type + condition )
fit0 <- fitNbinomGLMs( cds, count ~ type  )

## ----DESeq3,cache=TRUE---------------------------------------------------
res <- data.frame(
filterstat <- rowMeans(counts(cds)),
pvalue    <- nbinomGLMTest( fit1, fit0 ),
row.names <-featureNames(cds) )

## ----headres-------------------------------------------------------------
dim(res)
head(res)

## ----pass,echo=FALSE,cache=FALSE-----------------------------------------
theta = 0.4
pass = with(res, filterstat > quantile(filterstat, theta))

## ----figscatterindepfilt, fig.show = 'hide'------------------------------
with(res,
  plot(rank(filterstat)/length(filterstat), -log10(pvalue), pch=16, cex=0.45))

## ----figecdffilt, fig.show = 'hide'--------------------------------------
trsf = function(n) log10(n+1)
plot(ecdf(trsf(res$filterstat)), xlab=body(trsf), main="")

## ----badfilter1,cache=TRUE, fig.show = 'hide'----------------------------
badfilter = as.numeric(gsub("[+]*FBgn", "", rownames(res)))

## ----badfilter2,echo=FALSE-----------------------------------------------
stopifnot(!any(is.na(badfilter)))

## ----figbadfilter, fig.show = 'hide'-------------------------------------
plot(rank(badfilter)/length(badfilter), -log10(res$pvalue), pch=16, cex=0.45)

## ----pBH1,cache=TRUE-----------------------------------------------------
theta = seq(from=0, to=0.5, by=0.1)
pBH = filtered_p(filter=res$filterstat, test=res$pvalue, theta=theta, method="BH")

## ----pBH2----------------------------------------------------------------
head(pBH)

## ----figrejection,fig.width=5.5,fig.height=5.5, fig.show = 'hide'--------
rejection_plot(pBH, at="sample",
               xlim=c(0, 0.5), ylim=c(0, 2000),
               xlab="FDR cutoff (Benjamini & Hochberg adjusted p-value)", main="")

## ----filtered_R1,cache=TRUE----------------------------------------------
theta = seq(from=0, to=0.8, by=0.02)
rejBH = filtered_R(alpha=0.1, filter=res$filterstat, test=res$pvalue, theta=theta, method="BH")

## ----fignumreject,fig.width=5.5,fig.height=5.5, fig.show = 'hide'--------
plot(theta, rejBH, type="l",
     xlab=expression(theta), ylab="number of rejections")

## ----differentstats,cache=TRUE-------------------------------------------
filterChoices = data.frame(
  `mean`   = res$filterstat,
  `geneID` = badfilter,
  `min`    = rowMin(counts(cds)),
  `max`    = rowMax(counts(cds)),
  `sd`     = rowSds(counts(cds))
)
rejChoices = sapply(filterChoices, function(f)
  filtered_R(alpha=0.1, filter=f, test=res$pvalue, theta=theta, method="BH"))

## ----colours,results='hide'----------------------------------------------
myColours = brewer.pal(ncol(filterChoices), "Set1")

## ----figdifferentstats,fig.width=5.5,fig.height=5.5, fig.show = 'hide'----
matplot(theta, rejChoices, type="l", lty=1, col=myColours, lwd=2,
        xlab=expression(theta), ylab="number of rejections")
legend("bottomleft", legend=colnames(filterChoices), fill=myColours)

## ----Deseq2, message=FALSE, warning=FALSE--------------------------------

dds <- DESeqDataSet(se = airway, design = ~ cell + dex)
dds <- DESeq(dds)
de_res <- as.data.frame(results(dds))

## ------------------------------------------------------------------------
colnames(de_res)

## ----message=FALSE, warning=FALSE----------------------------------------
ihw_res <- ihw(pvalue ~ baseMean,  data=de_res, alpha = 0.1)

## ------------------------------------------------------------------------
rejections(ihw_res)

## ------------------------------------------------------------------------
head(adj_pvalues(ihw_res))
sum(adj_pvalues(ihw_res) <= 0.1, na.rm = TRUE) == rejections(ihw_res)

## ------------------------------------------------------------------------
padj_bh <- p.adjust(de_res$pvalue, method = "BH")
sum(padj_bh <= 0.1, na.rm = TRUE)

## ------------------------------------------------------------------------
head(weights(ihw_res))

## ------------------------------------------------------------------------
weights(ihw_res, levels_only=TRUE)

## ----plot_ihw_res, fig.width = 5, fig.height = 3-------------------------
plot(ihw_res)

## ------------------------------------------------------------------------
ihw_res_df <- as.data.frame(ihw_res)
colnames(ihw_res_df)

## ----geneExp-Exercise, echo = TRUE, results ='hide'----------------------
# a
#######################################################
1 - pnorm(1.2, 1.6, 0.42)

# b
#######################################################
pnorm(2, 1.6, 0.42) - pnorm(1.2, 1.6, 0.42)

# c
#######################################################
pnorm(2.4, 1.6, 0.42)  - pnorm(0.8, 1.6, 0.42)


# d
#######################################################
qnorm(0.025, 1.6, 0.42)
qnorm(0.975, 1.6, 0.42)

# e
#######################################################
test.sample = rnorm(1000, 1.6, 0.42)
mean(test.sample)
sd(test.sample)

## ----sol-fever,   echo = TRUE, results='hide'----------------------------
deaths <-c(210,122)
tot.cases <-c(747,661)

prop.test(deaths, tot.cases)

chisq.test(rbind(deaths, tot.cases-deaths))
fisher.test(rbind(deaths, tot.cases-deaths))

## ----sol-t-test,   echo = TRUE, results='hide'---------------------------

t.test(expALL[8197,] ~ groupsALL, var.equal=FALSE)

### strong difference between groups ...
### confirmed by wilcoxon test
wilcox.test(expALL[8197,] ~ groupsALL)

### the moderated t-test is also significant
modt.stat(t(expALL), groupsALL)[8197]
### p-value using the normal distribution
2 - 2*pnorm(abs(modt.stat(t(expALL), groupsALL)[8197]))

## ----sol gol multiple testing, results='hide', fig.keep='none'-----------
golub.qval <- Qvalue(pval.st)
significant <- as.factor(golub.qval$qValues < 0.05)
levels(significant) <- ifelse(levels(significant) , "yes", "no")

table(significant)

ggplot2::qplot(pval.st , xlab = "p-values of shrinkage t statistics",
               main = "Histogram of p-values,
       golub data, q-value < 0.05",
       fill = significant)

significant <- as.factor(pv.st.BH < 0.05)
levels(significant) <- ifelse(levels(significant) , "yes", "no")

(ggplot2::qplot(pval.st , xlab = "p-values of shrinkage t statistics",
                main = "Histogram of p-values,
       golub data, BH adjusted p-value < 0.05",
       fill = significant)  + scale_fill_brewer(type = "qual", palette = 8))

