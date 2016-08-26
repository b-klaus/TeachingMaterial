## ----style, echo=FALSE, results='asis'-----------------------------------
BiocStyle::latex()

## ----options, include=FALSE----------------------------------------------
options(digits=3, width=85, stringsAsFactors = FALSE )
opts_chunk$set(echo=TRUE,tidy=FALSE,include=TRUE,
               fig.path='DataHandling-',fig.show='hide',dev='pdf',
		 fig.width = 25, fig.height = 14, comment = '#>', dpi = 300,
		cache = TRUE, lazy.load = FALSE, background="grey93",message = FALSE )

## ----required packages and data, echo = TRUE, cache = FALSE--------------
library(TeachingDemos)
library(openxlsx)
library(multtest)
library(Biobase)
library(tidyr)
library(reshape2)
library(plyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(magrittr)
library(purrr)
library(readr)
library("DESeq2")
library("pasilla")
library("Biobase")
data("pasillaGenes")

countData <- counts(pasillaGenes)
colData <- pData(pasillaGenes)[,c("condition","type")]

pasilla_se <- DESeqDataSetFromMatrix(countData = countData,
                                      colData = colData,
                                      design = ~ condition)

library("TxDb.Dmelanogaster.UCSC.dm3.ensGene")



## ----accesRecap----------------------------------------------------------
sampleVector <- c("Alice" = 5.4, "Bob" = 3.7, "Claire" = 8.8)
sampleVector

## ----accessIndex, dependson="accesRecap"---------------------------------
sampleVector[1:2]
sampleVector[-(1:2)]

## ----accessBoolean, dependson="accesRecap"-------------------------------
sampleVector[c(TRUE, FALSE, TRUE)]
# or
subset(sampleVector, c(TRUE, FALSE, TRUE))

## ----accessBoolean2, dependson="accesRecap"------------------------------
sampleVector[sampleVector < 6]
# or
subset(sampleVector, sampleVector < 6 | names(sampleVector) == "Bob")

## ----multiDAEx-----------------------------------------------------------
pat<-read_csv("http://www-huber.embl.de/users/klaus/BasicR/Patients.csv")
pat
pat[1,c(1:3)]
pat["P1",]

## ----multiDAEx2, dependson="multiDAEx"-----------------------------------
pat$"Weight"
pat[["Weight"]][3]
# often acces via the dollar sign works without quotes (but for numbers!)
pat$Height

## ----loadBodyfat,   echo = TRUE------------------------------------------
load(url("http://www-huber.embl.de/users/klaus/BasicR/bodyfat.rda"))
dim   (bodyfat)    # how many rows and columns in the dataset?
names (bodyfat)    # names of the columns

## ----age-summary,   echo = TRUE------------------------------------------
## compute descriptive statistics for "age"
summary(bodyfat$age)
sd(bodyfat$age)
mean(bodyfat$age)
IQR(bodyfat$age)/1.349
## mean value of every variable in the bodyfat data set
apply(bodyfat, MARGIN = 2, FUN = mean)
## alternative : sapply(bodyfat, FUN = mean)

## ----age-sub, echo = TRUE------------------------------------------------
## all samples with age between 40 and 60 and height
##between 50 and 65
bodyfat[  bodyfat$age > 40 & bodyfat$age < 60 & bodyfat$height > 50 & bodyfat$height < 65, ]

### get the corresponding indices
which( bodyfat$age > 40 & bodyfat$age < 60 & bodyfat$height > 50 & bodyfat$height < 65)

### and samples
bodyfat[  which( bodyfat$age > 40 & bodyfat$age < 60 & bodyfat$height > 50 & bodyfat$height < 65)
, ]

## ----which-side-eff,  echo = TRUE----------------------------------------
max(bodyfat$age)
# 81
head(bodyfat[  !(bodyfat$age  > 81), ])

bodyfat[  !which(bodyfat$age  > 81), ]
## not equivalent!

## ----subset-example,   echo = TRUE---------------------------------------
## all samples with age between 40 and 60
## and height between 50 and 65
idx.age <- bodyfat$age > 40 & bodyfat$age < 60 & bodyfat$height > 50 & bodyfat$height < 65
subset(bodyfat, idx.age)

## only their bodyfat
subset(bodyfat, idx.age, select ="percent.fat")

## ----list_examples-------------------------------------------------------
a <- list(a = 1:3, b = "a string", c = pi, d = list(-1, -5))

a$a
a[["b"]]
a["b"]


## ----bodyfat_map, dependson = "loadBodyfat"------------------------------
head(map_dbl(bodyfat, mean))


## ----bodyfat_ex_element, dependson = "loadBodyfat"-----------------------
head(map_dbl(bodyfat, 5))

head(bodyfat[5, ])

## ----golub-read-in,   echo = TRUE, eval = TRUE---------------------------

# load the golub data
data(golub, package = "multtest")
dim(golub)
str(golub)
golub.gnames[1042,]
### expression values of CCND 3
golub[1042,]
## define group factor
gol.fac <- factor(golub.cl, levels=0:1, labels = c("ALL","AML"))

## ----sumexp, echo=FALSE, fig.show='asis'---------------------------------
par(mar=c(0,0,0,0))
plot(1,1,xlim=c(0,100),ylim=c(0,100),bty="n",
     type="n",xlab="",ylab="",xaxt="n",yaxt="n")
polygon(c(45,80,80,45),c(10,10,70,70),col=rgb(1,0,0,.5),border=NA)
polygon(c(45,80,80,45),c(68,68,70,70),col=rgb(1,0,0,.5),border=NA)
text(62.5,40,"assay(s)", cex = 3)
text(62.5,30,"e.g. 'exprs'", cex = 3)
polygon(c(20,40,40,20),c(10,10,70,70),col=rgb(0,0,1,.5),border=NA)
polygon(c(20,40,40,20),c(68,68,70,70),col=rgb(0,0,1,.5),border=NA)
text(30,40,"featureData", cex = 3)
polygon(c(45,80,80,45),c(75,75,90,90),col=rgb(.5,0,.5,.5),border=NA)
polygon(c(45,47,47,45),c(75,75,90,90),col=rgb(.5,0,.5,.5),border=NA)
text(62.5,82.5,"phenoData", cex = 3)

## ----pasilla-------------------------------------------------------------
pasilla_se

# phenoData
colData(pasilla_se)

# featureData
rowRanges(pasilla_se)

ids_pasilla <- names(rowRanges(pasilla_se))

# get the transcripts by genes
dm_genes <- transcriptsBy(TxDb.Dmelanogaster.UCSC.dm3.ensGene, by = "gene")

dm_genes

# add missing gene IDs to dm_genes
missing_ids <- setdiff(ids_pasilla, names(dm_genes))
dm_genes <- c(dm_genes, rowRanges(pasilla_se)[missing_ids])

# match the two sets
idx_se_in_dm_genes <- match(ids_pasilla, names(dm_genes))

# sanity check whether there is a missing match
any(is.na(idx_se_in_dm_genes))

# sanity check whether gene id_s match
all.equal(ids_pasilla, names(dm_genes[idx_se_in_dm_genes]))

# finally annotate the ranges
rowRanges(pasilla_se) <- dm_genes[idx_se_in_dm_genes]

# sample ranges to preview the data
rowRanges(pasilla_se)[sample(length(ids_pasilla), 5)]

# another sanity check
all.equal(ids_pasilla, names(rowRanges(pasilla_se)))


## ------------------------------------------------------------------------
## all samples with age between 40 and 60

head(filter(bodyfat, age > 40, age < 60 ))

tail(filter(bodyfat, age > 40, age < 60 ))

## ------------------------------------------------------------------------
## all samples with age of 40 or 60

head(filter(bodyfat, age == 40 | age == 60 ), 3)

tail(filter(bodyfat, age == 40 | age == 60 ), 3)

## ------------------------------------------------------------------------
## arrange by age and bodyfat
head(arrange(bodyfat, age, percent.fat),3 )

## ------------------------------------------------------------------------
## descending
head(arrange(bodyfat, desc(age), percent.fat),3 )

## ------------------------------------------------------------------------
## select fact age and heigth only
head(select(bodyfat, age, height, percent.fat ))

## select all body measures
head(select(bodyfat, weight:wrist.circum))

## exclude all body measures
head(select(bodyfat, -(weight:wrist.circum)))

## ----mutate-example------------------------------------------------------
pat <- read_csv("http://www-huber.embl.de/users/klaus/BasicR/Patients.csv")
pat$Weight[2] <- mean(pat$Weight, na.rm=TRUE)
mutate(pat, BMI <- Weight / Height^2)

## ------------------------------------------------------------------------
summarize(bodyfat, mean.age = mean(age, na.rm = TRUE),
              mean.BMI = mean( (weight*0.454) / (height*.0254)^2 ) )

## ----HTS-data,   echo = TRUE---------------------------------------------
load(url("http://www-huber.embl.de/users/klaus/BasicR/HTSdata.RData"))
head(HTSdata)

## ----split-example-------------------------------------------------------
split_HTS <- group_by(HTSdata, plate)
split_HTS

## ----split-example-2-----------------------------------------------------
split_HTS_rep  <-  group_by(HTSdata, plate, replicate)
split_HTS_rep

## ----summarize-example---------------------------------------------------
HTS_cellNumbers  <- summarize(split_HTS_rep, mean_CN = mean(CellNumber, na.rm = T))
HTS_cellNumbers

## ----custom function-example---------------------------------------------
HTS_skew <- summarize(split_HTS_rep,
HTS_CN_skew = mean(CellNumber, na.rm = T) / median(CellNumber, na.rm = T))

HTS_skew

## ----peeling example-----------------------------------------------------
plate_HTS_skew <- summarize(HTS_skew, mean_CN_skew = mean(HTS_CN_skew))

plate_HTS_skew

## ----helper functions example--------------------------------------------
## nice plotting
HTSdata <- tbl_df(HTSdata )

sample_n(HTSdata, 3)

## overview of variable
glimpse(HTSdata)

## number of different plate designs
summarize(HTSdata, n_distinct(plate) )

## number of replicates per plate for  plates
## with number greater than 15
filter(
  summarize(
      group_by(HTSdata, plate),
      rep_per_plate = n_distinct(replicate)
      )
  , plate > 15)

## ----chainingSimpleExample-----------------------------------------------

# create two vectors and calculate Euclidian distance between them
x1 <- 1:5; x2 <- 2:6

# usual way
sqrt(sum((x1-x2)^2))

# chaining method
(x1-x2)^2 %>%
sum() %>%
sqrt()


## ----chaining example----------------------------------------------------
## number of plates
summarize(HTSdata, n_distinct(plate) )

## number of replicates per plate for  plates
## plate number greater than 15

HTSdata %>%
group_by(plate)  %>%
summarize(rep_per_plate = n_distinct(replicate)) %>%
filter(plate > 15)


## ------------------------------------------------------------------------
HTSdata %>%
group_by(plate, replicate)  %>%
do( head(.,2))

## ------------------------------------------------------------------------
Preview <- HTSdata %>%
group_by(plate, replicate)  %>%
do(plate_preview = head(.,2))

Preview$plate_preview[[1]]

## ----load-protein-data---------------------------------------------------
proteins <- read_csv("http://www-huber.embl.de/users/klaus/BasicR/proteins.csv")
sample_n(proteins, 4)
proteins_pMek <- subset(proteins, Target == "pMEK")
proteins_pMek_sub <- subset(proteins_pMek, Condition == "10ng/mL HGF")

## ----plot-protein-data, fig.show='asis'----------------------------------
proteins_pMek_sub
qplot(min, Signal, data = proteins_pMek_sub, geom = "line")

## ----spread_example_tidyr,   echo = TRUE, eval = TRUE--------------------
proteins_spread <- proteins %>%
  select(-Sigma) %>%
  spread(key = min, value = Signal)

 sample_n(proteins_spread, 4)

## ----gather_example,   echo = TRUE, eval = TRUE--------------------------
proteins_gathered <- proteins_spread %>%
  gather(key = min, value = Signal, -Target, -Condition)

sample_n(proteins_gathered, 4)

## ----TSSplotPrep---------------------------------------------------------
covs <- read_csv(url("http://www-huber.embl.de/users/klaus/BasicR/DataTSS/covs.csv"),
                progress = FALSE)
names(covs) <- c("geneID", setdiff(seq(-3000, 3000), 0), "time")
sample_n(covs,10)[,1:5]

## ----TSSplot, dependson="TSSplotPrep", fig.show='asis'-------------------

data_gathered <- covs %>%
    gather(key = "pos_rel_to_tss", value = "coverage", -geneID, -time)

data_gathered$pos_rel_to_tss <- as.integer(data_gathered$pos_rel_to_tss)

sample_n(data_gathered, 10)

covs_collapsed <- data_gathered %>%
        group_by(time, pos_rel_to_tss) %>%
                  summarize(coverageC = mean(coverage))

qplot(pos_rel_to_tss, coverageC,  color = time,
      data=covs_collapsed, geom="smooth")


## ----stringHandling------------------------------------------------------
fName <- "tau138MGFP_Glu.lif - Series004 green_measure.tif"

## ----stringSplit, dependson="stringHandling"-----------------------------
f_name_split <- str_split(string=fName, pattern = "[ - _ .]")
f_name_split

## ----stringJoin, dependson="stringSplit"---------------------------------
f_name_split[[1]][c(5,2,6)]
paste0(f_name_split[[1]][c(5,2,6)], collapse = "--")

## ----readDataMichele-----------------------------------------------------
cell_imaging <- read_csv(url("http://www-huber.embl.de/users/klaus/BasicR/DataCellImaging/cellImaging.csv"))

head(cell_imaging)


## ----getDataProteomics---------------------------------------------------
data_dir <- file.path("DataProteomics")
list.files(data_dir)

## ----handleP1, dependson="getDataProteomics"-----------------------------

# get the proteomics data files
prot_data <- file.path(data_dir, list.files(data_dir, pattern="ratio"))
prot_data

namesExp <- str_sub(sapply(str_split(prot_data, "_"), "[", 1), 16)
namesExp



## ----importP, dependson="handleP1"---------------------------------------

prot_data  <- map_df(prot_data, read_csv, .id = "experiment"
, col_names = c("S1", "S2", "S3"))

prot_data <- mutate(prot_data, experiment = mapvalues(experiment,
                   from = unique(experiment), to = namesExp))



## ----importMetadata, dependson="joinP"-----------------------------------

sample_metadata <- read_csv(file.path(data_dir, "SampleLegend.csv"),
                           skip=1)

names(sample_metadata)[c(1,3:5)] <- c("sample_name", "S1", "S2", "S3")
sample_n(sample_metadata, 10)

# load protein identifiers
protein_list <-  read_csv(file.path(data_dir, "ProteinList.csv"),
                         skip=0)
sample_n(protein_list, 10)

## ----joinMeta, dependson="importMetadata"--------------------------------

## add protein names as rownames
prot_data$Accession <- rep(protein_list$Accession, length(namesExp))

# modify sample metadata to have one line per sample
sample_metadata <- gather(sample_metadata, key = "sample", value= "type", S1:S3)
sample_metadata <-  mutate(sample_metadata,
                                sample_id = paste0(sample_name, "_", sample))
sample_metadata <- arrange(sample_metadata, sample_id)



sample_metadata

## ----showProt, dependson="joinMeta"--------------------------------------
sample_n(prot_data, 10)

## ----sol-bodyfat,   echo = TRUE, results='hide'--------------------------
#a
########################################################################
sapply(bodyfat, mean)
sapply(bodyfat, sd)

# or

map_dbl(bodyfat, mean)
map_dbl(bodyfat, sd)

#b
########################################################################
small.idx <-  bodyfat$height < mean(bodyfat$height) - sd(bodyfat$height)
subset(bodyfat, small.idx )
## or
filter(bodyfat, small.idx)

#c
########################################################################
light.idx <-  bodyfat$weight < mean(bodyfat$weight) - sd(bodyfat$weight)
subset(bodyfat, light.idx)
## or
filter(bodyfat, light.idx)

#d small and light people
########################################################################
small.and.light.idx  <- small.idx & light.idx
subset(bodyfat, small.and.light.idx)
# or
filter(bodyfat, small.idx, light.idx)

## ----HTSscreen-solution,   echo = TRUE, results='hide'-------------------
#a
########################################################################
load(url("http://www-huber.embl.de/users/klaus/BasicR/HTSdata.RData"))

#b
########################################################################
HTS_n_cs <-  HTSdata %>%
group_by(plate, replicate) %>%
summarize(mean_CN = mean(CellNumber, na.rm = T),
          sd_CN = sd(CellNumber, na.rm = T))
HTS_n_cs

#c
########################################################################
TransR <- HTSdata

TransR$TransR <- ifelse(is.finite(TransR$TransR), TransR$TransR ,0)

TransRlist <-  TransR %>%
group_by(plate, replicate)  %>%
dplyr::select(TransR)  %>%
do(TRlist = identity(.))

TransRlist$TRlist[[1]]


## ----imagingDataSol,  dependson="readDataMichele",  echo = TRUE, results='hide'----
#a

label_Info <- str_split(string=cell_imaging$Label, pattern = "[ - _ .]")
extractedInfo <- map_df(label_Info,
                        function(x){data.frame(series = x[5],
                                               medium = x[2],
                                               label = x[6])})

head(extractedInfo)

#b

## add columns to the data thaht give the categories

cell_imaging <- mutate(cell_imaging, green_red=ifelse(is.na(str_match(Label, "green")[,1]),
                                                   "Red", "Green"),
              gal_glu=ifelse(is.na(str_match(Label, "Gal")[,1]), "Glu", "Gal"))
## check variable types
glimpse(cell_imaging)

#c

## group by category and get the mean
cell_imaging %>%
  group_by(green_red, gal_glu) %>%
  summarize(mean(Mean))



## ----inspect-golub,   echo = TRUE, results='hide'------------------------
#a
########################################################################
golub[1042,gol.fac=="AML"]
#b
########################################################################
meanALL <- apply(golub[,gol.fac=="ALL"], 1, mean)
meanAML <- apply(golub[,gol.fac=="AML"], 1, mean)

#c
########################################################################
head(meanALL[order(meanALL, decreasing = TRUE)], 100)
head(golub[order(meanALL, decreasing = TRUE),  ], 20)


(golub.gnames[order(meanALL, decreasing = TRUE),2])[1:3]

## ----handle-expression-sets,   echo = TRUE, results='hide'---------------
#a
########################################################################
library(Biobase)
data(sample.ExpressionSet)
sample.ExpressionSet
exprs(sample.ExpressionSet)
slotNames(sample.ExpressionSet)

#b
########################################################################
varLabels(sample.ExpressionSet)
sampleNames(sample.ExpressionSet)
varMetadata(sample.ExpressionSet)

#c
########################################################################
annotation(sample.ExpressionSet)
featureNames(sample.ExpressionSet)

#d
########################################################################
controls <- grep("AFFX", featureNames(sample.ExpressionSet), value = TRUE)
str(controls)

#e
########################################################################
phenoData(sample.ExpressionSet)
head(pData(sample.ExpressionSet))

## ----sessionInfo, cache=FALSE, results='asis'----------------------------
toLatex(sessionInfo())

