#install.packages("miscTools")
library(dplyr)
library(miscTools)
setwd("E:/Bioninformatics/year4,sem1/Biostatistics/Assignments/Assignment 3")
cancerDataProteome <- read.csv("77_cancer_proteomes_CPTAC_itraq.csv")
cancerDataClinical <- read.csv("clinical_data_breast_cancer.csv")
PAMProteins <- read.csv("PAM50_proteins.csv")

#------------------------------------1--------------------------------------

names(cancerDataProteome)[names(cancerDataProteome)=="RefSeq_accession_number"] <-"RefSeqProteinID"
cancerDataProteome <-na.omit(cancerDataProteome)
newClinical <- cancerDataClinical[,c("Complete.TCGA.ID","HER2.Final.Status")]

proteomeUpdated <- merge(x=cancerDataProteome,y=PAMProteins,by="RefSeqProteinID")
proteomeNames <- colnames(proteomeUpdated)
proteomeNamesUpdated <- proteomeNames

for (i in c(1:length(proteomeNames))){
  if(i>=4 && i<=length(proteomeNames)-3){
    tempSeq <- substr(proteomeNames[i],11,14)
    tempi <- substr(proteomeNames[i],1,2)
    tempStr <- substr(proteomeNames[i],4,7)
    fullName <- paste(tempSeq,tempi,tempStr,sep = "-")
    proteomeNamesUpdated[i] = fullName
  }
  else if(i<=3 || i>86){
    proteomeNamesUpdated[i]= proteomeNames[i]
  }
  else{
    if(i==84){
      tempSeq <- substr(proteomeNames[i],10,15)
      tempi <- substr(proteomeNames[i],9,9)
      tempStr <- substr(proteomeNames[i],1,7)
      fullName <- paste(tempSeq,tempi,tempStr,sep = "-")
      proteomeNamesUpdated[i] = fullName
    }
    else{
      tempSeq <- substr(proteomeNames[i],10,14)
      tempi <- substr(proteomeNames[i],8,8)
      tempStr <- substr(proteomeNames[i],1,6)
      fullName <- paste(tempSeq,tempi,tempStr,sep = "-")
      proteomeNamesUpdated[i] = fullName
    }
  }
}

names(proteomeUpdated) <- proteomeNamesUpdated
cancerDataProteomeUpdated <- proteomeUpdated
# finalData <- data.frame(data <- data.frame(matrix(NA, nrow = nrow(cancerDataProteomeUpdated),)))

refSeqProteinID <- proteomeUpdated$RefSeqProteinID
proteomeUpdatedT <- t(proteomeUpdated[,4:86])
tempName <- as.matrix(rownames(proteomeUpdatedT))
rownames(proteomeUpdatedT) <-  NULL
proteomeUpdatedT = as.data.frame(proteomeUpdatedT)
names(proteomeUpdatedT) <- refSeqProteinID
proteomeUpdatedT <- as.matrix(proteomeUpdatedT)
proteomeUpdatedT <- insertCol(proteomeUpdatedT,24,tempName,cName = "Complete.TCGA.ID")
final <- merge(x=newClinical, y=proteomeUpdatedT, by="Complete.TCGA.ID") 

#-----------------------------------2---------------------------------------

finalM <- final[,3:25]
herEn <- as.numeric(factor(final[,"HER2.Final.Status"], exclude = NULL))
finalDataEn <- insertCol(as.matrix(finalM),1,herEn,cName = "HER2.Final_Encoded")

finalDataEn <- as.matrix(finalDataEn)
corrRes <- c()

for (i in c(2:ncol(finalDataEn)) ){
  res <- cor(as.numeric(finalDataEn[,i]), y = as.numeric(finalDataEn[,1]), use = "everything", method = "pearson")
  corrRes <-c(corrRes,res)
}

corrResPostitive <- abs(corrRes)
orderedCorrRes <- order(-corrResPostitive)
finalDataOrdered <- as.data.frame(finalM[,orderedCorrRes])

# corrMean <- mean(corrResPostitive)
res <- sum(ifelse(corrResPostitive[orderedCorrRes]>0.05,1,0))
finalDataOrdered <- finalDataOrdered[,1:res]

colnames(finalDataOrdered)
ncol(finalDataOrdered)

#------------------------------3--------------------------------------------

final <- as.data.frame(final)
nrow(final)
positiveHer <- as.data.frame(final[final$HER2.Final.Status == "Positive", ])
negativeHer <- as.data.frame(final[final$HER2.Final.Status == "Negative", ])
posNeg <- rbind(positiveHer,negativeHer)
tRes = c()
pRes = c()

for (i in c(3:ncol(posNeg))) {
  tRes <- t.test(as.numeric(positiveHer[,i]), as.numeric(negativeHer[,i]), data=final)
  # statValue <- t.test(as.numeric(posNeg[,i]) ~ posNeg$HER2.Final.Status, data=final)$statistic
  # pValue <- t.test(as.numeric(posNeg[,i]) ~ posNeg$HER2.Final.Status, data=final)$p.value
  print(tRes)
  # tRes <- c(tRes, statValue)
  # pRes <- c(pRes, pValue)
}


posNegP <- posNeg[,3:ncol(posNeg)]
cp <- sum(ifelse(pRes<0.05,1,0))
cpFeature <- ifelse(pRes<0.05,names(posNegP),NA)
cpFeature <- as.matrix(cpFeature)
cpFeature <- na.omit(cpFeature)

cat("It is significantly different for some columns,",cp,"columns, the columns are: \n",
    cpFeature)

tResPostitive <- abs(tRes)
orderedTRes <- order(-tResPostitive)
tOrdered <- as.data.frame(finalM[,orderedTRes])
res <- sum(ifelse(tResPostitive[orderedTRes]>0.05,1,0))
tOrdered <- tOrdered[,1:res]

colnames(tOrdered)
ncol(tOrdered)

cat("Through using t-test with threshold 0.05, we have found that; \n 
    Only 21 columns got selected which are: NP_004439, NP_001159403, NP_008950, \n
    NP_058518, NP_001116539, NP_004487, NP_000116, NP_000917, NP_054895, NP_077006, \n
    NP_001035932, NP_003003, NP_000415, NP_000517, NP_065178, NP_005931, NP_005219\n
    NP_006614, NP_002408, NP_001153651, NP_006836")

cat("Through using correlation with threshold 0.05, we have found that; \n 
    Only 18 columns got selected which are: NP_004439, NP_001159403, NP_008950, \n
    NP_001116539, NP_058518, NP_000917, NP_000116, NP_065178, NP_077006, NP_001035932, \n
    NP_000415, NP_000517, NP_004487, NP_003003, NP_054895, NP_005219, NP_005931, NP_006614")

commonNames <- intersect(names(tOrdered), names(finalDataOrdered))

cat("Common features between both ways are: \n", commonNames[1:9], "\n", commonNames[9:18])


# library(waldo)
# difference <- compare(tOrdered, finalDataOrdered)
