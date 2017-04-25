source("estimation_helper.R");
library("plyr");
TESTING = 0;
##########
##HELPER##
##########
readTrainingData = function(strPath) {
  colClass = c("integer","character","character","numeric","numeric","numeric","numeric","character");
  data = read.csv(strPath,sep = "\t",header=TRUE,colClasses = colClass);
  data$Family = paste(data$Order,data$Family,sep="_")
  data$catg = rep("",dim(data)[1]);
  return(data);
}
#Format of area data: Sample\t Area so need to subset the relevent sample and to extract the areas
readAreaData = function(strPath,sampleNumber) {
  data = read.csv(strPath,sep = "\t",header=TRUE);
  data = data[which(data$Sample == sampleNumber),];
  return(data$Area/ sum(data$Area));
}
#Foramt Sample\tOrder\tFamily\tGenus\tSpecies\tReads
readSequencingData = function(strSeqPath,sampleNumber) {
  colClass = c("character","character","character","character","character","numeric");
  data = read.csv(strSeqPath,sep = "\t",header=TRUE,colClasses = colClass);
  data$Family = paste(data$Order,data$Family,sep="_")
  seqData = data[which(data$Sample == sampleNumber),];
  seqData$catg = rep("", dim(seqData)[1])
  seqData$fullName = paste(seqData$Family,seqData$Genus,seqData$Species);
  return(seqData);
}
countFamiliesInTrainingData = function(dfTrainingData){
  dfData = dfTrainingData;
  dfData$Cat = dfData$Family;
  families_catg = unique(dfData$Cat);
  families_occ = matrix(rep(0,length(families_catg)),ncol=length(families_catg));
  colnames(families_occ) = families_catg;
  samples = unique(dfData$Sample);
  for (i in 1:length(samples)){
    sample = samples[i];
    dfData.sample = dfData[which(dfData$Sample==sample),];
    families.in.sample = unique(dfData.sample$Cat);
    for(f in families.in.sample) {
      families_occ[grep(f,families_catg)] = families_occ[grep(f,families_catg)] + 1;
    }
  }
  return(families_occ)
}
updateCatg = function(dfData,catg.family) {
  for(i in 1:(dim(dfData)[1])) {
    if(dfData[i,]$Family %in% catg.family) {
      dfData[i,]$catg = dfData[i,]$Family;
    }
    else {
      dfData[i,]$catg = dfData[i,]$Order;}
  }
  return(dfData)
}

########
##DEFS##
########
trainigDataPath = "trainingSet.txt";
strAreaDataPath = "allBatch.areas.txt";
strSeqPath = "allBatch.seq.31.3.15";
EM.trails = 100000;
isL2 = 0;
########
##MAIN##
########

#Parse arguments
args <- commandArgs(trailingOnly = TRUE);
dirPath = args[1];
suffixName = args[2];
sampleNumber = args[3];
initalSeed = as.integer(args[4]);
seedsNumber = as.integer(args[5]);

#Read files
trainingData = readTrainingData(trainigDataPath);
areas = readAreaData(strAreaDataPath,sampleNumber);
sequncingData = readSequencingData (strSeqPath,sampleNumber);
#Start the learing process - choose the families with sufficent #occurences and learn coefficents 
familtCountTH = 4; #Should be fixed,minimal #Occurences for a family in order to be learned.
familyCount = countFamiliesInTrainingData(trainingData);
catg.family = colnames(familyCount)[which(familyCount>familtCountTH)];
#Update Catg fildes in training and sequencing data
trainingData = updateCatg(trainingData,catg.family);
sequncingData = updateCatg(sequncingData,catg.family);
all.catg = unique(c(trainingData$catg,sequncingData$catg));
#Learn data
lData = convertRawDataToTrainingData(trainingData,all.catg)
sampels.areas = lData[[1]];
sampels.reads = lData[[2]];
W = learnW(sampels.reads,sampels.areas);
#Predict
if(isL2) {
  sequncingData$normalizedReads = predictReads(sequncingData,W,all.catg);
}else {
  sequncingData$normalizedReads = predictReads.round(sequncingData,W,all.catg);
}
readsDataForEM = sequncingData$normalizedReads
areaDataForEM = areas;

lAssignmentResults = par.searchBestAssignemt(reads=readsDataForEM,fish_areas=areaDataForEM,max_trials=EM.trails,initalSeed=initalSeed,seedsNumber=seedsNumber,isL2 = isL2);
assignment = lAssignmentResults[[1]];
lhood = lAssignmentResults[[2]];
seed = lAssignmentResults[[3]];
steps = lAssignmentResults[[4]];
d.assignment = as.data.frame(assignment);
d.assignment$IDs = sequncingData$fullName[assignment[,1]];
d.assignment.areas = ddply(d.assignment,"IDs",summarise,areas = sum(fish_areas));
assignment[,1] = sequncingData$fullName[assignment[,1]];
sequncingData$IDs = sequncingData$fullName;
sequncingData$readsFrac = sequncingData$normalizedReads / sum(sequncingData$normalizedReads);
comb = merge(sequncingData,d.assignment.areas,by="IDs");
dataToReport = cbind(comb$IDs,comb$areas,comb$readsFrac,comb$normalizedReads,comb$Reads);
colnames(dataToReport) = c("IDs","areas","reads_frac","normalized_reads","raw_reads");
writeLines(text=paste("lhood = ",lhood,"\nseed = ",seed,"\ninitalSeed = ",initalSeed,"\nsteps = ",steps,sep =""),con=paste(paste(dirPath,sampleNumber,sep=""),suffixName,"scroing","seed",initalSeed,sep="_"));
write.table(dataToReport,paste(paste(dirPath,sampleNumber,sep=""),suffixName,"reads_areas","scroing","seed",initalSeed,sep="_"));
write.table(assignment,paste(paste(dirPath,sampleNumber,sep=""),suffixName,"fullAssignment","scroing","seed",initalSeed,sep="_"));
write.table(table(assignment[,1]),paste(paste(dirPath,sampleNumber,sep=""),suffixName,"_compositionAssignment","scroing","seed",initalSeed,sep="_"))
