# run fisrt part this script typing in your terinal:
#R CMD BATCH 1.Subsampling_GridSearch_Post.R > 1.Subsampling_GridSearch_Post.R.out

## Create output dirs
	system(paste0("mkdir ", getwd(),  "/results"))
	system(paste0("mkdir ", getwd(),  "/results/RoutFiles"))
	

## Load libraries and functions
library(ape)
library(partitions)
library(phrapl)

## migrationArray
load(url("https://github.com/ariadnamorales/phrapl-manual/raw/master/data/sensitivityAnalyses/example/input/MigrationArray_2pop_3K.rda"))

## 
currentAssign<-read.table(file="https://raw.githubusercontent.com/ariadnamorales/phrapl-manual/master/data/sensitivityAnalyses/example_with_outputFiles/input/Pleth_align.txt", head=TRUE)

currentTrees<-ape::read.tree("https://raw.githubusercontent.com/ariadnamorales/phrapl-manual/master/data/sensitivityAnalyses/example_with_outputFiles/input/Pleth_bestTree.tre")

########################
### 1. Subsampling  ####
########################

## Define arguments
subsamplesPerGene<-10 
nloci<-5
popAssignments<-list(c(3,3))

    
## Do subsampling
observedTrees<-PrepSubsampling(assignmentsGlobal=currentAssign,observedTrees=currentTrees,
popAssignments=popAssignments,subsamplesPerGene=subsamplesPerGene,outgroup=FALSE,outgroupPrune=FALSE)

### You will see this error message:
### Warning messages:
#		1: In PrepSubsampling(assignmentsGlobal = currentAssign, observedTrees = currentTrees,  :
#	  Warning: Tree number  1 contains tip names not included in the inputted assignment file. These tips will not be subsampled.
### This is because we exclude a few individuals from the Assignment File to reduce computational time for the sake of this tutorial.

## Get subsample weights for phrapl
subsampleWeights.df<-GetPermutationWeightsAcrossSubsamples(popAssignments=popAssignments,observedTrees=observedTrees)

## Save subsampled Trees and weights
save(list=c("observedTrees","subsampleWeights.df"),file=paste0(getwd(),"/phraplInput_Pleth.rda"))
######################
### 2. GridSearch  ###
######################

## Search details
modelRange<-c(1:5)
popAssignments<-list(c(3,3))  
nTrees<-100      
subsamplesPerGene<-10
totalPopVector<-list(c(4,4))     			## total number of indvs per pop
popScaling<-c(0.25, 1, 1, 1, 1)					## IMPORTANT: one of these loci is haploid

## Run search and keep track of the time
startTime<-as.numeric(Sys.time())

result<-GridSearch(modelRange=modelRange,
	migrationArray=migrationArray,         
	popAssignments=popAssignments,
	nTrees=nTrees,
	observedTree=observedTrees,
	subsampleWeights.df=subsampleWeights.df,
	print.ms.string=TRUE,
	print.results=TRUE,
	debug=TRUE,return.all=TRUE,
	collapseStarts=c(0.30,0.58,1.11,2.12,4.07),
	migrationStarts=c(0.10,0.22,0.46,1.00,2.15),
	subsamplesPerGene=subsamplesPerGene,
	totalPopVector=totalPopVector,
	popScaling=popScaling,					## IMPORTANT: one of these loci is haploid
	print.matches=TRUE)

# Grid list for Rout file
gridList<-result[[1]]

#Get elapsed time
stopTime<-as.numeric(Sys.time()) #stop system timer
elapsedSecs<-stopTime - startTime #elapsed time in hours
elapsedHrs<-(elapsedSecs / 60) / 60 #convert to hours
elapsedDays<-elapsedHrs / 24 #convert to days


#Save the workspace and cleaning
save(list=ls(), file=paste0(getwd(),"/results/Pleth_",min(modelRange),"_",max(modelRange),".rda"))
system(paste0("mv ", getwd(), "/1.Subsampling_GridSearch_Post.Rout ", getwd(), "/results/RoutFiles/1.Subsampling_GridSearch_Post.Rout"))
system(paste0("rm ", getwd(), "/1.Subsampling_GridSearch_Post.R.out"))	


##----------> In tutorial we restart here  ---> now the files are assumed to be in your computer

################################
### 4. Sensitivity analyses  ###
################################

## Load libraries
library(phrapl)
library(gplots)

## Load modified function
source(paste(getwd(), "/scripts/GenerateSetLoci_cum8_v3.R", sep=""))

## Load phrapl input (subsampled dataset)
load(paste(getwd(),"/input/phraplInput_Pleth.rda", sep=""))

## Load file with models and store them in an object named "migrationArray"
load(paste(getwd(),"/input/MigrationArray_2pop_3K.rda",sep=""))

## Create output folders
system(paste("mkdir ", getwd(), "/results/subsets",sep=""), ignore.stderr=TRUE)

#########################################################
### Define arguments (some of these were already defined)
nloci<-5
NintervalLoci<-1
modelRange<-c(1:5)
lociRange<-c(1:5)
subsamplesPerGene<-10
popAssignments<-list(c(3,3))
collapseStarts<-c(0.3, 0.58, 1.11, 2.12, 4.07)   
migrationStarts<-c(0.1, 0.22, 0.46, 1, 2.15) 
n0multiplierStarts<-NULL
setCollapseZero<-NULL
nTrees<-100
dAIC.cutoff<-2
nEq<-100
subsampleWeightsVec=subsampleWeights.df[[1]][,1]
migrationArrayShort<-migrationArray[modelRange]

#####################################
### Load files with match vectors ###
rdaFilename<-paste(getwd(),'/results/Pleth_1_5.rda', sep="")

system(paste("grep matches -A1 ",getwd(),"/results/RoutFiles/1.Subsampling_GridSearch_Post.Rout | grep -v matches | grep 0 > ",getwd(),"/results/RoutFiles/matches_2.RunGridSearch.Rout.txt", sep=""))
RoutFilename<-read.table(file=paste(getwd(),'/results/RoutFiles/matches_2.RunGridSearch.Rout.txt', sep="") ,skip=1)

#####################################################
### Recalculate model averages in subsets of loci ###

# Command line to run function GenerateSetLoci
GenerateSetLoci(lociRange=lociRange,NintervalLoci=NintervalLoci,RoutFilename,rdaFilename,migrationArray=migrationArray,modelRange=modelRange,subsamplesPerGene=subsamplesPerGene,collapseStarts=collapseStarts,migrationStarts=migrationStarts,n0multiplierStarts=n0multiplierStarts,setCollapseZero=setCollapseZero,cumulative=TRUE,nTrees=nTrees,dAIC.cutoff=2,nEq=nEq, subsampleWeightsVec=subsampleWeightsVec)


#########################################
### Postprocess RDA files per subset  ###

# Set paths to input and output files
pathSubsets<-paste(getwd(), "/results/subsets/", sep="")              # Files created by GenerateSetLoci. Do not forget "/" at the end
pathtotalData_subsets<-paste(getwd(), "/results/subsets/", sep="") # Do not forget "/" at the end

# Loop to analyze output files ----> Will generate a totalData file per subset per model
# You may need to create one folder per subset and move files to each folder
for(i in 1:5){                    ## ------> Be careful to specify the total number of subsets (in this case 101)
    totalData<-list()
    totalData<-ConcatenateResults(rdaFilesPath=paste(getwd(), "/results/subsets/subset",i,"/",sep=""),rdaFiles=NULL, migrationArray, rm.n0=TRUE, longNames = TRUE, outFile=NULL,addAICweights=TRUE,addTime.elapsed=FALSE, nonparmCols = 4)
	#modelAverages<-CalculateModelAverages(totalData,parmStartCol=9)
    write.table(totalData, file=paste(pathtotalData_subsets,"totalData_subset",i,".txt",sep=""), sep="\t", row.names=FALSE)
}

########################
### 5. Plot subsets  ###
########################

pathInputFiles<-paste(getwd(), "/results/subsets/",sep="")

totalData_list<-list.files(pathInputFiles, pattern="*.txt", full.names=FALSE)

totalData_files<-c()
for(x in 1:length(totalData_list)){
    totalData_files<-append(totalData_files,list(read.table(paste(pathInputFiles,totalData_list[x],sep=""),header=TRUE, fill=TRUE, sep="\t")))
}

#Sort each dataframe by model
totalData_filesSorted<-c()
for(x in 1:length(totalData_list)){
	totalData_filesSorted<-append(totalData_filesSorted,list(totalData_files[[x]][order(totalData_files[[x]]$models),]))
}

#Create a dataframe with all models and their wAIC
model_wAIC<-data.frame("models"=c(1:5),row.names=NULL)
for(x in 1:length(totalData_list)){
    model_wAIC[,paste("Subset ",x,sep="")]<-totalData_filesSorted[[x]][,"wAIC"]
}

########################
## Define plot settings

## Row Labels
row.label<-c()
setLoci<-0
for(i in 1:length(totalData_list)){
    setLoci<-setLoci+8      ##loci interval
    row.label<-append(row.label, c(paste(setLoci, " Loci",sep="")))
}

## Row labels ---> 
row.labels <- rep("", 5)

## Col Labels
col.label<-seq(1:5)

## Color scale
bks<-c(seq(0,0.99,by=0.05))
lm.palette <- colorRampPalette(c("light yellow","orange", "tomato"), space = "rgb")

## Possition of the plot and legend
lmat <- rbind(c(0,3,3,3),c(2,1,1,1),c(0,0,0,4))
lwid <- c(0.2,1,1,1)
lhei <- c(0.3,2.1,0.5)

pdf(file=paste(getwd(),"/results/wAICinSubsets.pdf", sep=""), width=9, height=7)
heatmap.2(as.matrix(t(model_wAIC[2:length(model_wAIC)])), lmat = lmat, lhei=lhei,lwid=lwid,
	Rowv=FALSE,Colv=FALSE, dendrogram="none", trace="none",
	margins=c(2,12),
	col=lm.palette(length(bks)-1), breaks=bks,
	#colsep=1:5,
	rowsep=1:(length(model_wAIC)-1), sepcolor="white",
	#labRow=row.labels,
	#labCol=col.label,
	key = TRUE,
	#keysize =1.0,
	density.info=c("density"),
	denscol=tracecol,
	key.xlab="wAIC",
	key.par=list(mar=c(4,4,4,4)),
	main="wAIC for each model per locus")
dev.off()




