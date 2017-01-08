setwd("/working_path/sensitivityAnalyses")


########################
### 1. Subsampling  ####
########################


## Load libraries and functions
library(ape)
library(partitions)
library(phrapl)


## Define arguments

## Name of files(without path)
inputAssignment<-"Pleth_align.txt"
inputTrees<-"Pleth_bestTree.tre"       ##Tree file in NEWICK format
subsamplesPerGene<-10 


nloci<-5
popAssignments<-list(c(3,3))

currentAssign<-read.table(paste(getwd(),"/input/",inputAssignment,sep=""), header=TRUE)
currentTrees<-read.tree(paste(getwd(),"/input/",inputTrees,sep=""))
#currentAssign<-read.table(paste(path.package("phrapl"),"/extdata/Pleth_align.txt.txt",sep=""), header=TRUE,stringsAsFactors=FALSE)
#currentTrees<-read.tree(paste(path.package("phrapl"),"/extdata/Pleth_bestTree.tre",sep=""))

    
## Do subsampling
observedTrees<-PrepSubsampling(assignmentsGlobal=currentAssign,observedTrees=currentTrees,
popAssignments=popAssignments,subsamplesPerGene=subsamplesPerGene,outgroup=FALSE,outgroupPrune=FALSE)

## Save subsampled Trees
save(list="observedTrees",file=paste(getwd(),"/input/phraplInput_Pleth.rda",sep=""))


## Get subsample weights for phrapl
subsampleWeights.df<-GetPermutationWeightsAcrossSubsamples(popAssignments=popAssignments,observedTrees=observedTrees)

load(paste(getwd(),"/input/MigrationArray_2pop_3K.rda",sep=""))
migrationArrayMap<-GenerateMigrationArrayMap(migrationArray)

save(list=c("observedTrees","subsampleWeights.df","migrationArrayMap"),file=paste(getwd(),"/input/phraplInput_Pleth.rda",sep=""))


######################
### 2. GridSearch  ###
######################
#######################Subsampling was done before

## Load libraries
library(ape)
library(partitions)
library(lattice)
library(polynom)
library(gmp)
library(rgenoud)
library(parallel)
library(optimx)
library(igraph)
library(numDeriv)
library(nloptr)
library(Matrix)
library(rgl)
library(RColorBrewer)
library(base)
library(phrapl)

## Load input files (from subsampling) if the run is starting here
load(paste(getwd(),"/input/phraplInput_Pleth.rda",sep=""))                 #load the input file 
load(paste(getwd(),"/input/MigrationArray_2pop_3K.rda",sep=""))
#load(paste(path.package("phrapl"),"/extdata/phraplInput_Pleth.rda",sep=""), header=TRUE,stringsAsFactors=FALSE)
#load(paste(path.package("phrapl"),"/extdata/MigrationArray_2pop_3K.rda",sep=""))


## Specify search details

#modelRange<-c(1:length(migrationArray))        
modelRange<-c(1:5)
popAssignments<-list(c(3,3))  
nTrees<-100      
subsamplesPerGene<-10
totalPopVector<-list(c(4,4))     ## total number of indvs per pop

## Run search and keep track of the time
startTime<-as.numeric(Sys.time())

result<-GridSearch(modelRange=modelRange,
	migrationArray=migrationArray,
	migrationArrayMap=migrationArrayMap,           # do not required in the new version CHECK
	popAssignments=popAssignments,
	nTrees=nTrees,
	#	msPath="/Users/ariadnamoralesgarcia/msdir/ms",
	#	comparePath="/Library/Frameworks/R.framework/Versions/3.1/Resources/library/phrapl/extdata/comparecladespipe.pl",
	observedTree=observedTrees,
	subsampleWeights.df=subsampleWeights.df,
	print.ms.string=TRUE,
	print.results=TRUE,
	debug=TRUE,return.all=TRUE,
	collapseStarts=c(0.30,0.58,1.11,2.12,4.07),
	migrationStarts=c(0.10,0.22,0.46,1.00,2.15),
	subsamplesPerGene=subsamplesPerGene,
	totalPopVector=totalPopVector,
	popScaling=c(0.25, 1, 1, 1, 1),
	print.matches=TRUE)

#Print summary results to end of Rout file
print(result[[1]])

#Make dedicated grid list
gridList<-result[[1]]

#Get elapsed time
stopTime<-as.numeric(Sys.time()) #stop system timer
elapsedSecs<-stopTime - startTime #elapsed time in hours
elapsedHrs<-(elapsedSecs / 60) / 60 #convert to hours
elapsedDays<-elapsedHrs / 24 #convert to days

#Save the workspace from first grid analysis and create output folders
system(paste("mkdir ", getwd(),  "/results", sep=""))
save(list=ls(), file=paste(paste(getwd(),"/results/Pleth_",min(modelRange),"_",max(modelRange),".rda",sep="")))
system(paste("mkdir ", getwd(),  "/results/RoutFiles", sep=""))
system(paste("mv ", getwd(), "/scripts/1.Subsampling_GridSearch_Post.Rout ", getwd(), "/results/RoutFiles/1.Subsampling_GridSearch_Post.Rout", sep=""))
system(paste("rm ", getwd(), "/scripts/1.Subsampling_GridSearch_Post.R.out", sep=""))


########################
### 3. Post-process  ###
########################

## If different models (from the same migrationArray file) were run separatly, the output can be concatenated.
totalData<-ConcatenateResults(rdaFilesPath=paste(getwd(),"/results/",sep=""),rdaFiles=NULL,addAICweights=TRUE,addTime.elapsed=FALSE)
modelAverages<-CalculateModelAverages(totalData, parmStartCol=9)

## Save output to a txt file
write.table(totalData, file=paste(getwd(),"/results/totalData.txt",sep=""), sep="\t", row.names=FALSE, quote=FALSE)

## Plot the best model
#PlotModel(migrationArray[[1]], taxonNames=c("S","N"))
#PlotModel(migrationArray[[6]], taxonNames=c("S","N"))
#PlotModel(migrationArray[[3]], taxonNames=c("S","N"))