# run this script typing in your terinal:
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
