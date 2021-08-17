# download and run this script typing in your terinal:
#
# wget https://github.com/ariadnamorales/phrapl-manual/raw/master/data/tutorial2.2/scripts/2.2.buildModels_loadData_Subsampling_GridSearch_Post.R
# R CMD BATCH 2.2.buildModels_loadData_Subsampling_GridSearch_Post.R > 2.2.buildModels_loadData_Subsampling_GridSearch_Post.Rout


## R script with all steps required for Tutorial2-2
##
## How to create and compare models of ancestral gene flow vs secondary contact. 
##
## 
##
################################
######---> R TUTORIAL STARTS HERE


## Create output dirs
	system(paste0("mkdir ", getwd(),  "/results"))

## Load libraries
library(phrapl)	

###########################
### 1. Generate models ####
###########################

## get migrationArray with 2 models (no gene flow and asymmetric gene flow)
load(url("https://github.com/ariadnamorales/phrapl-manual/raw/master/data/tutorial2.2/models/MigrationArray_2Pops_2models_noMig_assymMig.rda"))

## Verify length of migrationArray ---> We start with two models
length(migrationArray)

## Add model of ancenstral migration
migrationArray[[3]]<-AddEventToMigrationArray(migrationArray[[2]], 1, migrationMat=0)[[1]]

## Add model of secondary contact
migrationArray[[4]]<-AddEventToMigrationArray(migrationArray[[1]], 1, migrationMat=1)[[1]]

## Verify length of migrationArray ---> We should have now four models
length(migrationArray)

## Save new model set
migrationArrayMap<-GenerateMigrationArrayMap(migrationArray) 
save(migrationArray,migrationArrayMap, file="MigrationArray_2Pops_4models_noMig_assymMig_ancMig_secCont.rda")



#########################################################
### 2. load exaomple inputData  -- aleady subsampled ####
#########################################################

## Load input files & sumsampled data --> required objects: "observedTrees", "subsampleWeights.df"
load(url("https://github.com/ariadnamorales/phrapl-manual/raw/master/data/tutorial2.2/Data/tutorial_2.2_phraplInput.rda"))


##########################
### 3. Run GridSearch  ###    ----> IMPORTANT: input arguments are not the same for all models 
##########################					    you cannot use exaclty the same arguments in 'GridSearch'
##
## Ideally you should run 'GridSearch' with one script per model, and can be run in parallel
## which will reduce computational time, specially when you have hundreds or thousands of loci.
##
## ------- *** If you run everything in one script, like in this example (easier to explain the tutorial)
## ------- *** 
## ------- *** in R, the output of 'GridSearch' should be stored in an object with the same name for all models, "result" in this case
## ------- *** 
## ------- *** BUT each model output should be saved in a separate "rda" file (with different names)
## ------- *** 
## ------- *** AND you should remove that object before you run a new 'GridSearch'
## ------- *** 
## ------- *** this will avoid any confussion if a run for a scecific model crashes
## ------- *** because if the object "result" is returned at the end, that means the 'GridSearch'
## ------- *** finish succesfully.


## ************************************************************************************ ##
## **** First we run models 3 and 4 ---> the ones with an additional time interval **** ##


#####################################
### model 3 ---> Ancestral migration

modelRange=3			## IMPORTANT: make sure you specify the correct model

result<-GridSearch(modelRange=modelRange,
migrationArray=migrationArray,
popAssignments=list(c(2, 2)),
nTrees=3,
observedTree=observedTrees,
subsampleWeights.df=subsampleWeights.df,
print.ms.string=TRUE,
print.results=TRUE,
debug=TRUE,return.all=TRUE,
collapseStarts=c(0.3, 0.58, 1.11, 2.12),
migrationStarts=c(0.46, 1, 2.15, 4.64),
setCollapseZero=NULL,
subsamplesPerGene=100,
totalPopVector=list(c(28, 50)),
print.matches=TRUE,
addedEventTimeAsScalar=TRUE,		## IMPORTANT: indicate if time interval is relative or absolute
addedEventTime=0.5)					## IMPORTANT: indicate when the event (e.g. gene flow) started or stopped.

#Make dedicated grid list
gridList<-result[[1]]

#Save the workspace from first grid analysis
save("gridList", "result", file=paste0("results/tutorial_2.2_output_phrapl_model", modelRange, ".rda"))

## remove output objects to start clean in the next 'GridSearch'
rm(modelRange)
rm(result)
rm(gridList)


###################################
### model 4 ---> Secondary contact

modelRange=4			## IMPORTANT: make sure you specify the correct model

result<-GridSearch(modelRange=modelRange,
migrationArray=migrationArray,
popAssignments=list(c(2, 2)),
nTrees=3,
observedTree=observedTrees,
subsampleWeights.df=subsampleWeights.df,
print.ms.string=TRUE,
print.results=TRUE,
debug=TRUE,return.all=TRUE,
collapseStarts=c(0.3, 0.58, 1.11, 2.12),
migrationStarts=c(0.46, 1, 2.15, 4.64),
setCollapseZero=NULL,
subsamplesPerGene=100,
totalPopVector=list(c(28, 50)),
print.matches=TRUE,
addedEventTimeAsScalar=TRUE,		## IMPORTANT: indicate if time interval is relative or absolute
addedEventTime=0.5)					## IMPORTANT: indicate when the event (e.g. gene flow) started or stopped.

#Make dedicated grid list
gridList<-result[[1]]

#Save the workspace from first grid analysis
save("gridList", "result", file=paste0("results/tutorial_2.2_output_phrapl_model", modelRange, ".rda"))

## remove output objects to start clean in the next 'GridSearch'
rm(modelRange)
rm(result)
rm(gridList)

## ******************************************************************** ##
## **** Now we run models 1 and 2 ---> no additional time interval **** ##
## ****															   **** ##
## **** This is the "standard" way of running GridSearch           **** ##
## **** you may already be familiar. See tutorial 1                **** ##
## ****															   **** ##
## ******************************************************************** ##


##############################
### model 1 ---> No gene flow

modelRange=1 		## IMPORTANT: make sure you specify the correct model

result<-GridSearch(modelRange=modelRange,
migrationArray=migrationArray,
popAssignments=list(c(2, 2)),
nTrees=3,
observedTree=observedTrees,
subsampleWeights.df=subsampleWeights.df,
print.ms.string=TRUE,
print.results=TRUE,
debug=TRUE,return.all=TRUE,
collapseStarts=c(0.3, 0.58, 1.11, 2.12),
migrationStarts=c(0.46, 1, 2.15, 4.64),
setCollapseZero=NULL,
subsamplesPerGene=100,
totalPopVector=list(c(28, 50)),
print.matches=TRUE,
addedEventTimeAsScalar=TRUE,
addedEventTime=NULL)				## NOTE: In this type of model, this argument does not exist.

#Make dedicated grid list
gridList<-result[[1]]

#Save the workspace from first grid analysis
save("gridList", "result", file=paste0("results/tutorial_2.2_output_phrapl_model", modelRange, ".rda"))

## remove output objects to start clean in the next 'GridSearch'
rm(modelRange)
rm(result)
rm(gridList)

#####################################
### model 2 ---> Constant gene flow

modelRange=2	## IMPORTANT: make sure you specify the correct model

result<-GridSearch(modelRange=modelRange,
migrationArray=migrationArray,
popAssignments=list(c(2, 2)),
nTrees=3,
observedTree=observedTrees,
subsampleWeights.df=subsampleWeights.df,
print.ms.string=TRUE,
print.results=TRUE,
debug=TRUE,return.all=TRUE,
collapseStarts=c(0.3, 0.58, 1.11, 2.12),
migrationStarts=c(0.46, 1, 2.15, 4.64),
setCollapseZero=NULL,
subsamplesPerGene=100,
totalPopVector=list(c(28, 50)),
print.matches=TRUE,
addedEventTimeAsScalar=TRUE,
addedEventTime=NULL)				## NOTE: In this type of model, this argument does not exist.

#Make dedicated grid list
gridList<-result[[1]]

#Save the workspace from first grid analysis
save("gridList", "result", file=paste0("results/tutorial_2.2_output_phrapl_model", modelRange, ".rda"))

## remove output objects to start clean in the next 'GridSearch'
rm(modelRange)
rm(result)
rm(gridList)


########################
### 4. Post-process  ###  ----> This step is basically the same in all tutorials
########################

## concatenate results and calculate model averaged
totalData<-ConcatenateResults(rdaFilesPath=paste0(getwd(),"/results/"),rdaFiles=NULL,addAICweights=TRUE,addTime.elapsed=FALSE)
modelAverages<-CalculateModelAverages(totalData, parmStartCol=9)

## Save output table with wAIC values
write.table(totalData, file=paste0(getwd(),"/results/totalData.txt"), sep="\t", row.names=FALSE, quote=FALSE)

## Save relevant output objects in "rda" file
save(list=ls(), file=paste0(getwd(),"/tutorial_2.2_outputSummary.rda"))

######---> R TUTORIAL ENDS HERE