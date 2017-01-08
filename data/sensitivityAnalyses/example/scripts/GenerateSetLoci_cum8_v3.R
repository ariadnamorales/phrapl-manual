#########FUNCTION

#Outputs rda files for each subset of loci
GenerateSetLoci<-function(lociRange,NintervalLoci,RoutFilename,rdaFilename,migrationArray,modelRange=c(1:length(migrationArray)),subsamplesPerGene,collapseStarts,migrationStarts,n0multiplierStarts,setCollapseZero=NULL,cumulative=FALSE,nTrees,dAIC.cutoff=2,nEq=nEq,subsampleWeightsVec){

	#Get matching vectors from Rout
    matchesTEMP<-RoutFilename[2]
    matches<-c()
	for(i in 1:nrow(matchesTEMP)){
		matches<-append(matches,as.vector(matchesTEMP[i,]))
    }
		
				
	#Load rda file to get grid
	load(rdaFilename)
	migrationArray<-migrationArrayShort			##Make sure the correct migrationArray list is used (a complete list is loaded with rdaFilename)
	
	##Sort GridList (sorted by AIC in rda, but want to be sorted by grid order)
	for(i in 1:length(gridList)){
		gridList[[i]]<-gridList[[i]][order(as.numeric(row.names(gridList[[i]]))),]
	}

    ##Creates a vector that repeats grid length per model
    sep_intervalsGrid<-c()
    for(modelID in 1:length(migrationArray)){
        sep_intervalsGrid<-append(sep_intervalsGrid,c(as.numeric(rep(modelID,nrow(gridList[[modelID]])))))
    }

    ##Splits the matching vector into a list such that each model is an item in the list
    tableVectorByGridLength<-split(matches,sep_intervalsGrid)
	

    ##################Print RDAs for each subset of loci

    #For each subset of loci...
    treesVec<-sequence(max(lociRange) * subsamplesPerGene)
    counterBegin<-1
    counterEnd<-NintervalLoci * subsamplesPerGene
    gridListOriginal<-gridList
    for(locusSubset in 1:(max(lociRange)/NintervalLoci)){ #Get positions of each locus subset in the matching vector
        gridList<-gridListOriginal
        currentLociRange<-treesVec[counterBegin:counterEnd]
        if(cumulative==FALSE){
            counterBegin<-counterBegin + (NintervalLoci * subsamplesPerGene)
        }
    counterEnd<-counterEnd + (NintervalLoci * subsamplesPerGene)

    #For each model...
    for(model in 1:length(migrationArray)){

        #For each parameter combination in the grid
        for(gridpoints in 1:length(tableVectorByGridLength[[model]])){
            totalVector<-as.numeric(strsplit(tableVectorByGridLength[[model]][gridpoints]," ")[[1]])
            totalVectorSubsampled<-totalVector[currentLociRange]

            lnLValue<-ConvertOutputVectorToLikelihood.1sub(outputVector=totalVectorSubsampled,
            popAssignments=popAssignments,nTrees=nTrees,subsamplesPerGene=subsamplesPerGene,totalPopVector=totalPopVector,
            summaryFn="mean",nEq=nEq,subsampleWeightsVec=subsampleWeightsVec)

            #Replace the old AIC value in the grid with the new value (note: subtract from K for n0multiplier and any
            #collapses set to be zero.
            gridList[[model]][gridpoints,1]<-2*(-lnLValue + (KAll(migrationArray[[model]]) - length(setCollapseZero) - 1))
        }
    }

    #Re-sort gridList by AIC
    for(i in 1:length(gridList)){
        gridList[[i]]<-gridList[[i]][order(gridList[[i]]$AIC),]
    }

    #Now recalculate overall results object and parameter values object
    overall.results<-ExtractGridAICs(result=gridList,migrationArray=migrationArray,modelRange=modelRange)   #,setCollapseZero=setCollapseZero
    parameters<-ExtractGridParameters(migrationArray=migrationArray,result=gridList,popVector=popAssignments[[1]]) #,dAIC.cutoff=dAIC.cutoff
    result<-list("AIC.Grid"=gridList,"overall.results"=overall.results,"parameters"=parameters[[1]],"parameterIndexes"=parameters[[2]])

    #Save all this to new rda - you get a new rda for each subset of loci
    save(list=c("gridList","overall.results","parameters","result","nTrees"),file=paste(rdaFilename,"_1c_subset",locusSubset,".rda",sep=""))

    system(paste("mkdir ", getwd(), "/results/subsets",sep=""), ignore.stderr=TRUE)
    
    for(locusSubset in 1:(max(lociRange)/NintervalLoci)){ #Get positions of each locus subset in the matching vector
        system(paste("mkdir ", getwd(), "/results/subsets/subset",locusSubset,sep=""), ignore.stderr=TRUE)
  		system(paste("mv ",getwd(), "/results/*_1c_subset",locusSubset,".rda ", getwd(), "/results/subsets/subset",locusSubset, sep=""),ignore.stderr=TRUE)
	}
    }
}

