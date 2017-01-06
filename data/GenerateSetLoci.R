#This function takes output from a phrapl run performed using multiple loci and re-calculates

#AICs and parameter estimates for different subsets of loci. This can either be done for independent

#subsets of loci (cumulative=FALSE) or for accumulating subsets of loci (cumulative=TRUE, in which

#the second subset is added to the first subset, and then the third subset is added to these, etc).

#lociRange gives the number of loci in the original analysis and NintervalLoci gives the desired number 

#of loci in each subset, which must be a multiple of the total number of loci. If the models in

#the original analysis are only a subset of models in the migrationArray, this modelRange must be

#specified as well. Output of this function is an rda file for each locus subset (numbered in order

#with a numerical suffix, i.e., _1, _2, _3, etc). Note that subsets are taken in order from the original

#output files.

GenerateSetLoci<-function(lociRange,NintervalLoci,RoutFilename,rdaFilename,migrationArray,

	modelRange=c(1:length(migrationArray)),subsamplesPerGene,popAssignments,collapseStarts,

	migrationStarts,n0multiplierStarts,setCollapseZero=NULL,cumulative=FALSE,nTrees,

	dAIC.cutoff=2,nEq=nEq,totalPopVector=totalPopVector){



	#Get matching vectors from Rout

    matchesTEMP<-system(paste("grep matches -A1 ",RoutFilename," | grep [0123456789]",sep=""),intern=TRUE)

    matches<-c()

    for(i in 1:length(matchesTEMP)){

    	matches<-append(matches,strsplit(matchesTEMP[i],"\"")[[1]][2])

    }



	#Load rda file to get grid

	load(rdaFilename)

	

	##Sort GridList (sorted by AIC in rda, but want to be sorted by grid order)

	for(i in 1:length(gridList)){

		gridList[[i]]<-gridList[[i]][order(as.numeric(row.names(gridList[[1]]))),]

	}



    ##Creates a vector that repeats grid length per model

    sep_intervalsGrid<-c()

    for(modelID in 1:length(modelRange)){

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

                	summaryFn="mean",nEq=nEq)

                

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

        overall.results<-ExtractGridAICs(result=gridList,migrationArray=migrationArray,

        	modelRange=modelRange,setCollapseZero=setCollapseZero)

        parameters<-ExtractGridParameters(migrationArray=migrationArray,result=gridList,

			popVector=popAssignments[[1]],dAIC.cutoff=dAIC.cutoff)

		result<-list("AIC.Grid"=gridList,"overall.results"=overall.results,

			"parameters"=parameters[[1]],"parameterIndexes"=parameters[[2]])

				

		#Save all this to new rda - you get a new rda for each subset of loci

		save(list=c("gridList","overall.results","parameters","result","nTrees"),

			file=paste(rdaFilename,"_subset",locusSubset,".rda",sep=""))

	}

}
