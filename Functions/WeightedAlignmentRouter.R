


WeightedAlignementRouter <- function(Blastlist){


Datafinal <- Blastlist


RecursiveListlevel <- function(x) {
  #Must input a non-empty list
  if( identical(length(x),as.integer(0)) ){
    stop("Empty call, length of input is 0")
  }
  
  xy <- 1
  y <- FALSE
  while(y == FALSE){
    
    bb <- eval(parse(text=paste0("x", paste0(eval(parse(text="rep(paste0(\"[[1]]\"),xy)")),collapse="") )))
    if( is.character(bb) ){
      y <- TRUE
      
    }
    xy <- xy + 1 
  }
  xy-1
  
}



Bleh <- list()
for(i in 1:length(levels(as.factor(Datafinal$Subjacc))) ){ 
  
  
  Data <- subset(Datafinal, Datafinal$Subjacc == levels(as.factor(Datafinal$Subjacc))[i] )
  
  #Add a sense check first and flip everything to pos sense, then add a flag for anti sense /sense
  #Sense = 1 , antisense = 0
  
  Data$Sense <- as.numeric(Data$Sstart < Data$Send)
  
  
  for( ij in 1:dim(Data)[1]){
    #Rowwise, flip sstart and ssend if data$sense is equal to 0
    TempVar1 <- cbind(Data$Sstart[ij], Data$Send[ij]  )
    
    if( !identical(Data$Sense[ij],1) ){
      Data$Sstart[ij] <- TempVar1[2] 
      Data$Send[ij] <- TempVar1[1]
      
    }
    
    
  }
  
  
  
  #QC: are there overlapping alignments? 
  WorkingD <- Data
  # WorkingD <- Data[,c(seq(1,9,1),11)]
  WorkingD <- WorkingD[c(order(WorkingD$Sstart,WorkingD$Evalue)),]
  #Are there multiple results for Data?
  if( dim(Data)[1]>=2 ){ 
    
    
    
    
    #Check if results are from multiple contigs or etc
    
    
    #For loop that checks first object against all, then checks 2nd object against all etc
    
    #We will assume everything is pos sense
    
    #Create matrix of overlaps, 1 is overlap, 0 is not
    OLMat <- data.frame()
    for(z in 1:(dim(WorkingD)[1]-1) ){
      
      #Determine overlap between object z and others:
      
      
      OLcheck <- !c( (WorkingD$Sstart[z] > WorkingD$Send[(z+1):length(WorkingD$Send)]) |
                       (WorkingD$Send[z] < WorkingD$Sstart[(z+1):length(WorkingD$Sstart)]) )
      # 
      #Be less strict, allow some slack (2% ?)
      #Instead of binary checks, calculate percentage overlaps, and dissallow if overlap exceeds 2% of alignment length
      
      #This check is not sufficient ?
      if( any( (WorkingD$Send[z] >= WorkingD$Sstart[(z+1):length(WorkingD$Send)])&(WorkingD$Send[z] < WorkingD$Send[(z+1):length(WorkingD$Send)])  )  ){
        #uhhh does this cover every possible overlap type? surely not
        
        # WorkingD$Sstart[z]<WorkingD$Sstart[(z+1):length(WorkingD$Send)] 
        a <- which((WorkingD$Send[z] >= WorkingD$Sstart[(z+1):length(WorkingD$Send)])&(WorkingD$Send[z] < WorkingD$Send[(z+1):length(WorkingD$Send)]) )
        b <- 100*((WorkingD$Send[z]-WorkingD$Sstart[(z+1):length(WorkingD$Send)])[c(a)]/WorkingD$Length[z])
        OLcheck[c(a)] <- !(b < 5)
      }
      
      
      #Sanity check: Are we comparing alignments of the same contig or not
      
      # OLcheck[ which( WorkingD$Queryacc[z] == WorkingD$Queryacc[(z+1):length(WorkingD$Send)] ) ] <- rep(1, dim(WorkingD)[1]-z)
      
      #
      # OLcheck[ which( WorkingD$Queryacc[z] == WorkingD$Queryacc[(z+1):length(WorkingD$Send)] ) ] <- rep(1, length( which( WorkingD$Queryacc[z] == WorkingD$Queryacc[(z+1):length(WorkingD$Send)] ) ) )
      
      #Only results to the right  / or above the diagonal are computed, which is why we can create a n-1 by n matrix of overlaps 
      OLMat <- rbind(OLMat, c( rep(1,z), as.numeric(OLcheck)) )
      
    }
    
    if( identical(TRUE,any(OLMat == 0)) ){
      #Compute best path (ie path of non-overlapping results that give highest sum of bitscore)
      # BestP <- which(WorkingD$Bitscore==(max(WorkingD$Bitscore))) #fix this
      # print(paste(i, "hey"))
      # stop()
      
      
      
      nOLroute <- list()
      nOLroute[[1]] <- (as.character(seq(1,dim(WorkingD)[1])))
      #Compute all routes, by row of OLMat
      for( ja in 1:dim(OLMat)[1]){
        
        #If the sum is 0, the best route is 1 step long
        if( identical(sum(!OLMat[ja,]),as.integer(0)) ){
          
          # TempVar3 <- as.character(ja) #not really useful, since we add 1-step routes later
          # TempVar3 <- as.character("0")
          next
          
          #Sum is 1, the best route is 2 steps long:
        } else if( identical(sum(!OLMat[ja,]),as.integer(1)) ){
          
          TempVar3 <- as.character(paste0(ja,".",which(!OLMat[ja,])))
          #Sum is 2 or above, there may be multiple routes, or >=3 step route(s)
        } else{
          #Max iterations of jaz is sum(!OLMat[ja,])-1  
          
          #Which alignments don't overlap
          TempVar4 <- which(!OLMat[ja,])
          
          #Check overlap between var 1 and 2 of tempvar4
          #Sum of non-overlapping alignments is 2: 
          if( identical(length(TempVar4),as.integer(2)) ){ 
            #Check if theres overlap or non-overlap between the variables of tempvar4
            if( identical(OLMat[TempVar4[1],TempVar4[2]],as.numeric(1)) ){
              #If ^is equal to 1, there is no 3 step path
              TempVar3 <- as.character(paste0(ja,".",which(!OLMat[ja,])))
              
            } else {
              #Previous check is equal to 0, which means there is a 3 step path
              TempVar3 <- paste0(c(ja,which(!OLMat[ja,])),collapse=".")
              
            }
            
            
            
            
          } else if(identical(length(TempVar4)>2,TRUE) ){
            # stop(paste("what",i,ja))
            
            
            TempVar5 <- list()
            # TempVar6 <- list()
            TempVar7 <- list()
            #For every row, compute all "routes"
            #Highest step path is max(rowSums(!OLMat))
            TempVar5[[1]] <- list(as.character(paste0(ja,".",which(!OLMat[ja,]))))
            
            
            # eval(parse(text = "paste0(\"hey\",jl)"  ))
            #
            #Highest possible order of paths
            
            for( jl in 1:(sum(!OLMat[ja,])) ){ 
              # print(paste( "jl =", jl ))
              #any(grepl(0,TempVar5[[jl]]))
              if( identical(length(TempVar5[[jl]][[1]])==0,TRUE) | identical(TempVar5[[jl]][[1]],0)  ){
                # print("hey")
                break
                #print(jl,ji,jc)
              }
              
              TempVar7 <- list()
              #Returns the first level that is not a list:
              
              #This does not work. 
              # length(TempVar5[[jl]])
              abb <- c()
              for( zji in 1:length(TempVar5[[jl]])){
              # abb[[zji]] <- list(TempVar5[[jl]][[zji]])
              abb <- c(abb,TempVar5[[jl]][[zji]])
              }
              if( identical(jl,as.integer(1) )){
              abb <- list(abb)
              } 
              # abb <- paste0("TempVar5[[jl]]",rep("[[1]]",RecursiveListlevel(TempVar5[[jl]])-1))
              # 
              # abb <- eval(parse(text=eval(parse(text="abb"))))
              # 
              #Length of TempVar5[[jl]] doesn't work when tempvar5[[jl]] is 1 list (top level) length((TempVar5[[jl]][[1]])), need tot determine how many layers we need to go before we take length
              for(ji in 1:length( abb )){
                # print(paste( "ji =", ji ))
                TempVar6 <- list()
                
                
                
                #Use Recursivelistlevel to find repeats of [[1]]
                # if( identical(typeof(TempVar5[[jl]][[1]]),"list")  ){
                #   
                #   #str(TempVar5[[2]],max.level = 1)
                #   TempTemp5 <- TempVar5[[jl]][[1]][[ji]]
                # } else { 
                TempTemp5 <- abb[[ji]]
                # }
                
                #Check length of temptemp5 before running this
                
                if( !identical(length(TempTemp5),as.integer(1)) ){
                  # TempVar6 <- list() #already defined as list for every iteration of ji
                  
                  for( jc in 1:(length(TempTemp5)-1) ){
                    # print(paste( "jc =", jc ))
                    
                    #eg. If 1 paths to 2, and 1 paths to 3, check if we can make a 1-2-3 path
                    
                    
                    # (TempVar5[[1]][[1]][jc])
                    #Find the last element of the jc'th tempvar5
                    jcz <- unlist(strsplit((TempTemp5[[jc]]),".",fixed=TRUE))[length(unlist(strsplit((TempTemp5[[jc]]),".",fixed=TRUE)))]
                    
                    #Combine jcz with 
                    
                    #Last element can combine with the following in a 3-step path:
                    # which(!OLMat[jcz,])
                    if( !identical(which(!OLMat[jcz,]),integer(0)) ){
                      # TempVar6[jc] <- (list(as.character(paste0(ja,".",jcz,".",which(!OLMat[jcz,]) ))))
                      TempVar6[[jc]] <- paste0(TempTemp5[jc],".",which(!OLMat[jcz,]))
                      
                    } else {
                      TempVar6[[jc]] <- "0"
                    }
                    
                  }
                }
                
                
                if( identical(length(TempVar6)>=1,TRUE) ){ 
                  #any(grepl(0,TempVar6)) 
                  #
                  if( any(unlist(lapply((TempVar6=="0"),as.logical))) ){
                    #This one removes too much shit??
                    # TempVar6 <- TempVar6[ -which(grepl(0,TempVar6))]
                    TempVar6 <- TempVar6[ c(-which((unlist(lapply((TempVar6=="0"),as.logical))))) ]
                  }
                  if( length(TempVar6)>=1 ){ 
                    TempVar7[ji] <- list(TempVar6)
                  } else {
                    TempVar7[ji] <- (as.numeric(0))
                  }
                  # TempVar7[ji] <- list(TempVar6)
                  
                } else {
                  
                  TempVar7[ji] <- (as.numeric(0)) #next? 
                }
              }
              
              
              
              if( identical(length(TempVar7[[ji]])==0,TRUE) ){
                # TempVar5[[jl+1]] <- as.numeric(0)
                break
                #break?
              } else{ 
                #
                #any(grepl(0,TempVar7))
                if( any(unlist(lapply(TempVar7[[ji]]=="0",as.logical))) ){
                  #What am i doing here??? this is fucking stupid what the fuck this check is fucking bad what the fuck
                  if( identical( sum(unlist(lapply(TempVar7=="0",as.logical))) , length(unlist(TempVar7)))  ){
                    # print("hi")
                    break
                  }else {
                    #TempVar5[[jl+1]] <- TempVar7[ -c(grep(0,TempVar7)) ] #this does not work
                    TempVar5[[jl+1]] <- TempVar7[ c(-which(unlist(lapply(TempVar7=="0",as.logical)))) ]
                  }
                  
                } else { 
                  TempVar5[[jl+1]] <- TempVar7
                }
              }
              
            }
            
            
            TempVar3 <- unlist(TempVar5)
            
            if( any(TempVar3=="0") ){
              TempVar3 <- TempVar3[TempVar3!="0"]
              
            }
            
            #We need to include single "letter" routes, if we haven't already
            #Check and remove non-list elements from TempVar5 
            # any(!unlist(lapply(TempVar5,is.list)))
            
            #todo Check this function on datasets and datasets with randomised order
            
            # # max(rowSums(!OLMat))
            # jc <- length(TempVar4) 
            # jz <- 1
            # TempVar5 <- list()
            # TempVar6 <- list()
            # #jcz <- 1
            # #First find all 2-step paths: (number of paths: p1)
            # #TempVar5[] <- list(as.character(paste0(ja,".",which(!OLMat[ja,]))))
            # 
            # 
            # #Pull data from TempVar5 results 
            # while( jc >= 2 ){
            #   
            #   jc <- sum(as.numeric(!OLMat[TempVar4[jz],]))
            #   
            #   if( jc == 0  ){
            #     # TempVar5 <- list(as.character(unlist(TempVar5)))
            #     break
            #   
            #   }
            #   TempVar5[jz] <- list(paste0(ja,".",TempVar4[jz],".",which((!OLMat[TempVar4[jz],])) ))
            #   
            #   jz = jz+1
            # }
            # 
            # #2-step paths
            # TempVar6[[1]] <- list(as.character(paste0(ja,".",which(!OLMat[ja,]))))
            # #3-step paths
            # TempVar6[[2]] <- list(as.character(unlist(TempVar5)))
            # 
            # 
            # if( identical(length(TempVar4)==as.integer(3),TRUE) ){
            #   for( jc in 1:(length(TempVar4)-1)){
            #     
            #     # TempVar5 <- TempVar4[jc] == TempVar4[ c( (jc+1):length(TempVar4) ) ]
            #     which(!OLMat[TempVar4[jc], ]) %in% TempVar4[c( (jc+1):length(TempVar4)   )]
            #     
            #     
            #     
            #   }
            #   
            # }
            
            #Afterwards run for loop (p1-1) to check for 3-step paths between 2-step paths
            
            #iterate untill all possible n-paths have been exhausted 
            #Mb use while loop and countdown  
            
            
            
            
            
            
            
          }
          
          
          
          
          
          
        }
        
        nOLroute[[ja+1]] <- TempVar3
        
        
      }
      
      #Remove NULL elements
      if( any(unlist((lapply(nOLroute,is.null)))) ){
        nOLroute <- nOLroute[-c(which(unlist((lapply(nOLroute,is.null)))))]
      }
      
      
      #Create a matrix of bitscores using negated OLMat, then calculate row sum
      OLMatsum <- list()
      for( jb in 1:length(nOLroute) ){ 
        #Very simplified, only looks at "1 step" paths
        # TempVar2 <- if( identical(numeric(0), WorkingD$Bitscore[ !OLMat[jb,] ]) ){
        #   0 } else {  max(WorkingD$Bitscore[ !OLMat[jb,] ])  }
        
        # nOLroute[[jb]]
        
        TempVar2 <- strsplit(nOLroute[[jb]],".",fixed=TRUE)
        OLMatsum[[jb]] <- list()
        for( jbz in 1:length(TempVar2)){
          
          OLMatsum[[jb]][[ as.character(nOLroute[[jb]][jbz]) ]] <- list()
          OLMatsum[[ jb ]][[ as.character(nOLroute[[jb]][jbz]) ]] <- sum(WorkingD$Bitscore[ c(as.numeric(TempVar2[[jbz]])) ])
          
          
        }
        #unlist(lapply(OLMatsum,names))
        
        # WorkingD$Bitscore[ ]
        # OLMatsum <- cbind(OLMatsum, sum( WorkingD$Bitscore[  as.numeric(unlist(strsplit(nOLroute[[jb]],".",fixed=TRUE)))   ]) )
        # OLMatsum <- cbind(OLMatsum, sum(WorkingD$Bitscore[jb], TempVar2))
        # OLMatsum <- cbind(OLMatsum, max(WorkingD$Bitscore[jb] + WorkingD$Bitscore[ !OLMat[jb,] ]) )
        # nOLroute <- lks
        
        
      }
      BestP <- unlist(lapply(OLMatsum,as.character))
      names(BestP) <- unlist(lapply(OLMatsum,names))
      # which(as.numeric(BestP)==max(as.numeric(BestP)))
      
      BestP <- c(as.numeric(unlist(strsplit(names(BestP[   which(as.numeric(BestP)==max(as.numeric(BestP)))   ]),".",fixed=TRUE))))
      # BestP <- as.numeric(unlist(strsplit( nOLroute[[ as.numeric(which(OLMatsum == max(OLMatsum))) ]], ".",fixed=TRUE)))
      
      #Find best match, then retrace what the route was
      #Route:
      
    } else {
      #Chose highest bitscore hit
      BestP <- which(WorkingD$Bitscore==(max(WorkingD$Bitscore)))
      
      
      
    }
    
    
    WorkingD <- WorkingD[c(BestP),]
    
    
    #Redefine b
    # b <- dim(WorkingD)[1]
  }
  
  
  ### Remove non base symbols
  #Determine number of rows
  # if( identical(length( dim(WorkingD)[1] ),as.integer(1)) ){
  #   #True for only 1 row ( sequences)
  #   # print(paste("type 1 match for",i))
  #   
  #   test <- sapply(strsplit(Data2$qseq,""),as.list)[ as.numeric(b) ]
  #   test <- lapply(test,as.character)
  #   
  #   a <- unlist(lapply(grepl("-",test),any))
  #   
  #   if(any(a) ){ 
  #     for(j in 1:sum(a)){
  #       
  #       #Which to only get TRUE objects, 
  #       # select the j'd true object
  #       
  #       #Remove hyphens (gaps)
  #       test[[ which(a)[j] ]] <- test[[ which(a)[j] ]][-c(grep("-",test[[ which(a)[j] ]])) ]
  #       
  #       
  #       
  #     }
  #   }
  #   
  # } else if( identical(length( b )>=2,TRUE)   ){
  #   #There's two or more sequences 
  #   
  #   
  #   # print(paste("type 2 match for",i))
  #   test <- sapply(strsplit(Data2$qseq,""),as.list)[ c(as.numeric(b)) ]
  #   test <- lapply(test,as.character)
  #   
  #   a <- unlist(lapply(grepl("-",test),any))
  #   
  #   #Only remove hyphens from sequences that contain hyphens. 
  #   if(any(a)){ 
  #     for(j in 1:sum(a)){
  #       
  #       #Which to only get TRUE objects, 
  #       # select the j'd true object
  #       
  #       #Remove hyphens (gaps)
  #       test[[ which(a)[j] ]] <- test[[ which(a)[j] ]][-c(grep("-",test[[ which(a)[j] ]])) ]
  #       
  #       
  #       
  #     }
  #   }
  #   
  # } else {
  #   #couldn't find match between datafinal and seqnames most likely?
  #   # print(paste("no match for",i))
  #   
  # }
  # 
  # 
  # 
  # 
  # # tester <- lapply(test,unlist)
  # #another for loop to apply reverse complementary to antisense alignments 
  # for( ia in 1:dim(WorkingD)[1]){
  #   
  #   #For every sequence, check if WorkingD$Sense == 1, if yes: do reverse complementary on sequence in test[[ia]] 
  #   if(identical(WorkingD$Sense[ia],1)){
  #     
  #     test[[ia]] <- rev(comp(test[[ia]],forceToLower = FALSE))
  #     
  #   }
  #   
  # }
  # 
  # # tester[[2]] <- rev(comp(tester[[2]]  ))
  # #Stitch sequences together in pos sense:
  # 
  # 
  # # !WorkingD$Sense
  # 
  # #Determine order of WorkingD rows by Sstart 
  # test <- test[c(order(WorkingD$Sstart))]
  # #rev(comp(seq))
  # 
  # 
  # 
  # 
  # 
  # Finaldata[[ Seqnms[i] ]] <- (test)
  Bleh[[i]] <- list(WorkingD, c( levels(as.factor(WorkingD$Subjacc)) ,paste0(rownames(WorkingD),collapse="."), as.character(sum(WorkingD$Bitscore)) ), c(dim(Data)[1],dim(WorkingD)[1]) )
  # write.fasta(help, file.out = paste0(dOutput,"elp"), names = c())
  
  
  
}
# rm(i);rm(ij);rm(j);rm(z);rm(a);rm(b);rm(BestP);rm(test);rm(WorkingD);rm(TempVar1);rm(OLMat);rm(Data);rm(Data2);rm(OLcheck)
# rm(nOLroute);rm(OLMatsum);rm(ja);rm(jb);rm(TempVar2);rm(TempVar3);rm(TempVar4)
# rm(abb);rm(TempVar5);rm(TempVar6);rm(TempVar7);rm(ia);rm(jbz);rm(jc);rm(jcz);rm(ji);rm(jl);rm(TempTemp5)

# 
# Finaldata2 <- lapply(Finaldata, unlist)
# 
# for(i in 1:length(Finaldata2)){
#   Finaldata2[[i]] <- paste0(Finaldata2[[i]],collapse="")
#   
# }
# 
# 
Bleh

}
# 
# 
# write.fasta(Finaldata2, file.out = paste0(dOutput, basename(dOutput), "_elp.fasta"),names=Seqnms, as.string = TRUE)
