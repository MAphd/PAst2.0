#### PAst 
#the Pseudomonas aeruginosa serotyper (PAst) v2

#R conversion

PAst <- function(Output, Input, OSAdb, Blasted){
  
  #todo: Check if all inputs contain directories, and they do not contain backslashes
  
  # source("./WeightedAlignmentRouter.R")

  library(readr)
  #accepted extension types .fasta, .fas, .fa, .seq, .fsa 
  #todo: Record a variable for the extension since we use it for naming later on and right now, just assume .fasta

  #todo: only read .fasta files
  Seqnms <- list.files(Input) #,pattern = c("*fasta$","*fna$","*fsa$","*fas$","*fa$","*seq$"))
  
  
  #Check if files have valid extension names 
  fileval <- lapply(c("*.fasta$","*.fna$","*fsa$","*fas$","*fa$","*seq$"), grepl, list.files(Input))
  
  temphere <- data.frame(matrix(ncol=length(fileval),nrow=length(list.files(Input)) ))
  for( i in 1:length(fileval) ){
    temphere[,i] <- fileval[[i]]
  }
  
  #Remove files that aren't allowed 
  Seqnms <- Seqnms[ apply(temphere,1,any) ] 
  

  ###### Check for existing blast results if flag is set (one more input) and use these instead of blasting again
  
  if(length(Seqnms)>=1){ 
  print( paste("Succesfully loaded ",length(Seqnms)," sequences", sep=""))
  } else {
    
    stop("No input files could be loaded")
  }
  
  if( !Blasted ){ 
  
  Cmdline <- function(Seqfile){
    paste("blastn -query ", Input, Seqfile, " -subject ", OSAdb, " -out ",
          Output, "PAst_", Seqfile,
          ".txt -outfmt ", "\"", "6 qaccver saccver pident sstart send length evalue bitscore nident", "\"" , sep="")

  }


  Cmds <- mapply(FUN = Cmdline, Seqfile = Seqnms) 
  

  sapply(Cmds, system)  
  
  }
  
  #Blastn spews an error if there are any spaces "Error:  (CArgException::eSynopsis) Too many positional arguments (1)"
  #in the path or file name.
  

  
  #Format: (default): qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore
  #Format (used): qaccver saccver pident sstart send length evalue bitscore nident"
  
  #For each sequence file, we load the corresponding output
  
  Blastnms <- paste("PAst_",Seqnms,".txt",sep="")
  # Blastnms <- c("PAst_GCA_000795365.1_AZPAE14827_genomic.fna.txt","PAst_GCA_001450005.1_WH-SGI-V-07237_genomic.fna.txt")
  #paste(Output, Blastnms,sep="")
  
  #Serogroup reference sizes: O1-O13, and WzyB, O15, O17
  Refsize <- data.frame(matrix(nrow=14))
  Refsize <- Refsize[,-1]
  Refsize$Serogroup <- c("O1","O2","O3","O4","O6","O7","O9","O10","O11","O12","O13","WzyB","O15","O17")
  #Define serotype cluster sizes, according to order given above for the purpose of calculating coverage
  Refsize$Size <- c(18368, 23303, 20210, 15279, 15649, 19617, 17263, 17635, 13868, 25864, 14316, 1140, 15676, 22236)
  
  Finaldata <- data.frame()
  for(i in 1:length(Seqnms)){
 # i <- 15
    Data <- data.frame(read_delim(paste(Output, Blastnms,sep="")[i], delim="\t", col_names = c("Queryacc","Subjacc","Identity","Sstart","Send","Length","Evalue","Bitscore", "nident")))
    #assume false to reset after each loop
    WzyB <- FALSE

    
    #First we check for presence of WzyB
    if(any(grepl("WzyB",Data$Subjacc))){
      
      #Check for multiple fragments in different contigs 
      
      #If present, we check the coverage of WzyB
      #WzyB returns true if our query covers over 95% of the subject WzyB
      WzyB <- 100*(1-(Refsize$Size[12]-as.numeric(Data$Length[grep("WzyB", Data$Subjacc)]) )/Refsize$Size[12]) >95
      # stop(paste(i))
    } 
    
    #Remove row with WzyB if present
    if(WzyB){
     Data <- Data[-grep("WzyB",Data$Subjacc),]

    }
    
    #Calculate best routes for each subjacc
    Data2 <- WeightedAlignementRouter(Data)
    
    TempData1 <- data.frame(matrix(ncol=6,nrow=length(Data2)))
    for( j in 1:length(Data2)){
      
      TempData1[j,] <- c( Data2[[j]][[2]][c(1,2)] ,dim(Data2[[j]][[1]])[1], colSums(Data2[[j]][[1]][c(6,8,9)]) ) 
      
    }
    colnames(TempData1) <- c("Subjacc","Row","Fragments","Length","Bitscore","nident")
    
    
    Data <- TempData1
    
    Data$Cover <- NA
    #Create vector of coverage for each result
    for(z in 1:dim(Data)[1]){
    
      Data$Cover[z] <- 100*(1-((Refsize$Size[which(Refsize$Serogroup%in%Data$Subjacc[z])] - as.numeric(Data$Length[z]))/ (Refsize$Size[which(Refsize$Serogroup%in%Data$Subjacc[z])])))
      
      # Data$Cover[z] <- (100-100* (Refsize$Size[which(Refsize$Serogroup%in%Data$Subjacc[z])]-as.numeric(Data$Length[z])) /Refsize$Size[which(Refsize$Serogroup%in%Data$Subjacc[z] )])
      if(Data$Cover[z] > 100){
        Data$Cover[z] <- 100
      }
    
    }
    
    Data$Identity <- (as.numeric(Data$nident)/as.numeric(Data$Length))*100
    Data$Bitscore <- as.numeric(Data$Bitscore)
    #Type <- NA
    #Serotyping starts here:
   
    #Does any of our hits have >95% coverage?
    if( any(Data$Cover>95) ){
      #Does multiple results have >95% coverage?
      if(sum(Data$Cover>95)>=2 ){
        #Multitypeable?
        Type <- "Multiple"
        } else {
          #Is WzyB present at >95% coverage in our BLAST hits?
          if(WzyB){
            #Is the query with >95% coverage serogroup O2?
            if(identical("O2",Data$Subjacc[which(Data$Cover>95)])){
              Type <- "O2"
              } 
            } else { #No wzyB
              #Is the query with >95% coverage serogroup O2?
              if(identical("O2",Data$Subjacc[which(Data$Cover>95)])){
                Type <- "O5"
              } else {
                #Else we base it on the best result (highest coverage) that has a cover > 95%
                Type <- Data$Subjacc[Data$Cover==max(Data$Cover) ]
              }
              
            }
        }
      } else { 
        # If , The length of the alignment for the result with the largest coverage is over 4000, type based on result with largest coverage
        # if( max(as.numeric(Data$Length[Data$Cover == max(Data$Cover)])) > 4000 ){
        #   Type <- Data$Subjacc[Data$Cover==max(Data$Cover)  ]
        # } else { 
        Type <- "NT"
        # }
      }
    
    #Done typing
    # 
    # if( Type %in% c("NT","Multiple") ){
    #   
    #   Finaldata <- rbind(Finaldata,
    #                      cbind(paste(unlist(strsplit(Seqnms,".fasta")),sep="")[i],
    #                            Type,
    #                            0,
    #                            0,
    #                            0,
    #                            0,
    #                            0,
    #                            0))
    #   
    #   
    # } else {
    #   
     
    #Special if clause for multi typable to prevent it from creating multiple rows per result
    if( identical(Type,"Multiple")  ){
      
      #Chose first element in the unlikely event that the top alignments share identical bitscores
      Finaldata <- rbind(Finaldata,
            cbind(paste(unlist(strsplit(Seqnms,".fasta")),sep="")[i],
                  Type,
                  Data$Cover[ Data$Bitscore==max(Data$Bitscore)  ][1],
                  Data$Identity[ Data$Bitscore==max(Data$Bitscore) ][1],
                  Data$Queryacc[ Data$Bitscore==max(Data$Bitscore) ][1],
                  Data$Row[ Data$Bitscore==max(Data$Bitscore) ][1],
                  Data$Subjacc[ Data$Bitscore==max(Data$Bitscore) ][1],
                  Data$Fragments[ Data$Bitscore==max(Data$Bitscore) ][1]
            ))
      
    } else { 
      #Only works for the case where the Type is the same as the subjacc (which wont be the case for O5 )
      #workaround:
    if( Data$Subjacc[ Data$Cover==max(Data$Cover) ] == "O2" & !WzyB ){
      Data$Subjacc[ Data$Cover==max(Data$Cover) ] <- "O5"
    }
    
    
      Finaldata <- rbind(Finaldata,
                         cbind(paste(unlist(strsplit(Seqnms,".fasta")),sep="")[i],
                               Type,
                               Data$Cover[ Data$Cover==max(Data$Cover)  ],
                               Data$Identity[ Data$Cover==max(Data$Cover) ],
                               Data$Queryacc[ Data$Cover==max(Data$Cover) ],
                               Data$Row[ Data$Cover==max(Data$Cover) ],
                               Data$Subjacc[ Data$Cover==max(Data$Cover) ],
                               Data$Fragments[ Data$Cover==max(Data$Cover) ]
                                              ))
      
    }
      
      

    
    
    
    
    
  }
  colnames(Finaldata) <- c("Name","Serotype","Coverage","Identity","Row","Besthit","Fragments")
  Finaldata$Identity <- round(as.numeric(as.vector(Finaldata$Identity)), digits= 2 )
  Finaldata$Coverage <- round(as.numeric(as.vector(Finaldata$Coverage)), digits= 2 )
  
  write.table(Finaldata, file=paste(Output,gsub("/","",strsplit(Input,dirname(Input))[[1]][2]),"_Serotyping.txt",sep=""),row.names=FALSE)
       
   
  # print(Finaldata)
  print(paste("Report saved to ",Output,gsub("/","",strsplit(Input,dirname(Input))[[1]][2]),"_Serotyping.txt",sep=""))
  
  
  Finaldata 
 }
  
  ##
  # OSA clusters with >95% coverage in the query genome represent a positive hit
  # Multiple results at >95% represent multi-typeable
  # No hits at >95% represents non-typeable
  #
  # Serogroups o1-o13 can be defined
  # wzybeta used to distinguish serogroup O2(O2 and O16) and O5(o5, O18, O20)
  # Check for coverage of each

