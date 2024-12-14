#!/usr/bin/env Rscript
#perl TargetPerla.pl SourceTarget.txt S04380110_Covered.bed MyTarget_w10k_231024 10000 hg38 
excavator1.f <- function(path.case, path.out,batch.case){
  
  dir.create(path.out)
  path.cont = "[path_to_control_bam_files]"
  setwd(path.cont) 
  bam.list.cont <- dir(pattern=".bam$") #
  bam.list.cont = bam.list.cont[1:10]  
  
  
  setwd(path.case) 
  bam.list.case <- dir(pattern=".bam$")
  
  bam.l = length(bam.list.cont)+length(bam.list.case)
  
  setwd("[path_to_EXCAVATOR2_package]/EXCAVATOR2_Package_v1.1.2")
  
  # case
  k=1
  for(i in 1:length(bam.list.case)){
    SampleId = unlist(strsplit(bam.list.case[i],split=".bam",fixed=T))[1]
    write.table(paste(path.case, "/",unlist(bam.list.case)[i]," ", path.out, "/",SampleId, " ", SampleId, sep=""),
                paste(path.out, "/", batch.case, "_ExperimentalFilePrepare_w20k_", length(bam.list.case), "_", length(bam.list.cont), ".txt", sep=""), col.names=F, row.names=F, quote=F, append=T)
    k = k + 1
  }
  
  
  
  # control
  k=1
  for(i in 1:length(bam.list.cont)){
    ContrId = unlist(strsplit(bam.list.cont[i],split=".bam",fixed=C))[1]
    write.table(paste(path.cont, "/",unlist(bam.list.cont)[i], " ", path.out, "/",ContrId, " ", ContrId, sep=""), 
                paste(path.out, "/", batch.case, "_ExperimentalFilePrepare_w20k_", length(bam.list.case), "_", length(bam.list.cont), ".txt", sep=""), col.names=F, row.names=F, quote=F, append=T)
 	   k = k + 1
  }
  
  # Running ReadPerla.pl
  ## --mode somatic: case & control= 1:1 matching
  ## --mode pooling
  
  system(paste("perl", " ", "EXCAVATORDataPrepare.pl", " ", path.out,"/", batch.case, "_ExperimentalFilePrepare_w20k_", length(bam.list.case), "_", 
               length(bam.list.cont), ".txt ", " ", "--processors 6 ", "--target MyTarget_w20k_240404", " ", "--assembly hg38",sep=""))
}


# patient_wes: 9 samples
excavator1.f(path.case = "[bam_file_dir]", path.out = "[output_dir]", batch.case = "[output_name]")

####################################################################################################################


excavator2.f <- function(path.case, path.out,batch.case){
  
  dir.create(path.out)
  path.cont = "[path_to_control_bam_files]"
  setwd(path.cont) 
  bam.list.cont <- dir(pattern=".bam$") #
  bam.list.cont = bam.list.cont[1:10]  
  
  
  setwd(path.case) 
  bam.list.case <- dir(pattern=".bam$")
  bam.l = length(bam.list.cont)+length(bam.list.case)
  
  setwd("[path_to_EXCAVATOR2_package]/EXCAVATOR2_Package_v1.1.2")
  
  # case
  k=1
  for(i in 1:length(bam.list.case)){
    SampleId = unlist(strsplit(bam.list.case[i],split=".bam",fixed=T))[1]
    print(bam.list.case[i])
    write.table(paste("T", k, " ", path.out, "/",SampleId, " ", SampleId, sep=""),
                paste(path.out, "/", batch.case, "_ExperimentalFilePrepare_w20k_v2_", length(bam.list.case), "_", length(bam.list.cont), ".txt", sep=""), col.names=F, row.names=F, quote=F, append=T)
    k = k + 1
  }
  
  
  
  # control
  k=1
  for(i in 1:length(bam.list.cont)){
    ContrId = unlist(strsplit(bam.list.cont[i],split=".bam",fixed=C))[1]
    write.table(paste("C", k, " ",  path.out, "/",ContrId, " ", ContrId, sep=""), 
                paste(path.out, "/", batch.case, "_ExperimentalFilePrepare_w20k_v2_", length(bam.list.case), "_", length(bam.list.cont), ".txt", sep=""), col.names=F, row.names=F, quote=F, append=T)
    k = k + 1
  }
  
  # Running ReadPerla.pl
  ## --mode somatic: case & control= 1:1 matching
  ## --mode pooling
  

  system(paste("perl", " ", "EXCAVATORDataAnalysis.pl"," ", path.out, "/", batch.case, "_ExperimentalFilePrepare_w20k_v2_", length(bam.list.case), "_", length(bam.list.cont), ".txt ", "--processors 6  ", " ", "--target MyTarget_w20k_240404", " ", "--assembly hg38", " ", "--output", " ", path.out, "/Excavator_w20K", " ", " --mode pooling", sep=""))
}

# patient_wes
excavator2.f(path.case = "[bam_file_dir]", path.out = "[output_dir]",  batch.case = "[output_name]")
