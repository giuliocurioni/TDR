###tdr.cal.bec.R <- function(...) {
# Giulio Curioni
# 13 October 2017 

## perform two steps calibration (Huisman et al., 2008, Bechtold et al., 2010) for conductivity
 # note: the script works only if the file names contain information on the solutions used (e.g. sol1, sol2 etc. and their conductivity, e.g. 44.4 (mS/m) etc.), see below, sample details should be e.g. sol4-22.5
 
      # check ### before running the script
 
     # waveform should be saved with the following format, keep '-' as the separator EXCEPT temperature (if present). Use '_' to contain temperature, if temperature is measured. See examples below for conductivity measurements used for calibration. The file formats can be different, depending on which software was used to take measurements (e.g. '.mat', '.DAT'):
     # i.e. probe.ID  -  multiplexer level  -  start plot +m  -  window length +m  -  sample details  _  temperature +C  _  repetition
     # probe ID (i.e. probe name, as present in the file 'tdr.probe.par.txt')
     # mux ID (i.e. mux0 = no multiplexers, mux2 = two multiplexers etc.)
     # start parameter +m (m), i.e. where the plot starts  (see example below)
     # window length parameter +m (m), i.e. where the plot ends  (see example below)
     # sample ID + additional info separated by '-'
     # if present, temperature with 1 decimal point +C, contained within two '_'    (see example below)
     # repetition information  (always add an entry if temperature was taken, see examples)
     
     # examples:
     # e.g.1 lab1.2-mux0-7.6m-2.6m-SC.mix1-16gwc_23.1C_rep1
     # e.g.2 p15.1-mux2-0m-500m-water_20.0C_1  # IMPORTANT: ALWAYS add something after temperature (e.g. use _1, _rep1 or _r1 etc. for no repetitions, i.e. one measurement only. Note: There is no need to add anything if temperature was not taken)

     # additional examples specific to conductivity measurements used for calibration:
     # lab1.3-mux0-0m-500m-sol2-1124_16.8C_rep4.DAT     # note that 'sol1', 'sol2' etc. and the ref. cond in mS/m MUST be present
     # lab1.5-mux2-0m-500m-air.open_21.0C_rep3.DAT      # always use 'air.open' for open measurements in air

library(dplyr)
library(readr)

loadall()
graphics.off()
options(digits=8)
filter <- dplyr::filter    

### first run tdr.cond.R to obtain the reflection coefficient at long distances 
# NOTE1: It is possible to run tdr.cond.R on ALL the measurements taken with different probes, including the open measurements in air.
# NOTE2: if calibrating multiple probes with different TDR units, run the script separately for each TDR unit used (tdr.name needs to be fixed, see below).
 
### make sure to add an entry (e.g. row) in tdr.probe.par.txt corresponding to the probes analysed. This is necessary to obtain L.cable.


# import analysed measurements with tdr.cond.R (only RC is needed, not BEC values so the probe name only needs to match the one used in tdr.cond.R)
d0 <- read.delim("lab1.4.BEC.txt", header=T)     ### change file name if necessary

       tdr.name <- "ATU.red"     ### choose the TDR unit used for the measurements. Allowed options are: TDR19200, Dan.PhD, Giulio.PhD, ATU.lab, ATU.red, ATU.blue              # DO NOT analyse probes calibrated using different TDR units. Re-run the script separately for each TDR unit.
       
       #L.cable <- 6.2    # input cable length (m) corresponding to the probe used
                  


additional.headers <- which(d0[,1] == "file.name")
if(length(additional.headers) != 0) { d0 <- d0[-additional.headers,] }
d0 <- data.frame(lapply(d0[,1:length(d0)], as.character), stringsAsFactors=FALSE)
file.name <- d0$file.name
d0 <- data.frame(lapply(d0[,2:length(d0)], as.numeric), stringsAsFactors=FALSE)
d0 <- cbind(file.name, d0)
d0$file.name <- as.character(d0$file.name)
d0$BEC <- d0$BEC.mS.m / 1000

d0.file.name <- d0$file.name

### were the files taken with the Matlab script? if so set matlab.format = TRUE, otherwise set it to FALSE

matlab.format <- FALSE

if(matlab.format == T) {d0.file.name <- paste(d0.file.name, ".mat", sep='')}


# find temperature, repetition etc.
ee <- data.frame(d0.file.name)
ee$d0.file.name <- as.character(ee$d0.file.name)

# general split
split1 <- separate(ee, col=d0.file.name, into=c("probe", "mux", "start.plot", "win.Lapp", "other"), sep= "-")    ### this assumes that the file names are saved in a constant format. Modify "into" if necessary.


# probe
probe <- split1$probe
mux <- split1$mux
probe_mux <- paste(probe, mux, sep="_")

# find temperature
temp <- split1 %>%
        select(other) %>%
        separate(col=other, into=c("medium", "temp", "rep"), sep= "_") %>%
        select(temp)        

# find repetition
repetition <- split1 %>%
              select(other) %>%
              separate(col=other, into=c("medium", "temp", "rep"), sep= "_") %>%
              select(rep) %>%
              separate(col=rep, into=c("rep", "rep2"), sep= "[.]")        

# find medium
medium <- split1 %>%
              select(other) %>%
              separate(col=other, into=c("medium", "temp", "rep"), sep= "_") %>%
              select(medium)


#              Splits <- unlist(strsplit(d0.file.name[1], "-"))  # it assumes that each waveform is named in a standard way
#              split.name <- matrix(NA, nrow=nrow(d0), ncol=length(Splits) + 1)
#              split.name <- data.frame(split.name)
#              for (ll in 1:nrow(d0)) {
#                  #split.name[ll,] <- unlist(strsplit(file.name.col[ll], "[.]"))
#                  ll.split <- unlist(strsplit(d0.file.name[ll], "-"))
#                  split.name[ll,1:length(ll.split)] <- ll.split
##                  if("metal" %in% ww[ll,] == TRUE) {mould.type[ll,1] <- "metal"}
##                  if("plastic" %in% ww[ll,] == TRUE) {mould.type[ll,1] <- "plastic"}
#                  }
#              colnames(split.name) <- c("probe", "mux", "start.plot", "win.Lapp", "medium", "temp.rep") # change if necessary
#              
###              
#              # correct issue with air.open measurements (because theere is no ref conductivity information the column with "medium" is not correct):
#              temp.rep.NA <- which(is.na(split.name$temp.rep))
#              split.name$temp.rep[temp.rep.NA] <- split.name$medium[temp.rep.NA]
#              split.name$medium[temp.rep.NA] <- "open.air"       ### change if necessary (sometimes open.air is used...)
              
              
#              # select temperature
#              Splits <- unlist(strsplit(d0.file.name[1], "_"))  # it assumes that each waveform is named in a standard way
#              temp <- matrix(nrow=nrow(d0), ncol=length(Splits) + 1)
#              temp <- data.frame(temp)
#              for (ll in 1:nrow(d0)) {
#                  #split.name[ll,] <- unlist(strsplit(file.name.col[ll], "[.]"))
#                  ll.split <- unlist(strsplit(d0.file.name[ll], "_"))
#                  temp[ll,1:length(ll.split)] <- ll.split
##                  if("metal" %in% ww[ll,] == TRUE) {mould.type[ll,1] <- "metal"}
##                  if("plastic" %in% ww[ll,] == TRUE) {mould.type[ll,1] <- "plastic"}
#                  }
#              if(ncol(temp) > 1) {    
#              temp <- temp[,2]  ## ## important: this assumes that temperature is contained within '_' and that the filename continues after that (i.e. has a repetition value)
#              # temperature must be specified with 1 decimal value
#              # select numeric value for temperature
#              for (ll in 1:nrow(d0)) {
#                  temp[ll] <- substring(temp[ll], 1, (nchar(temp[ll])-1))         # this assumes that temperature is specified with 1 decimal value
#                  }
#              temp <- as.numeric(temp)
#              } else {temp <- NA}           
              
              
#              # select repetition
#              Splits <- unlist(strsplit(d0.file.name[1], "_"))  # it assumes that each waveform is named in a standard way
#              split0 <- matrix(nrow=nrow(d0), ncol=length(Splits) + 1)
#              split0 <- data.frame(split0)
#              for (ll in 1:nrow(d0)) {
#                  #split.name[ll,] <- unlist(strsplit(file.name.col[ll], "[.]"))
#                  ll.split <- unlist(strsplit(d0.file.name[ll], "_"))
#                  split0[ll,1:length(ll.split)] <- ll.split
##                  if("metal" %in% ww[ll,] == TRUE) {mould.type[ll,1] <- "metal"}
##                  if("plastic" %in% ww[ll,] == TRUE) {mould.type[ll,1] <- "plastic"}
#                  }
#              Rep <- split0[, length(Splits)]    
#              Splits <- unlist(strsplit(Rep[1], "[.]"))  # it assumes that each waveform is named in a standard way
#              repetition <- matrix(nrow=nrow(d0), ncol=length(Splits) + 1)
#              repetition <- data.frame(repetition)
#              for (ll in 1:nrow(d0)) {
#                  #split.name[ll,] <- unlist(strsplit(file.name.col[ll], "[.]"))
#                  ll.split <- unlist(strsplit(Rep[ll], "[.]"))
#                  repetition[ll,1:length(ll.split)] <- ll.split                  
##                  if("metal" %in% ww[ll,] == TRUE) {mould.type[ll,1] <- "metal"}
##                  if("plastic" %in% ww[ll,] == TRUE) {mould.type[ll,1] <- "plastic"}
#                  }
#              repetition <- repetition[,1]
                  
              
#              # select reference conductivity from conductivity meter
#              Splits <- unlist(strsplit(split.name$temp.rep[1], "_"))  # it assumes that each waveform is named in a standard way
#              split0 <- matrix(nrow=nrow(d0), ncol=length(Splits) + 1)
#              split0 <- data.frame(split0)
#              for (ll in 1:nrow(d0)) {
#                  #split.name[ll,] <- unlist(strsplit(file.name.col[ll], "[.]"))
#                  ll.split <- unlist(strsplit(split.name$temp.rep[ll], "_"))
#                  split0[ll,1:length(ll.split)] <- ll.split
##                  if("metal" %in% ww[ll,] == TRUE) {mould.type[ll,1] <- "metal"}
##                  if("plastic" %in% ww[ll,] == TRUE) {mould.type[ll,1] <- "plastic"}
#                  }
#              ref.cond <- split0[, 1]
#              ref.cond.NA <- which(ref.cond == "open.air") ### change if necessary (e.g. sometimes open.air is used...
#              ref.cond[ref.cond.NA] <- NA
              
              
              ### ### IF THE REF COND VALUES WERE NOT SAVED IN THE FILE NAME IMPORT THEM SEPARATELY (comment this section if the ref cond is saved in the file names):
              path_ref_cond <- data.frame(choose.files(caption = "Select reference conductivity TXT file"))
              ref_cond <- read_delim(as.character(path_ref_cond[1,1]), delim="\t", col_names = TRUE)
              ref.cond <- medium$medium
              ref.cond[which(ref.cond == "air.open")] <- NA
              ref.cond[which(ref.cond == "open.air")] <- NA 
              ref.cond[which(ref.cond == "sol1")] <- ref_cond$conductivity_S_m[which(ref_cond$solution == "sol1")]  
              ref.cond[which(ref.cond == "sol2")] <- ref_cond$conductivity_S_m[which(ref_cond$solution == "sol2")]
              ref.cond[which(ref.cond == "sol3")] <- ref_cond$conductivity_S_m[which(ref_cond$solution == "sol3")]
              ref.cond[which(ref.cond == "sol4")] <- ref_cond$conductivity_S_m[which(ref_cond$solution == "sol4")]
              ref.cond[which(ref.cond == "sol5")] <- ref_cond$conductivity_S_m[which(ref_cond$solution == "sol5")]
              ref.cond[which(ref.cond == "sol6")] <- ref_cond$conductivity_S_m[which(ref_cond$solution == "sol6")]
              ref.cond[which(ref.cond == "sol7")] <- ref_cond$conductivity_S_m[which(ref_cond$solution == "sol7")]
              ref.cond[which(ref.cond == "sol8")] <- ref_cond$conductivity_S_m[which(ref_cond$solution == "sol8")]
              
              
              ref.cond <- as.numeric(ref.cond)                   
              #ref.cond <- ref.cond/1000       # convert mS/m to S/m 
                  
data_BEC <- cbind(d0, probe, mux, probe_mux, medium, temp, repetition, ref.cond)                  



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 



## summarise results
                 
               
by_groups <- group_by(data_BEC, probe, mux, medium, add=TRUE)
summary_RC <- summarise(by_groups, 
               mean.RC = round(mean(RC, na.rm = TRUE), digits=8),
               min.RC = round(min(RC, na.rm = TRUE), digits=8),
               max.RC = round(max(RC, na.rm = TRUE), digits=8),
               max_min.RC = round(max(RC, na.rm = TRUE) - min(RC, na.rm = TRUE), digits=8), 
               sd.RC = round(sd(RC, na.rm = TRUE), digits=8)
               )
            
#ref.cond <- dplyr::filter(data_BEC, probe == "lab1.3" & mux == "mux0")

#by_groups <- group_by(data_BEC, probe, mux, medium, add=TRUE)
summary_ref_cond <- summarise(by_groups, 
               mean.ref.cond = round(mean(ref.cond, na.rm = TRUE), digits=5),
               min.ref.cond = round(min(ref.cond, na.rm = TRUE), digits=5),
               max.ref.cond = round(max(ref.cond, na.rm = TRUE), digits=5),
               max_min.ref.cond = round(max(ref.cond, na.rm = TRUE) - min(ref.cond, na.rm = TRUE), digits=5), 
               sd.ref.cond = round(sd(ref.cond, na.rm = TRUE), digits=5)
               )
               
data.cal.tdr <- select(summary_RC, medium, mean.RC) 
data.cal.tdr <- arrange(data.cal.tdr, probe, medium)
data.cal.hannah <- select(summary_ref_cond, medium, mean.ref.cond)
data.cal.hannah <- arrange(data.cal.hannah, probe, medium)

data.cal <- bind_cols(data.cal.tdr, data.cal.hannah)
data.cal <- select(data.cal, -probe1, -mux1, -medium1)
data.cal$probe.mux <- paste(data.cal$probe, data.cal$mux, sep='_')

list.cal <- split(data.cal, c(data.cal$probe.mux))
 
    

 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
                        

## * ##    
for(ii in 1:length(list.cal)) { 

       aa <- data.frame(list.cal[ii])
       colnames(aa) <- c("probe", "mux", "medium", "RC", "ref.cond", "probe.mux")
       
       
       
       probe.name <- aa$probe[1]  ## ## choose probe (it must match one of the existing probe.ID in tdr.probe.par.txt)
       mux.level <- aa$mux[1]     ## ## choose the multiplexer level used for the measurements. Allowed options are: mux0, (mux1), mux2   




## import probe parameters:
    tdr.probe.par <- read_delim("C:\\R.wd\\scripts\\tdr\\tdr.probe.par.txt", delim='\t', col_names=T)
    #tdr.probe.par <- format(tdr.probe.par, digits=8)
    tdr.unit <- tdr.probe.par$tdr.unit
#    probe.ID <- tdr.probe.par$probe.ID
#    mux <- tdr.probe.par$mux
#    tdr.probe.par <- data.frame(lapply(tdr.probe.par[,4:length(tdr.probe.par)], as.character), stringsAsFactors=FALSE)
#    tdr.probe.par <- data.frame(lapply(tdr.probe.par[,1:length(tdr.probe.par)], as.numeric), stringsAsFactors=FALSE)
#    tdr.probe.par <- cbind(tdr.unit, probe.ID, mux, tdr.probe.par)
#    tdr.probe.par$tdr.unit <- as.character(tdr.probe.par$tdr.unit)
#    tdr.probe.par$probe.ID <- as.character(tdr.probe.par$probe.ID)
#    tdr.probe.par$mux <- as.character(tdr.probe.par$mux)  
#
## define probe parameters corresponding to the chosen probe (see probe.name above)
    n.probe.ID.diff.units <- which(tdr.probe.par$probe.ID == probe.name)
    n.probe.ID <- which(tdr.probe.par$tdr.unit == tdr.name & tdr.probe.par$probe.ID == probe.name & tdr.probe.par$mux == mux.level)   # choose the probe identified by probe.name and tdr.name
        
    if(length(n.probe.ID) == 0) {stop("Double check that you chose the correct tdr.name in the script and make sure that there is an entry in tdr.probe.par.txt corresponding to the probes analysed.")}
    
    L.cable <- tdr.probe.par$L.cable[n.probe.ID]             # cable length (m)  
    n0 <- tdr.probe.par$n.data.points[n.probe.ID]            # number of data points
    start.plot.bec <- tdr.probe.par$start.plot.bec[n.probe.ID]  # apparent length (m) corresponding to the parameter start in PCTDR (start of the plot) for conductivity measurements
    win.Lapp.bec <- tdr.probe.par$win.Lapp.bec[n.probe.ID]      # apparent length (m) of the window (entire TDR plot = length parameter in PCTDR) for conductivity measurements
    long.dist.points <- tdr.probe.par$long.dist.points[n.probe.ID]  # number of long distance data points to be used for conductivity analysis


       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       # correct for open measurement (Lin et al., 2007, 2008)          
       Ropen <- aa$RC[which(aa$medium == "open.air")]   ### change if necessary. Open measurement (in theory it should be 1, but it never is in practice, so scale the measurements so that an open would be 1)
       aa <- filter(aa, !is.na(ref.cond))      ### change if necessary to air.open
       aa$RC.corr <- ((2*(aa$RC+1))/(Ropen+1))-1        # corrected reflection coefficient using open measurements in air
       
       # calculate load resistance and sample conductance
       aa$RL <- 50*((1+aa$RC.corr)/(1-aa$RC.corr))   # load resistance
       aa$Gs <- 1/aa$RL      # sample conductance
       
       # find probe constant (slope of ref.cond vs sample conductance)
       index.Kp <- which(aa$medium == "sol5" | aa$medium == "sol6" | aa$medium == "sol7" | aa$medium == "sol8")
       Gs.low <- aa$Gs[index.Kp]
       lm.Kp <- lm(Gs.low ~ aa$ref.cond[index.Kp])
       Kp <- 1/coef(lm.Kp)[2]
       
       probe.name <- aa$probe.mux[1]
       
#       write.table(Kp, paste(probe.name, '_', "Kp.txt", sep=''), row.names=F, col.names="Kp", sep='\t')
 
    

 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
                        

    
## simplex nedler and mead, optimisation to find Rc and R0 (extra series resistance parameters)

      Cond <- function(r){
            # input ref cond (e.g. HI8733)
            y <<- aa$ref.cond
            
            # input RL 
            RL <<- aa$RL
            
            Rc <- r[1]
            R0 <- r[2]
            zmin <- (sum((y - (Kp / (RL - (L.cable*Rc+R0))) )^2)) ^0.5	#I need to minimise this quantity by optimising Rc and Ro
            }

      #find the best fit combining a number of possible Rc and R0 values (from 0 to 0.5, with steps 0.01)
      seqR <- seq(0,0.5,0.01)
      seqrc <- rep(seqR, each=length(seqR))	#repeats a number of times the elements in seqR!!!
      seqr0 <- rep(seqR, length(seqR))
      Niter <- length(seqR)*length(seqR)
      bah <- matrix(nrow=Niter, ncol=5)
      bah <- data.frame(bah)
      for (i in 1:Niter){ 
      	prova <- optim(c(seqrc[i],seqr0[i]), Cond)
      	bah[i,1] <- seqrc[i]	
      	bah[i,2] <- seqr0[i]
      	bah[i,3] <- prova$value[1]
      	bah[i,4] <- prova$par[1]
      	bah[i,5] <- prova$par[2]
      	colnames(bah) <- c("Rc.ini", "R0.ini", "zmin", "Rc", "R0")
        }
      	
      	
      write_delim(bah, paste(probe.name, "_best.optim.txt", sep=""), delim='\t', append=F)
      minbah <- which.min(bah[,3])
      #print("global minimum:", quote=F)
      #print(bah[minbah,], quote=F)
      
      Rc1 <- bah[minbah,4] 
      R01 <- bah[minbah,5]
      
      ##in case something is wrong (e.g. Rc values are too high) use this:
      bah1 <- subset(bah, bah[,4]<0.5)     ### change <0.3 value if there are no good results
      minbah1 <- which.min(bah1[,3])
      #print("NOT USED IN PLOTS: minimum with Rc < 0.1, see Bechtold et al., 2010")
      #print(bah1[minbah1,], quote=F)
      
      #plot simplex model results
      par(mar=c(5, 4.5, 4, 2) + 0.1)
      plot(RL,y, lwd=2, cex.axis=1.3, cex.lab=1.3, font=1, font.lab=1, xlab='Load resistance [ohm]', ylim=c(0,2), ylab='Electrical conductivity [S/m]')	#main='Modelled vs measured'	
      curve((Kp / (x - (L.cable*Rc1+R01))), lwd=2, from=min(RL)-1, to=max(RL), add=TRUE)	#if curve is not complete change N. in min(RL)- N.)
      par(font=1)
      box(lwd=1)
      savePlot(filename = paste(probe.name, ".simplex_model.png", sep=""), type = c("png"))    # save plot	
      
      
      windows()	#opens a new figure
      
      z1 <- Kp/(RL[1]-(L.cable*Rc1+R01))
      z2 <- Kp/(RL[2]-(L.cable*Rc1+R01))
      z3 <- Kp/(RL[3]-(L.cable*Rc1+R01))
      z4 <- Kp/(RL[4]-(L.cable*Rc1+R01))
      z5 <- Kp/(RL[5]-(L.cable*Rc1+R01))
      z6 <- Kp/(RL[6]-(L.cable*Rc1+R01))
      z7 <- Kp/(RL[7]-(L.cable*Rc1+R01))
      z8 <- Kp/(RL[8]-(L.cable*Rc1+R01))
      
      z <- c(z1,z2,z3,z4,z5,z6,z7,z8)
      
      par(mar=c(5, 4.5, 4, 2) + 0.1)
      plot(y,z, log='xy', lwd=2, cex.axis=1.3, cex.lab=1.3, font=1, font.lab=1, xlab='Reference conductivity [S/m]', ylab='Modelled conductivity [S/m]',xaxt="n", yaxt="n")
      ##curve((Kp / (x - (L.cable*Rc1+R01))), lwd=2, from=min(RL)-5, to=max(RL), add=TRUE)
      axis(1, at=c(0,0.01,0.05,0.5,1,2), labels=c(0,0.01,0.05,0.5,1,2), cex.axis=1.3, font=1, font.lab=1)
      axis(2, at=c(0,0.01,0.05,0.5,1,2), labels=c(0,0.01,0.05,0.5,1,2), cex.axis=1.3, font=1, font.lab=1)
      par(font=1)
      box(lwd=1)
      savePlot(filename = paste(probe.name, ".modelled vs ref.png", sep=""), type = c("png"))     # save plot	
      
      
      windows()
      par(mar=c(5, 4.5, 4, 2) + 0.1)
      plot(bah[,3], type='h', lwd=2, cex.axis=1.3, cex.lab=1.3, font=1, font.lab=1, xaxt="n", xlab='Combinations of Rc, Ro', ylab='Function minimum')
      axis(1, at=c(0,511,1021,1531,2041,2551), labels=c("0, 0", "0.1, 0","0.2, 0","0.3, 0","0.4, 0","0.5, 0"), cex.axis=1.3, font=1, font.lab=1)
      box(lwd=1)
      savePlot(filename = paste(probe.name, ".combinations Rc,Ro.png", sep=""), type = c("png"))  # save plot	
      
      #print("NOT USED IN PLOTS: minimum with Rc < 0.1, see Bechtold et al., 2010", quote=F)
      
      optim.results <- bah1[minbah1,]
      optim.results[length(optim.results)+1] <- L.cable 
      
      colnames(optim.results) <- c("initial.Rc", "initial.R0", "min", "Rc", "R0", "L.cable")
      print(optim.results, quote=F)
      
      BEC.cal <- c(probe.name, Ropen, round(Kp, digits=4), round(optim.results$Rc, digits=4), round(optim.results$R0, digits=4))
      BEC.cal <- t(BEC.cal)
      BEC.cal <- data.frame(BEC.cal)
      colnames(BEC.cal) <- c("probe.mux", "Ropen", "Kp", "Rc", "R0")
       
      write_delim(BEC.cal, paste(probe.name, '_', "BEC.cal.txt", sep=''), delim='\t', append=F) 
       
       
} ## * ## end of for(ii in 1:length(list.cal))