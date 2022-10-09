##tdr.cond <- function(...) {
# Giulio Curioni
# 03 March 2017


## TDR WAVEFORM ANALYSIS FOR BULK ELECTRICAL CONDUCTIVITY (S/m)
   # the analysis is based on TDR waveforms taken with Campbell Scientific TDR100 and Campbell Scientific 3-wire probes
   # it should work with any probe provided that there is a drop in the reflection coefficient in the probe head
   
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
      
         source("C:\\R.wd\\loadall.R")        # load packages and functions
         loadall()
         graphics.off()
         options(warn=-1)       # suppress warnings

#check ### before running the script

### choose number of data points
n0 <- 2048

### 1) import calibration and measurement parameters. Note: names are case sensitive!
       probe.name <- "plab75.redtape4"    ## ## choose probe (it must match one of the existing probe.ID in tdr.probe.par.txt)
       mux.level <- "mux0"       ## ## choose the multiplexer level used for the measurements. Allowed options are: mux0, mux1, mux2   
       tdr.name <- "ATU.lab"     ## ## choose the TDR unit used for the measurements. Allowed options are: TDR19200, Dan.PhD, Giulio.PhD, ATU.lab, ATU.red, ATU.blue  
      

## import probe parameters:
    tdr.probe.par <- read.table("C:\\R.wd\\scripts\\tdr\\tdr.probe.par.txt", sep='\t', header=T)
    tdr.probe.par <- format(tdr.probe.par, digits=8)
    tdr.unit <- tdr.probe.par$tdr.unit
    probe.ID <- tdr.probe.par$probe.ID
    mux <- tdr.probe.par$mux
    tdr.probe.par <- data.frame(lapply(tdr.probe.par[,4:length(tdr.probe.par)], as.character), stringsAsFactors=FALSE)
    tdr.probe.par <- data.frame(lapply(tdr.probe.par[,1:length(tdr.probe.par)], as.numeric), stringsAsFactors=FALSE)
    tdr.probe.par <- cbind(tdr.unit, probe.ID, mux, tdr.probe.par)
    tdr.probe.par$tdr.unit <- as.character(tdr.probe.par$tdr.unit)
    tdr.probe.par$probe.ID <- as.character(tdr.probe.par$probe.ID)
    tdr.probe.par$mux <- as.character(tdr.probe.par$mux)
     
## define probe parameters corresponding to the chosen probe (see probe.name above)
    n.probe.ID.diff.units <- which(tdr.probe.par$probe.ID == probe.name)
    n.probe.ID <- which(tdr.probe.par$tdr.unit == tdr.name & tdr.probe.par$probe.ID == probe.name & tdr.probe.par$mux == mux.level)   # choose the probe identified by probe.name and tdr.name
        
    if(length(n.probe.ID) == 0) {stop("Double check that you chose the correct probe.name, mux.level and tdr.name in the script")
                          } else { 
                          if(is.na(n.probe.ID) == T) {stop("Double check that you chose the correct probe.name, mux.level and tdr.name in the script")} }
                          
    L.cable <- tdr.probe.par$L.cable[n.probe.ID]                # cable length (m)
    start.plot <- tdr.probe.par$start.plot[n.probe.ID]      # apparent length (m) corresponding to the parameter start in PCTDR (start of the plot)
    win.Lapp <- tdr.probe.par$win.Lapp[n.probe.ID]          # apparent length (m) of the window (entire TDR plot = length parameter in PCTDR)                                                 
    
    #n0 <- tdr.probe.par$n.data.points[n.probe.ID]            # number of data points
    Lcal <- tdr.probe.par$Lcal[n.probe.ID]                  # calibrated probe length (m)
    L0 <- tdr.probe.par$L0[n.probe.ID]                      # calibrated length between reference point and start point (m)
    start.plot.bec <- tdr.probe.par$start.plot.bec[n.probe.ID]  # apparent length (m) corresponding to the parameter start in PCTDR (start of the plot) for conductivity measurements
    win.Lapp.bec <- tdr.probe.par$win.Lapp.bec[n.probe.ID]      # apparent length (m) of the window (entire TDR plot = length parameter in PCTDR) for conductivity measurements
    long.dist.points <- tdr.probe.par$long.dist.points[n.probe.ID]  # number of long distance data points to be used for conductivity analysis
    Ropen <- tdr.probe.par$Ropen[n.probe.ID] 
    Kp <- tdr.probe.par$Kp[n.probe.ID]     # probe constant 
    Rc <- tdr.probe.par$Rc[n.probe.ID]     # additional resistance due to cables (e.g. Huisman et al., 2008)
    R0 <- tdr.probe.par$R0[n.probe.ID]     # additional resistance due to multiplexers, connectors, attachments etc. (e.g. Huisman et al., 2008)       

    
    
    
### ### if using custom parameters (e.g. window length for TDR inversion) manually modify and uncomment the parameters (uncomment only the ones to be modified):
    # note: the start.plot value does not make a difference in the calculation of permittivity so it is not added here

#        win.Lapp <- 3.5       
#        start.plot <- 7.6
#        n0 <- 2048 
#        Lcal <- 0.14969	
#        L0 <- 0.04334
#        L <-  0.1498				
#        z <- 	0.3381 
    
    
    n <- n0
    pto <- win.Lapp/n     # apparent length (m) of a single point
    
    if(length(long.dist.points) == 0) {long.dist.points <- 100}    


                                                                   
### 5) import waveform(s) - uncomment only the method used. IMPORTANT: SELECT ONLY MEASUREMENTS TAKEN WITH THE SAME PROBE! 

      ## if using in-house Matlab script:  
              path.name.bec <- data.frame(choose.files(caption = "Select CONDUCTIVITY waveforms"))   # select conductivity (bec) waveforms
              path.name.bec[,1] <- as.character(path.name.bec[,1])
              n.selected.waves <- length(path.name.bec[,1])     # number of selected waveforms to be analysed 
              
      #waves.Ka <- matrix(NA, nrow=n.selected.waves, ncol=n0+1)         # create file for the Ka waveforms   ### change ncol if necessary!
      waves.BEC <- matrix(NA, nrow=n.selected.waves, ncol=n0+1)        # create file for the BEC waveforms  ### change ncol if necessary!
              
      #output.Ka <- matrix(nrow=n.selected.waves, ncol=5)         # create file for the Ka output
      #output.w.dens <- matrix(nrow=n.selected.waves, ncol=11)    # create file for the Jung et al. 2013 method
      #output.fft <- matrix(nrow=n.selected.waves, ncol=16)       # create file for the TDR inversion output (fft)
      output.BEC <- matrix(nrow=n.selected.waves, ncol=3)        # create file for the BEC output
      
#      # import one waveform to see if it was taken with PCTDR or Matlab
#       file.name <- data.frame(strsplit(path.name.bec[1,1], ""))
#              file.name <- as.character(file.name[,1])
#              find.backlash <- which(file.name == '\\')
#              start.name <- find.backlash[length(find.backlash)]
#              rm(find.backlash)
#              file.name <- substring(path.name.bec[1,1], start.name+1, length(file.name))
#              Splits <- unlist(strsplit(file.name[1], "[.]"))
#              if(Splits[length(Splits)] == "DAT") { format.wave <- "PCTDR2" }
#              if(Splits[length(Splits)] == "dat") { format.wave <- "PCTDR3" }
#              if(Splits[length(Splits)] == "mat") { format.wave <- "Matlab" }
      
      
                            
      for(ii in 1:n.selected.waves) {       # select one or multiple measurements to analyse 
              file.name <- data.frame(strsplit(path.name.bec[ii,1], ""))
              file.name <- as.character(file.name[,1])
              find.backlash <- which(file.name == '\\')
              start.name <- find.backlash[length(find.backlash)]
              rm(find.backlash)
              file.name <- substring(path.name.bec[ii,1], start.name+1, length(file.name))          
              file.name <- data.frame(strsplit(file.name, "[.]"))      
              # find extension
              extension <- as.character(file.name[nrow(file.name),1])
              # remove extension
              file.name <- data.frame(file.name[1:(nrow(file.name)-1),1])
              # recombine to find name without extension
              file.name <- paste(file.name[,1], collapse='.')
                                                     
                       ##print file.name to track where the script stopped in case of errors
                       print(paste("file.name currently analysed: ", file.name, sep=''), quote=F) 
                       



            ## method 1) if using in-house Matlab script:
                             
            #wave <- data.frame(readMat(choose.files())); wave <- t(wave)
            #wave <- read.cb(sep='\t', header=F)
            #wave <- t(wave)
                                                                               
            if(extension == "mat") {
                wave <- data.frame(readMat(path.name.bec[ii,1]))
                #n <- ncol(wave)    
                wave <- t(wave)  
                current.wave <- c(file.name, wave)
                current.wave <- data.frame(current.wave)
                waves.BEC[ii, ] <- t(current.wave)        
                }

              
            ## method 2) if using PCTDR2 software:
            
            if(extension == "DAT") {
              wave <- read.delim(path.name.bec[ii,1], header=F)  #read.table("water.DAT", sep='\t', dec=',', header=F)
              wave <- wave[-1:-7,1]  # remove first 7 values (the eighth is always a zero)
              #n <- nrow(wave)             
              current.wave <- c(file.name, t(wave))
              current.wave <- data.frame(current.wave)
              waves.BEC[ii, ] <- t(current.wave)
              }

              
            ## method 3) if using PCTDR3 software:
            
            if(extension == "dat") {
              wave <- read.delim(path.name.bec[ii,1], sep=',', header=F, skip=4)  #read.table("water.DAT", sep='\t', dec=',', header=F)
              wave <- wave[,-1:-(ncol(wave)-n0)]  
              #n <- nrow(wave)             
              current.wave <- c(file.name, t(wave))
              current.wave <- data.frame(current.wave)
              waves.BEC[ii, ] <- t(current.wave)
              }

wave <- as.numeric(waves.BEC[ii, -1])
               
      x <- 1:n
      wave <- smooth(wave[1:n])   # smooth initial waveform to remove some roughness and stop at n data points (just in case) [note: this smoothing does not do much or nothing at all]

      ## create a data frame with X = n. data points, Y = reflection coefficient
      M <- cbind(x, wave)
      M <- data.frame(M)    
      colnames(M) <- c("x", "wave")
      
      
      
cond <- M$wave[(n-(long.dist.points-1)):n]     # take the last long.dist.points data points

RC <- mean(cond)			# measured reflection coefficient at long distance      ## ## I removed <<-
RC <- round(RC, digits=8)

Rcorr <- ((2*(RC+1))/(Ropen+1))-1	# corrected reflection coefficient

RL <- 50*((1+Rcorr)/(1-Rcorr))	  # load resistance       ## ## I removed <<-

BEC <- Kp/(RL-(L.cable*Rc+R0))	    # bulk electrical conductivity in S/m     ## ## I removed <<-
BEC <- round(BEC, digits=4)         

# plot

windows(w=12, h=7)

par(mar=c(5, 5, 4, 2) + 0.1)
plot(x, M$wave, type='l', lwd=1, ylim=c(-1,1), xlab='data points', ylab='reflection coefficient', font=1, font.lab=1, cex.lab=1.3, cex.axis=1.3)
            
            ## add title with BEC and long.dist.RC
            print.RC <- paste("long dist. reflection coeff:", as.character(RC), sep=" ")
            print.BEC <- paste("BEC (mS/m):", as.character(BEC*1000), sep=" ")
            mtext(text=print.RC, side = 3, line = 2, cex=1.2, col='black')
            mtext(text=print.BEC, side = 3, line = 0.5, cex=1.2, col='black')



      print(c("long distance reflection coefficient:", RC), quote=F);
      print(c("bulk electrical conductivity (mS/m):", BEC*1000), quote=F);


savePlot(filename = paste(file.name, ".BEC.emf", sep=""), type = c("emf"))    # save plot


                       # save the output in an object (file.name, BEC, long. dist refl. coeff)
                       output.BEC[ii,] <- c(file.name, BEC*1000, RC)
                       
                       
                       graphics.off()



}    # end of for(ii in 1:n.selected.waves)

                       # export the output to a .txt file
                       write.table(output.BEC,  paste(probe.name, ".BEC.txt", sep=''), sep='\t', row.names=F, col.names=c("file.name", "BEC.mS.m", "RC"), quote=F, append=T)

                 

##}