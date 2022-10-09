##tdr.perm.tta <- function(...) {
# Giulio Curioni
# 12 September 2018

     
## TDR WAVEFORM ANALYSIS FOR PERMITTIVITY (TRAVEL TIME ANALYSIS)
   # the analysis is based on TDR waveforms taken with Acclima TDR315 probes
   
     # waveform should be saved with the following format, keep '-' as the separator EXCEPT temperature (if present). Use '_' to contain temperature, if temperature is measured. See examples below. It is advisable to measure first using the Acclima software, at least to get the temperature.
     # Note that Acclima probes do not require multiplexers and start is always 0 ps.
     # i.e. probe.ID  -  window length (ps)  -  sample details  _  temperature +C  _  repetition
     # probe ID (i.e. probe name, as present in the file 'tdr.probe.par.txt')
     # window length parameter (ps), i.e. where the plot ends  (see example below)  [for conductivity use 1667840 ps]
       # for consistency with Campbell Scientific probes use approx 500m winLapp for conductivity, corresponding to 1667840 ps
     # sample ID + additional info separated by '-'
     # if present, temperature with 1 decimal point +C, contained within two '_'    (see example below)
     # repetition information  (always add an entry if temperature was taken, see examples)
     
     
     # examples:
     # e.g.1 sn8501652-0-15000-air_23.1C_rep1
     # e.g.2 sn8501652-0-20000-water_20.0C_rep1  # IMPORTANT: ALWAYS add something after temperature (e.g. use _1, _rep1 or _r1 etc. for no repetitions, i.e. one measurement only. Note: There is no need to add anything if temperature was not taken)

     # additional examples specific to conductivity measurements used for calibration:
     # sn8501652-1667840-sol2-1124_16.8C_rep4.DAT     # note that 'sol1', 'sol2' etc. and the ref. cond in mS/m MUST be present
     # sn8501652-1667840-air.open_21.0C_rep3.DAT      # always use 'air.open' for open measurements in air
     
     
        
         source("C:\\R.wd\\loadall.R")        # load packages and functions
         loadall()
         graphics.off()
         options(warn=-1)       # suppress warnings
                                                       
#check ### before running the script (note: search ### ### for important points for TDR inversion)



### 1) import waveform(s)
  
              path.name <- data.frame(choose.files(caption = "Select PERMITTIVITY waveforms"))       # select permittivity (Ka) waveforms
              path.name[,1] <- as.character(path.name[,1])
              n.selected.waves <- length(path.name[,1])     # number of selected waveforms to be analysed
      
      n0 <- 4000 ### ### initially set n0 to 4000 (probably max allowed for Acclima probes in 2018)        
      waves.Ka <- matrix(NA, nrow=n.selected.waves, ncol=n0+1)         # create file for the Ka waveforms   ### change ncol if necessary!
      #waves.BEC <- matrix(NA, nrow=n.selected.waves, ncol=n0+1)        # create file for the BEC waveforms  ### change ncol if necessary!
              
      output.Ka <- matrix(nrow=n.selected.waves, ncol=9)         # create file for the Ka output
      #output.w.dens <- matrix(nrow=n.selected.waves, ncol=11)    # create file for the Jung et al. 2013 method
      #output.fft <- matrix(nrow=n.selected.waves, ncol=16)       # create file for the TDR inversion output (fft)
      #output.BEC <- matrix(nrow=n.selected.waves, ncol=3)        # create file for the BEC output

      ## define important parameters
      j <- sqrt(as.complex(-1))   # or simply 0+1i  [better to use j in case I use i in a 'for' loop]
      c0 <- 2.9979e8              # speed of light (m/s)
      E0 <- 8.8542e-12            # absolute dielectric permittivity of free space (F/m)
      mu0 <- 1.2566e-6            # absolute magnetic permeability of free space (H/m)
      mu.abs <- mu0               # note: this assumes that the material is non-magnetic (H/m)
      mu <- mu.abs/mu0            # relative magnetic permeability

      ##    ii <- 1 
      
for(ii in 1:n.selected.waves) {       # select one or multiple measurements to analyse 
              file.name <- data.frame(strsplit(path.name[ii,1], ""))
              file.name <- as.character(file.name[,1])
              find.backlash <- which(file.name == '\\')
              start.name <- find.backlash[length(find.backlash)]
              rm(find.backlash)
              file.name <- substring(path.name[ii,1], start.name+1, length(file.name))          
              file.name <- data.frame(strsplit(file.name, "[.]"))      
              # find extension
              extension <- as.character(file.name[nrow(file.name),1])
              # remove extension
              file.name <- data.frame(file.name[1:(nrow(file.name)-1),1])
              # recombine to find name without extension
              file.name <- paste(file.name[,1], collapse='.')
                                                     
                       ##print file.name to track where the script stopped in case of errors
                       print(paste("file.name currently analysed: ", file.name, sep=''), quote=F) 

file.name1 <- data.frame(file.name)            
probe.name <- separate(file.name1, col=file.name, into=c("time", "sn"), sep = " ")
time <- probe.name$time
probe.name <- probe.name$sn

## import probe parameters:
    tdr.probe.par <- read_delim("C:\\R.wd\\scripts\\tdr\\tdr.probe.par.txt", delim='\t', col_names = TRUE)
    #tdr.probe.par <- format(tdr.probe.par, digits=8)
#    tdr.unit <- tdr.probe.par$tdr.unit
#    probe.ID <- tdr.probe.par$probe.ID
#    start.plot <- tdr.probe.par$start.plot[n.probe.ID]      # start of the plot (ps). Always keep it set to 0.
#    win.Lapp <- tdr.probe.par$win.Lapp[n.probe.ID]          # apparent length (m) of the window (entire TDR plot = length parameter in PCTDR)                                   
    n.probe.ID <- which(tdr.probe.par$probe.ID == probe.name) 
    Lcal <- tdr.probe.par$Lcal[n.probe.ID]                  # calibrated probe length (m)
    L0 <- tdr.probe.par$L0[n.probe.ID]                      # calibrated length between reference point and start point (m)
    win.Lapp.bec <- tdr.probe.par$win.Lapp.bec[n.probe.ID]      # long term (ps) of the window used for conductivity measurements
    long.dist.points <- tdr.probe.par$long.dist.points[n.probe.ID]  # number of long term data points to be used for conductivity analysis
    Ropen <- tdr.probe.par$Ropen[n.probe.ID] 
    Kp <- tdr.probe.par$Kp[n.probe.ID]     # probe constant 
    Rc <- tdr.probe.par$Rc[n.probe.ID]     # additional resistance due to cables (e.g. Huisman et al., 2008)
    R0 <- tdr.probe.par$R0[n.probe.ID]     # additional resistance due to multiplexers, connectors, attachments etc. (e.g. Huisman et al., 2008)       
    a <- tdr.probe.par$wire.diam[n.probe.ID]          # outer diameter of INNER conductor (m) [i.e. rod diameter]       
    b <- tdr.probe.par$ext.rod.spacing[n.probe.ID]    # inner diameter of OUTER conductor (m) 
    L <- tdr.probe.par$L[n.probe.ID]       # calibrated L after fft (TDR inversion)
    z <- tdr.probe.par$z[n.probe.ID]       # calibrated z after fft (TDR inversion)
    alpha <- tdr.probe.par$alpha[n.probe.ID]          # calibrated alpha after fft (TDR inversion)
    t0 <- tdr.probe.par$t0[n.probe.ID]                # calibrated t0 after fft (TDR inversion) 
wave0 <- read_csv(path.name[ii,1], col_names=T, na = "NA", skip = 6, n_max = Inf)
colnames(wave0) <- c("time", "amplitude")
if(!is.na(Ropen)) {
                  wave0$RC <- ((2*wave0$amplitude) - Ropen) / Ropen     # convert amplitudes to reflection coefficients
                  }


                wave <- wave0$RC
                n <- length(wave)  # actual length of waveform 
                win.tt <- wave0$time[nrow(wave0)]     # "length" of the plot in time units (ps). 
                win.Lapp <- (win.tt*10^-12)*c0        # convert ps to s and then time to apparent length (m)
                
                current.wave <- c(file.name, wave)
                current.wave <- t(current.wave)
                waves.Ka[ii, c(1:(n+1))] <- current.wave       


   
                ## ## input chosen length parameter (this is important in case something different from calibration is used)  

### 1) is the waveform a shorted waveform in air? YES = TRUE;  NO = FALSE.  Choose accordingly. 
       shorted <- F
       
### 2) are ~ horizontal tangent lines (not perfectly horizontal lines) necessary (e.g. in dry soil)? Choose accordingly.
# NOTE: this is set to FALSE for shorted measurements.
       use.end.tangent <- F      # for the end point (at the end of the rods)
       use.ref.tangent <- F      ## ## for the reference point (in the probe head) [this should ALWAYS remain set to FALSE for Acclima probes]
      
       if(shorted == T) {use.end.tangent <- FALSE}

### 3) is an interactive solution necessary to better plot the tangents? YES = TRUE;  NO = FALSE.
       interactive.sol <- F

          ## NOTE ON TANGENTS:
           # tangent 1 is the rising tangent near the reference point
           # tangent 2 is the rising tangent near the end point
           # tangent 3 is the horizontal or nearly horizontal tangent line near the reference point
           # tangent 4 is the horizontal or nearly horizontal tangent line near the end point
 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
      
      # set interactive.tan <- FALSE by default (no interactive solution)
      interactive.tan1 <- F     # keep this set to FALSE
      interactive.tan2 <- F
      interactive.tan3 <- F
      interactive.tan4 <- F
      
      
      if(interactive.sol == T) {
      # type 1 for tan1, 2 for tan2, 3 for tan3, 4 for tan4
      # multiple entries are accepted in case the interactive solution is for more than one tangent  
      answer.tan <- readline("\n\nManually select range of data points to plot tangents:\n\nType 1 for rising tangent1 (for refpoint)\nType 2 for rising (or descending) tangent2 (for endpoint)\nType 3 for ~horizontal tangent1 (for refpoint)\nType 4 for ~horizontal tangent2 (for endpoint)\n")
      
      # take a max of 15 characters to check if answer.tan = 1, 2, 3 or 4 (or a combination of these)
      substr.matrix <- matrix(nrow=15, ncol=1)
      for(iii in 1:15) {
            substr.matrix[iii,1] <- substr(answer.tan, iii, iii)
            }
      substr.matrix<- as.numeric(substr.matrix) # discard warnings (NA could be produced here but do not affect the script)
            
      # find if there are 1, 2, 3, 4 in substr.matrix and set interactive.tan      
      if(1 %in% substr.matrix == T) {interactive.tan1 <- T} 
      if(2 %in% substr.matrix == T) {interactive.tan2 <- T}
      if(3 %in% substr.matrix == T) {interactive.tan3 <- T}     ## ## keep it always to FALSE for Acclima probes
      if(4 %in% substr.matrix == T) {interactive.tan4 <- T}
      
      if(interactive.sol == T & interactive.tan1 == F & interactive.tan2 == F & interactive.tan3 == F & interactive.tan4 == F) {stop("Value not accepted. Re-run the script and choose 1 (tan1), 2 (tan2), 3 (tan3) or 4 (tan4) or a combination of these.")}
      
      # automatically set use.ref.tangent or use.end.tangent to TRUE if interactive.sol is TRUE and the user types 3 or 4
      if(interactive.tan3 == T) { use.ref.tangent <- T }
      if(interactive.tan4 == T) { use.end.tangent <- T } 
      
      } # end of if(interactive.sol == T)



# backup the initial settings 
interactive.sol.ini <- interactive.sol
interactive.tan1.ini <- interactive.tan1
interactive.tan2.ini <- interactive.tan2
interactive.tan3.ini <- interactive.tan3
interactive.tan4.ini <- interactive.tan4
probe.name.ini <- probe.name
shorted.ini <- shorted
#tdr.name.ini <- tdr.name
use.end.tangent.ini <- use.end.tangent
use.ref.tangent.ini <- use.ref.tangent

      

     
#    # import input functions from FFT calibration in air
#    input.fun <- read.table("C:\\R.wd\\scripts\\tdr\\tdr.input.fun.txt", sep='\t', header=F, stringsAsFactors=FALSE)
#    input.fun.header <- input.fun[1:3,-1]
#    input.fun.values <- input.fun[-1:-3,-1]
#    input.fun.values <- data.frame(lapply(input.fun.values[,], as.numeric), stringsAsFactors=FALSE)    
#    input.wave <- input.fun.values[,which(input.fun.header[1,] == tdr.name & input.fun.header[2,] == probe.name & input.fun.header[3,]== mux.level)]                 

              
    
    
### ### if using custom parameters for inversion manually modify and uncomment the parameters (uncomment only the ones to be modified):

#        Lcal <- 0.14904	
#        L0 <- 0.0514
#        L <-  0.1498				
#        z <- 	0.3381 
    
    pto <- win.Lapp/n     # step between two points in apparent length units
    pto_ps <- win.tt/n        # time step between two points       



## ## increase or decrease number of data points by interpolation using approx (this could be a good option in case long win.Lapp are used)
        # HOWEVER, FROM CALIBRATION MEASUREMENTS IN WATER INCREASING POINTS DID NOT IMPROVE THINGS AND IS TIME-CONSUMING  
      new.n <- n              
#              new.n <- 256 #8192             
#              wave <- approx(x=1:n, y=wave, n= new.n)
#              wave <- wave$y
#              bec.wave <- approx(x=1:n, y=bec.wave, n= new.n)
#              bec.wave <- bec.wave$y
#              n <-  new.n
#              pto <- win.Lapp/n    
                  
     x <- 1:n
     ### ###   IMPORTANT! APPLY SMOOTHING? loess is better. HOWEVER, by smoothing the waveform the high frequencies are eliminated!!
     # this is like using long cables or multiplexers! Ka is generally overestimated, possibly by 0.5 - 1 units! BETTER TO NOT USE UNLESS WAVEFORM IS VERY NOISY. IMPORTANT: startpoint shifts significantly to the left after smoothing!
#     #wave <- smooth(wave[1:n])   # smooth initial waveform to remove some roughness and stop at n data points (just in case) [note: this smoothing does not do much or nothing at all]
#     #bec.wave <- smooth(bec.wave[1:n])
                     
#             # smooth waveform using loess (better option). 
#              aa1 <- loess(wave ~ x, enp.target=100, control=loess.control(trace.hat = "approximate"))
#              wave <- predict(aa1)  # smoothed waveform

      ## create a data frame with X = n. data points, Y = reflection coefficient
      M <- cbind(x, wave) #, bec.wave)
      M <- data.frame(M)    
      colnames(M) <- c("x", "wave") #, "bec.wave")




                 









# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

## TTA: travel time analysis (traditional method for analysing TDR waveforms using tangent lines).  See Heimovaara and Bouten, 1990; Menziani et al., 1996.
# IMPORTANT: some bits of it are necessary for the frequency domain analysis, DO NOT DELETE
        
                                                                                                               
## ## prova <- function(...){      # useful in monitoring applications, discards dodgy measurements with a shrinked y-axis with respect to normal

## calculate first and second derivatives of the waveform and smooth results with cubic smoothing spline
        
        # first derivative
        # note: a max in the first derivative corresponds to an inflection point in a rising poriton of the original waveform
        # note: a min in the first derivative corresponds to an inflection point in a descending poriton of the original waveform
        wave.der1 <- diff(wave)
        x.der1 <- 1:(n-1)      # the first derivative has one point less than the waveform
        der1.spline <- smooth.spline(x=x.der1, y=wave.der1)
        der1 <- der1.spline$y         # save smoothed first derivative as a vector

        # second derivative
        # note: a max in the second derivative roughly corresponds to a min in the original waveform (roughly because of the smoothing)
        # note: a min in the second derivative roughly corresponds to a max in the original waveform (roughly because of the smoothing)
        wave.der2 <- diff(der1)
        x.der2 <- 1:(n-2)      # the second derivative has one point less than the first derivative
        der2.spline <- smooth.spline(x=x.der2, y=wave.der2)
        der2 <- der2.spline$y         # save smoothed second derivative as a vector







## find points of interest



# notes on rollapply: use a large enough width to avoid finding small max and min. If width is changed, which.max(...)==N, where N is half width
# see: http://stackoverflow.com/questions/6836409/finding-local-maxima-and-minima

# empirically, on the TDR-315 probes, it is necessary to discriminate peaks separated max 100 ps
max.roll.width <- nrow(wave0[wave0$time < 100,])
roll.width <- max.roll.width 
if(roll.width < 5) {roll.width <- 3}    # in times, if time step between points is larger (or close to) 100 ps try a small roll.width  
#if(n <= 2048) {roll.width <- 20} else {roll.width <- 100}
#roll.width <-  100    ### ### only in case of problems, modify this value (use a smaller number for large win.Lapp; e.g. for win.Lapp=2.5 try 100, for win.Lapp=10 use 20)
                       # in times, if time step between points (i.e. pto_ps) is larger or close to 100 ps try a very small roll.width but it could mean that it is not possible to analyse the waveform, either because there are too few points or because the total apparent length (win.Lapp) is too long  
half.width <- round(roll.width/2, digits=0)


## ## IN CASE OF PROBLEMS NEAR THE REFERENCE POINT (PROBE HEAD) (e.g. first horizontal line and/or first rising tangent are not correct), try in sequence:
    # 1) set change.n.exclusion <- TRUE  AND  plot.Lapp <- FALSE (see below)
         # n.exclusion excludes some points before ref point to avoid false identification of max and min. In most cases it should remain set to 0. The script will ask to select a point in the descending part of the waveform just before the reference point in the probe head (do not select a point in the cable).
    
    # 2) modify roll.width (see above)
    
    change.n.exclusion <- FALSE

    if(change.n.exclusion == TRUE) {         
           plot(wave)
           title("choose a point in the descending bit\nbetween cable and reference point")
           loc <- locator(1)
           n.exclusion <- loc$x      
           } else {
           n.exclusion <- 0 
           }  


## find local maxima and minima of the waveform
       xz <- as.zoo(wave)
       roll.max <- rollapply(xz, width=roll.width, function(wave) which.max(wave)==half.width)
       roll.min <- rollapply(xz, width=roll.width, function(wave) which.min(wave)==half.width)
       wave.max <- index(roll.max)[coredata(roll.max)]        ; wave.max <- wave.max[wave.max > n.exclusion]
       wave.min <- index(roll.min)[coredata(roll.min)]        ; wave.min <- wave.min[wave.min > n.exclusion]


## find local maxima and minima of first derivative corresponding to inflections of the waveform
       xz <- as.zoo(der1)
       roll1.max <- rollapply(xz, width=roll.width, function(der1) which.max(der1)==half.width)
       roll1.min <- rollapply(xz, width=roll.width, function(der1) which.min(der1)==half.width)
       der1.max <- index(roll1.max)[coredata(roll1.max)]      ; der1.max <- der1.max[der1.max > n.exclusion]
       der1.min <- index(roll1.min)[coredata(roll1.min)]     ; der1.min <- der1.min[der1.min > n.exclusion]


## find local maxima and minima of second derivative corresponding to inflections of the first derivative
# these should correspond to max and min of the original waveform
       xz <- as.zoo(der2)
       roll2.max <- rollapply(xz, width=roll.width, function(der2) which.max(der2)==half.width)
       roll2.min <- rollapply(xz, width=roll.width, function(der2) which.min(der2)==half.width)
       der2.max <- index(roll2.max)[coredata(roll2.max)]     ; der2.max <- der2.max[der2.max > n.exclusion]
       der2.min <- index(roll2.min)[coredata(roll2.min)]      ; der2.min <- der2.min[der2.min > n.exclusion]



## find inflection points where tangents will be drawn:
        #first inflection point (rising part of the waveform in the probe head (ref point))
        #second inflection point (rising part of the waveform at the end of the probe (end point))

        ## the first local max of the first derivative (excluding small local max) corresponds to the first inflection point (rising part just after the reference point in the probe head)
        ## note: the local max in der1 < 20% of the global max of der1 are discarded (this avoids to assign inflection1 to small local max due to noise in the waveform)
        global.max.der1 <- max(der1)                        # global max of der1
        ten.per.global.max.der1 <- global.max.der1 * 0.20   # 20% of global max of der1
        der1.max.sub.ref <- which(der1[der1.max] > ten.per.global.max.der1)      # take only local max of der1 larger than 20% of the global max of der1
        n.first.der1.max <- der1.max.sub.ref[1]      # index in der1.max.sub corresponding to inflection1 (first local max of subset)
        inflection1 <- der1.max[n.first.der1.max]    # this is the x index corresponding to the first inflection point
                
        global.min.der1 <- which.min(der1)

        if(shorted == F) {
              #measurement is not shorted.

               ## in theory the second inflection point is the absolute max of the first derivative AFTER the first inflection point:
               der1.max.sub.end <- der1.max[der1.max > global.min.der1]       # use global min instead of inflection1. In Acclima probes there is always a local min after startpoint. This is safer.
               der1.max.sub.end <- der1.max.sub.end[-1]                       # additionally, remove first local max after global min. This is very unlikely to be inflection2.
               max.der1.after.inflection1 <- max(der1[der1.max.sub.end])      # this is the global max of der1 after inflection1
               index.abs.max <- which.max(der1[der1.max.sub.end])
               inflection2 <- der1.max.sub.end[index.abs.max]                 # this is the x index corresponding to the second inflection point

               } else {
               #measurement is shorted.
               #second inflection point (descending part of the waveform at the end of the probe (end point), min of the first derivative AFTER inflection1)
               inflection2 <- which.min(der1[inflection1:length(der1)]) + inflection1     # this is the x index corresponding to the second inflection point
               }


## check inflections
#plot(der1, t='o'); points(x=c(inflection1, inflection2), y=c(0,0), pch="X", col="red"); points(x=der1.max.sub.end, y=rep(0,length(der1.max.sub.end)), pch="|", col="green")
#par(new=T)
#plot(wave, t='l', col='darkgreen')


## find min of the waveform corresponding to the reference point (head of the probe)


            # take the global minimum among the local minima before inflection1 corresponding to the drop in the probe head
#            n.vec.min.ref <- wave.min[wave.min < inflection1]
#            n.min.ref <- which.min(wave[n.vec.min.ref])
#            min.ref <- n.vec.min.ref[n.min.ref]  # index corresponding to the min in the probe head
#            min.ref.value <- M$wave[min.ref]     # min value corresponding to the min in the probe head
# for acclima probes:
            min.ref <- round(inflection1/3, digits=0)              # take one point that is clearly before the step pulse. Maybe better to not use the first point.
            min.ref.value <- M$wave[min.ref]      # min value corresponding to min.ref


## find max of the waveform roughly corresponding to the actual start point (bump after the probe head, i.e. start of the rods)
# note1: approx.startpoint is the same or very close to the actual start point in case of shorted measurements
# note2: if there is not a bump (i.e. waveform keeps going up, for example in air) this does not correspond to the true start point

            # this is the first local max after the first inflection point
            wave.max.after.inflection1 <- wave.max[wave.max > inflection1]
            approx.startpoint <- wave.max.after.inflection1[1]         # index roughly corresponding to actual start point
            peak.start <- approx.startpoint                            # index roughly corresponding to actual start point
            approx.startpoint.value <- M$wave[approx.startpoint]       # max value roughly corresponding to the actual start point



## find min of the waveform roughly corresponding to the end point


      # set rising.wave and auto.select.range1 and 2 to FALSE by default 
      rising.wave <- F
      auto.select.range1 <- F
      auto.select.range2 <- F
      
      # find global min of der1 between 2 inflections and find local minima (if any) in the waveform between global.min.der1 and inflection2
            sub.der1.min.between.inflections <- der1.min[der1.min > inflection1 & der1.min < inflection2] # subset der1 between inflections
            der1.min.between.inflections <- which.min(sub.der1.min.between.inflections) # local index of minima
            der1.min.between.inflections <- sub.der1.min.between.inflections[der1.min.between.inflections]  # global min of der1 between inflections
            
            # find if there are local min of the waveform between der1.min.between.inflections and inflection2
            local.wave.min <- wave.min[wave.min > der1.min.between.inflections  &  wave.min < inflection2]
      
      if(use.end.tangent == F) {
            
            if(length(local.wave.min) != 0){ 
                              # a local min in the waveform exists
                              # take the last local min before inflection2
                              min.end <- local.wave.min[length(local.wave.min)]
                              y.min.end <- M$wave[min.end]
                              } else {
                              # CASE 2.1 (see below) a local min in the waveform DOES NOT exist and one or more der1.min.between.inflections exist before inflection2 (e.g. in air) [note: it is assumed that a local min of der1 between inflections exists]
                              # take the last min of der1 before inflection2 and use tangent over a range of data points after der1.min.between.inflections 
                              der1.min.between.inflections <- sub.der1.min.between.inflections[length(sub.der1.min.between.inflections)]
                              min.end <- der1.min.between.inflections  # call it min.end 
                              y.min.end <- M$wave[min.end]
                              
                              if(shorted == F){
                                         use.end.tangent <- T 
                                         rising.wave <- T
                                         } 
                              }
       } else {
       # use.end.tangent = TRUE
       
       # CASE 2.2 (see below) a local min in the waveform DOES exist AND one or more der1.min.between.inflections exist after this local min and before inflection2 (i.e. it is still preferable to use the tangent method and not a perfectly horizontal tangent line across end point)
                                                                      
            if(length(local.wave.min) != 0) {
                              # a local min in the waveform exists
                              # take the last local min of the waveform before inflection2
                              min.end <- local.wave.min[length(local.wave.min)]
                              y.min.end <- M$wave[min.end]
                              
                              # take the first local min of der1 after min.end and before inflection2 (CASE 2.2, see below)
                              sub.der1.min.after.min.end <- der1.min[der1.min > min.end & der1.min < inflection2]
                              
                        if(length(sub.der1.min.after.min.end) != 0) {
                                  # a local min in der1 exists between min.end and inflection2
                                  der1.min.after.min.end <- sub.der1.min.after.min.end[1]
                                  min.end <- der1.min.after.min.end       # redefine min.end
                                  y.min.end <- M$wave[min.end]
                                  auto.select.range1 <- T
                                  } else {
                                  # CASE 2.3 (see below) a local min in the waveform DOES exist but there are no min in der1 between min.end and inflection2 (but it is still preferable to use the tangent method and not a perfectly horizontal tangent line across end point) 
                                  auto.select.range2 <- T
                                  } # end of if(length(sub.der1.min.after.min.end) != 0)
                              
                              } else {
                              # repeat the above for use.end.tangent = TRUE [this is necessary]
                              # CASE 2.1 (see below) a local min in the waveform DOES NOT exist and one or more der1.min.between.inflections exist before inflection2 (e.g. in air) [note: it is assumed that a local min of der1 between inflections exists]
                              # take the last min of der1 before inflection2 and use tangent over a range of data points after der1.min.between.inflections 
                              der1.min.between.inflections <- sub.der1.min.between.inflections[length(sub.der1.min.between.inflections)]
                              min.end <- der1.min.between.inflections  # call it min.end 
                              y.min.end <- M$wave[min.end]
                              
                              if(shorted == F){
                                         use.end.tangent <- T  
                                         rising.wave <- T
                                         }  
                              } # end of if(length(local.wave.min) != 0)    
       
       } # end of if(use.end.tangent == F)                       
       

            
            

## ## in case of persistent problems, use this for diagnostics (plot waveform, first and second derivative)

#            windows(w=12, h=7, title="first derivative")
#            plot(wave, t='l', ylim=c(-1,1), ylab=NA, yaxt='n')
#            par(new=T)
#            plot(wave.der1, t='l', col='red', ylab=NA, yaxt='n')
#            par(new=T)
#            plot(der1, t='l', col='blue', ylab=NA, yaxt='n')
#            legend("topleft",
#                c("waveform",
#                "first derivative",
#                "smoothed first derivative"
#                ),
#                text.col=c("black", "red", "blue")
#                )
#           points(x=der1.min, y=rep(0, length(der1.min)), pch='X', col='blue')    
#           points(x=wave.min, y=rep(0, length(wave.min)), pch=19, col='black')     
#            
#            windows(w=12, h=7, title="second derivative")
#
#            plot(wave, t='l', ylim=c(-1,1), ylab=NA, yaxt='n')
#            par(new=T)
#            plot(wave.der2, t='l', col='red', ylab=NA, yaxt='n')
#            par(new=T)
#            plot(der2, t='l', col='blue', ylab=NA, yaxt='n')
#            legend("topleft",
#                c("waveform",
#                "second derivative",
#                "smoothed second derivative"
#                ),
#                text.col=c("black", "red", "blue")
#                )
#            
#            windows()


               



                    



## identify tangents and lines

      # to facilitate the choice of data points used to create the tangents set plot.points <- TRUE below and use locator() if necessary 



      ## FIRST 2 TANGENTS (~ VERTICAL)




      ## ## rising tangent 1 (i.e. near reference point in probe head)

      # define a range of points near the inflection point to run a linear model (i.e. tangent)
      # note: this is critical to create a good input function for frequency domain analysis

if(interactive.tan1 == F) {
     # non-interactive solution
     
     y.inflection1 <- der1[inflection1]

     der1.min.before.inflection1 <- which(der1.min < inflection1)  # take all the local min in der1 before inflection1
     der1.min.before.inflection1 <- der1.min.before.inflection1[length(der1.min.before.inflection1)]  # take the last min of der1 before inflection1
     der1.min.before.inflection1 <- der1.min[der1.min.before.inflection1]   # take the actual index in der1 (not a index of a subset)

     der1.min.after.inflection1 <- which(der1.min > inflection1)  # take all the local min in der1 after inflection1
     der1.min.after.inflection1 <- der1.min.after.inflection1[1]  # take the first min of der1 after inflection1
     der1.min.after.inflection1 <- der1.min[der1.min.after.inflection1]   # take the actual index in der1 (not a index of a subset)

     sub.der1.inflection1 <- der1[ der1.min.before.inflection1 : der1.min.after.inflection1 ]    # subset between last min of der1 before inflection1 and first min of der1 after inflection1

     sub.auto.inflection1 <- which(sub.der1.inflection1 > y.inflection1 * 0.8)  +  der1.min.before.inflection1 # subset of data points around inflection1 to fit tangent (take values that are 80% from 0 to Y value of inflection1)
     low.lim1 <- inflection1 - sub.auto.inflection1[1]    # lower limit to draw tangent
     upp.lim1 <-  sub.auto.inflection1[length(sub.auto.inflection1)] - inflection1    # upper limit to draw tangent
     } else {
     # interactive solution
            
               windows(h=7, w=12)
               plot(wave, t='l')
               title("choose one point BEFORE inflection1 to draw the tangent")
               points(x=inflection1, y=M$wave[inflection1], pch=4, col='red')
               loc.before1 <- locator(1)
               low.lim1 <- round(inflection1 - loc.before1$x, digits=0) 
               
               plot(wave, t='l')
               title("choose one point AFTER inflection1 to draw the tangent")
               points(x=inflection1, y=M$wave[inflection1], pch=4, col='red')
               loc.after1 <- locator(1)
               upp.lim1 <- round(loc.after1$x - inflection1, digits=0)     
     } # end of if(interactive.tan1 == F)
     
     ### TANGENT 1: rising tangent 1 (i.e. near reference point in probe head)
     # if the procedure above is not giving an appropriate range of data points AND an interactive solution is not acceptable (e.g. in monitoring applications with many waveforms to analyse), uncomment and manually modify low.lim and upp.lim values:
     # low.lim1 <- 12
     # upp.lim1 <- 12
     
     sub.tangent1 <- subset(M, M$x >= inflection1 - low.lim1  &  M$x <= inflection1 + upp.lim1)
     sub.tangent1 <- data.frame(sub.tangent1)
     
                            # run linear model (i.e. tangent)
                            # rising tangent 1 across inflection1
                            tryCatch(lm.tangent1 <- lm(sub.tangent1$wave ~ sub.tangent1$x, data=sub.tangent1), error = function(err) {
                            print("Error in lm.tangent1. Try modifying low.lim1 and upp.lim1 and check inflection1 value.")
                            }
                      )   




      ## ## rising (or descending in case of shorted measurements) tangent 2 (i.e. near the end point)

      # define a range of points near the inflection point to run a linear model (i.e. tangent)
     
if(interactive.tan2 == F) {       
     # non interactive solution
     
  if(shorted == F) {
  # measurement is not shorted
     
     y.der1.inflection2 <- der1[inflection2]

     der1.min.before.inflection2 <- which(der1.min < inflection2)  # take all the local min in der1 before inflection2
     der1.min.before.inflection2 <- der1.min.before.inflection2[length(der1.min.before.inflection2)]  # take the last min of der1 before inflection2
     der1.min.before.inflection2 <- der1.min[der1.min.before.inflection2]   # take the actual index in der1 (not a index of a subset)

     der1.min.after.inflection2 <- which(der1.min > inflection2)  # take all the local min in der1 after inflection2
     der1.min.after.inflection2 <- der1.min.after.inflection2[1]  # take the first min of der1 after inflection2
     der1.min.after.inflection2 <- der1.min[der1.min.after.inflection2]   # take the actual index in der1 (not a index of a subset)

     if(is.na(der1.min.after.inflection2 == T)) {der1.min.after.inflection2 <- n}
     
     sub.der1.inflection2 <- der1[ der1.min.before.inflection2 : der1.min.after.inflection2 ]    # subset between last min of der1 before inflection2 and first min of der1 after inflection2

     sub.auto.inflection2 <- which(sub.der1.inflection2 > y.der1.inflection2 * 0.8)  +  der1.min.before.inflection2 # subset of data points around inflection2 to fit tangent (take values that are 80% from 0 to Y value of inflection2)
     low.lim2 <- inflection2 - sub.auto.inflection2[1]    # lower limit to draw tangent
     upp.lim2 <- sub.auto.inflection2[length(sub.auto.inflection2)] - inflection2    # upper limit to draw tangent
   } else {
   # measurement is shorted

     y.der1.inflection2 <- der1[inflection2]                   
     der1.max.after.inflection2 <- which(der1.max > inflection2)  # take all the local max in der1 after inflection2
     der1.max.after.inflection2 <- der1.max.after.inflection2[1]  # take the first max of der1 after inflection2
     der1.max.after.inflection2 <- der1.max[der1.max.after.inflection2]   # take the actual index in der1 (not a index of a subset)

     sub.der1.inflection2 <- der1[ inflection1 : der1.max.after.inflection2 ]    # subset between inflection1 and first max of der1 after inflection2

     sub.auto.inflection2 <- which(sub.der1.inflection2 < y.der1.inflection2 * 0.8)  +  inflection1 # subset of data points around inflection2 to fit tangent (take values that are 80% from 0 to Y value of inflection2)
     low.lim2 <- inflection2 - sub.auto.inflection2[1]    # lower limit to draw tangent
     upp.lim2 <- sub.auto.inflection2[length(sub.auto.inflection2)] - inflection2    # upper limit to draw tangent
   
   
   } # end of if(shorted == F)  
     } else {
     # interactive solution
               
               windows(h=7, w=12)
               plot(wave, t='l')
               title("choose one point BEFORE inflection2 to draw the tangent")
               points(x=inflection2, y=M$wave[inflection2], pch=4, col='red')
               loc.before2 <- locator(1)
               low.lim2 <- round(inflection2 - loc.before2$x, digits=0) 
               
               plot(wave, t='l')
               title("choose one point AFTER inflection2 to draw the tangent")
               points(x=inflection2, y=M$wave[inflection2], pch=4, col='red')
               loc.after2 <- locator(1)
               upp.lim2 <- round(loc.after2$x - inflection2, digits=0)     
     } # end of if(interactive.tan2 == F)
          
     ### TANGENT 2: rising (or descending in case of shorted measurements) tangent 2 (i.e. near the end point)
     # if the procedure above is not giving an appropriate range of data points AND an interactive solution is not acceptable (e.g. in monitoring applications with many waveforms to analyse), uncomment and manually modify low.lim and upp.lim values:
     # low.lim2 <- 10
     # upp.lim2 <- 10
     
     sub.tangent2 <- subset(M, M$x >= inflection2 - low.lim2  &  M$x <= inflection2 + upp.lim2)
     sub.tangent2 <- data.frame(sub.tangent2)
     
                      # run linear model (i.e. tangent)
                      # rising tangent 2 across inflection2
                      tryCatch(lm.tangent2 <- lm(sub.tangent2$wave ~ sub.tangent2$x, data=sub.tangent2), error = function(err) {
                            print("Error in lm.tangent2. Try modifying low.lim2 and upp.lim2 and check inflection2 value.")
                            }
                      )     











      ## SECOND 2 TANGENTS (~ HORIZONTAL)




      # ref point
    
      # 1) by default use horizontal base lines (i.e. perfectly horizontal line) [calculate this anyway]
      x.ref <- rep.int(min.ref.value,n)                  # perfectly horizontal line at the ref point (probe head)
      
      
      # 2) ~ horizontal tangent line (not necessarily horizontal) for reference point
      if(use.ref.tangent == T) {

     ## horizontal TANGENT line at the ref point (as in Menziani) - interactive solution (the alternative is already calculated)
     if(interactive.tan3 == TRUE){
     # interactive solution
     
     # define a range of points near the waveform minimum (ref point) to run a linear model (i.e. tangent)
               
               windows(h=7, w=12)
               plot(wave, t='l')
               title("choose one point BEFORE min.ref to draw the tangent")
               points(x=min.ref, y=M$wave[min.ref], pch=3, col='red')
               loc.before3 <- locator(1)
               low.lim3 <- round(min.ref - loc.before3$x, digits=0) 
               
               plot(wave, t='l')
               title("choose one point AFTER min.ref to draw the tangent")
               points(x=min.ref, y=M$wave[min.ref], pch=3, col='red')
               loc.after3 <- locator(1)
               upp.lim3 <- round(loc.after3$x - min.ref, digits=0)     
               } # end of if(interactive.tan3 == T)
               
     ## ## TANGENT 3: ~ horizontal tangent line (not necessarily horizontal) for reference point
     # if the procedure above is not giving an appropriate range of data points (e.g. in monitoring applications with many waveforms to analyse), uncomment and manually modify low.lim and upp.lim values
      # low.lim3 <- 10
      # upp.lim3 <- 10
     
     sub.tangent3 <- subset(M, M$x >= min.ref - low.lim3  &  M$x <= min.ref + upp.lim3)
     sub.tangent3 <- data.frame(sub.tangent3)
     
                      # run linear model (i.e. tangent)
                      # nearly horizontal tangent across min in probe head
                      tryCatch(lm.tangent3 <- lm(sub.tangent3$wave ~ sub.tangent3$x, data=sub.tangent3), error = function(err) {
                            print("Error in lm.tangent3. Try modifying low.lim3 and upp.lim3 and check min.ref value.")
                            }
                      )

     } # end of if(use.ref.tangent == T)     






      # end point

      #  1) by default use horizontal base line (i.e. perfectly horizontal line) for end point [calculate this anyway]
           
           
            if(shorted == F) {
            # measurement is not shorted
            x.end <- rep.int(y.min.end,n)            # perfectly horizontal line at the end point (end of the probe)
            } else {
            # measurement is shorted
            x.end <- rep.int(approx.startpoint.value,n)        # horizontal line for shorted measurements (at the bump roughly corresponding to start point)
            }
            
            
            
      # 2) ~ horizontal tangent line (not necessarily horizontal) for end point. Useful for example in dry soils.
      
                  
      # CASE 2.1 - there is no local min in waveform between der1.min.between.inflections and inflection2 (i.e. waveform keeps rising).
      if(use.end.tangent == T & rising.wave == T) {
                  
      # define a range of points near der1.min.between.inflections (min.end, see above) to run a linear model (i.e. tangent):
                  
     if(interactive.tan4 == F) {
     # non interactive solution
     
     y.der1.min.between.inflections <- der1[min.end]   # Y value of min in der1
     y.der1.inflection2 <- der1[inflection2]           # Y value of inflection2
     diff.y.min.end.inflection2 <- abs(y.der1.inflection2 - y.der1.min.between.inflections)   # absolute difference between Y values of inflection2 and min in der1

     sub.der1.min.end.2.1 <- der1[ inflection1 : inflection2 ]    # subset der1 between inflections

     sub.auto.min.end.2.1 <- which(sub.der1.min.end.2.1 < (y.der1.min.between.inflections + diff.y.min.end.inflection2 * 0.10) )  +  inflection1 # subset of data points around min.end to fit tangent (take values that are 10% from y.der1.min.between.inflections to Y value of inflection2)
     low.lim4 <- min.end - sub.auto.min.end.2.1[1]    # lower limit to draw tangent
     upp.lim4 <- sub.auto.min.end.2.1[length(sub.auto.min.end.2.1)] - min.end    # upper limit to draw tangent
     } else {
     # interactive solution
            
               windows(h=7, w=12)
               plot(wave, t='l')
               points(x=min.end, y=y.min.end, pch=3, col='red')
               title("choose one point near min.end to draw the tangent for end point\n(typically on the left of min.end)")
               points(x=min.end, y=M$wave[min.end], pch=3, col='red')
               loc.before4 <- locator(1)
               low.lim4 <- round(min.end - loc.before4$x, digits=0) 
               
               plot(wave, t='l')
               points(x=min.end, y=y.min.end, pch=3, col='red')
               title("choose one point near min.end AFTER the previous one\nto draw the tangent for end point")
               points(x=min.end, y=M$wave[min.end], pch=3, col='red')
               loc.after4 <- locator(1)
               upp.lim4 <- round(loc.after4$x - min.end, digits=0)     
     } # end of if(interactive.tan4 == F)

     } # end of if(use.end.tangent == T & rising.wave == T)





      # CASE 2.2 - there is a local min in waveform between der1.min.between.inflections and inflection2 and there is a local min in der1 after this min and before inflection2
      if(use.end.tangent == T & auto.select.range1 == T) {
                  
      # define a range of points between min.end and der1.min.after.min.end (see above) to run a linear model (i.e. tangent):
                  
     if(interactive.tan4 == F) {
     # non interactive solution
                                                
     low.lim4 <- 5    # lower limit to draw tangent (5 is an arbitrary number of data points to use to the left of min.end)
     upp.lim4 <- (der1.min.after.min.end - min.end) + 15    # upper limit to draw tangent (5 is an arbitrary number of data points to use to the right of der1.min.after.min.end)
     } else {
     # interactive solution
            
               windows(h=7, w=12)
               plot(wave, t='l')
               points(x=min.end, y=y.min.end, pch=3, col='red')
               title("choose one point near min.end to draw the tangent for end point\n(typically on the left of min.end)")
               points(x=min.end, y=M$wave[min.end], pch=3, col='red')
               loc.before4 <- locator(1)
               low.lim4 <- round(min.end - loc.before4$x, digits=0) 
               
               plot(wave, t='l')
               points(x=min.end, y=y.min.end, pch=3, col='red')
               title("choose one point near min.end AFTER the previous one\nto draw the tangent for end point")
               points(x=min.end, y=M$wave[min.end], pch=3, col='red')
               loc.after4 <- locator(1)
               upp.lim4 <- round(loc.after4$x - min.end, digits=0)     
     } # end of if(interactive.tan4 == F)
     
     } # end of if(use.end.tangent == T & auto.select.range1 == T)





      # CASE 2.3 - there is a local min in waveform between der1.min.between.inflections and inflection2 but there is no local min in der1 after this min and before inflection2
      if(use.end.tangent == T & auto.select.range2 == T) {
                  
      # define a range of points around min.end to run a linear model (i.e. tangent):
                  
     if(interactive.tan4 == F) {
     # non interactive solution
     
     # if a use.end.tangent = T it means that the tangent should be drawn to the right of min.end:                                           
     low.lim4 <- 0     # lower limit to draw tangent (0 is an arbitrary number of data points to use to the left of min.end)
     upp.lim4 <- 25    # upper limit to draw tangent (25 is an arbitrary number of data points to use to the right of min.end)
     } else {
     # interactive solution
            
               windows(h=7, w=12)
               plot(wave, t='l')
               points(x=min.end, y=y.min.end, pch=3, col='red')
               title("choose one point near min.end to draw the tangent for end point\n(typically on the left of min.end)")
               points(x=min.end, y=M$wave[min.end], pch=3, col='red')
               loc.before4 <- locator(1)
               low.lim4 <- round(min.end - loc.before4$x, digits=0) 
               
               plot(wave, t='l')
               points(x=min.end, y=y.min.end, pch=3, col='red')
               title("choose one point near min.end AFTER the previous one\nto draw the tangent for end point")
               points(x=min.end, y=M$wave[min.end], pch=3, col='red')
               loc.after4 <- locator(1)
               upp.lim4 <- round(loc.after4$x - min.end, digits=0)     
     } # end of if(interactive.tan4 == F)
     
     } # end of if(use.end.tangent == T & auto.select.range2 == T)
          
     
     if(use.end.tangent == T) {     
     ### TANGENT 4: ~ horizontal tangent line (not necessarily horizontal) for end point. Useful for example in dry soils.
     # if the procedure above is not giving an appropriate range of data points AND an interactive solution is not acceptable (e.g. in monitoring applications with many waveforms to analyse), uncomment and manually modify low.lim and upp.lim values:
     # low.lim4 <- 10
     # upp.lim4 <- 10
     
     sub.tangent4 <- subset(M, M$x >= min.end - low.lim4  &  M$x <= min.end + upp.lim4)
     sub.tangent4 <- data.frame(sub.tangent4)
     
                      # run linear model (i.e. tangent)
                      # nearly horizontal tangent near the end point
                      tryCatch(lm.tangent4 <- lm(sub.tangent4$wave ~ sub.tangent4$x, data=sub.tangent4), error = function(err) {
                            print("Error in lm.tangent4. Try modifying low.lim4 and upp.lim4 and check min.end value.")
                            }
                      )  

     } # end of if(use.end.tangent == T)












      
      
      
      
      
      
      
      
      
      
      
      
                   
      
                   
      



             


## plots


windows(w=12, h=7, title="waveform TTA")


## ## if plot.Lapp = TRUE the x-axis shows the apparent length, if FALSE the x-axis shows the number of data points
plot.Lapp <- F

if(plot.Lapp == T){
## plot apparent length on x-axis

      lp <- win.Lapp/n
      lapp <- seq(lp,win.Lapp,lp)

      if(use.end.tangent == F) {
      # use horizontal base lines.

          par(mar=c(5, 5, 4, 2) + 0.1)
          plot(M$x, M$wave, type='l', pch='.', ylim=c(-1,1), axes=F, xlab=NA, ylab=NA)
          # plot ~ vertical tangent lines
          abline(coef(lm.tangent1))
          abline(coef(lm.tangent2))
          points(x=sub.tangent1$x, y=sub.tangent1$wave, col='red', pch='.', cex=1.5)        
          points(x=sub.tangent2$x, y=sub.tangent2$wave, col='red', pch='.', cex=1.5)        


          par(new=TRUE)


          ## choose this to plot horizontal TANGENT lines (i.e. not necessarily horizontal)
      
          if(use.ref.tangent == TRUE) {
                  abline(coef(lm.tangent3)) 	# to draw horizontal TANGENT for reference point 
                  points(x=sub.tangent3$x, y=sub.tangent3$wave, col='red', pch='.', cex=1.5)
                  } else {
          plot(M$x, x.ref, type='l', pch='.', ylim=c(-1,1), axes=F, xlab=NA, ylab=NA)		# for plotting horizontal line
          } # end of if(use.ref.tangent == TRUE)
          par(new=T)
          plot(lapp, x.end, type='l', pch='.', ylim=c(-1,1),  xlab='apparent length (m)', ylab='reflection coefficient')		# for plotting horizontal line
                 

          } else {
          # use ~horizontal tangent line for the end point (an option is included also for ref point).

          par(mar=c(5, 5, 4, 2) + 0.1)
          plot(M$x, M$wave, type='o', pch='.', ylim=c(-1,1), axes=F, xlab=NA, ylab=NA)
          # plot ~ vertical tangent lines
          abline(coef(lm.tangent1))
          abline(coef(lm.tangent2))
          points(x=sub.tangent1$x, y=sub.tangent1$wave, col='red', pch='.', cex=1.5)        
          points(x=sub.tangent2$x, y=sub.tangent2$wave, col='red', pch='.', cex=1.5)
   

          par(new=TRUE)


          ## choose this to plot horizontal TANGENT lines (i.e. not necessarily horizontal)
      
          if(use.ref.tangent == TRUE) {
                  abline(coef(lm.tangent3)) 	# to draw horizontal TANGENT for reference point 
                  points(x=sub.tangent3$x, y=sub.tangent3$wave, col='red', pch='.', cex=1.5)
                  } else {
          plot(M$x, x.ref, type='l', pch='.', ylim=c(-1,1), axes=F, xlab=NA, ylab=NA)		# for plotting horizontal line
          } # end of if(use.ref.tangent == TRUE)
      
          abline(coef(lm.tangent4)) 	  # to draw horizontal TANGENT for end point
          points(x=sub.tangent4$x, y=sub.tangent4$wave, col='red', pch='.', cex=1.5)
          

          par(new=TRUE)


          plot(lapp, x.ref, type='n', pch='.', ylim=c(-1,1), axes=F, xlab=NA, ylab=NA)		# leave it with type='n' (suppress plot)
          par(new=T)
          plot(lapp, x.end, type='n', pch='.', ylim=c(-1,1),  xlab='apparent length (m)', ylab='reflection coefficient')		# leave it with type='n' (suppress plot)

          } # end of if(use.end.tangent == F)


} else {
## to plot data points on x-axis
      
      if(use.end.tangent == F) {
      # use horizontal base line for end point (an option is included also for ref point).

          par(mar=c(5, 5, 4, 2) + 0.1)
          plot(M$x, M$wave, type='l', pch='.', ylim=c(-1,1), axes=F, xlab=NA, ylab=NA)
          # plot ~ vertical tangent lines
          abline(coef(lm.tangent1))
          abline(coef(lm.tangent2))
          points(x=sub.tangent1$x, y=sub.tangent1$wave, col='red', pch='.', cex=1.5)        
          points(x=sub.tangent2$x, y=sub.tangent2$wave, col='red', pch='.', cex=1.5)        
          

          par(new=TRUE)


          ## choose this to plot horizontal TANGENT lines (i.e. not necessarily horizontal)
      
          if(use.ref.tangent == TRUE) {
                  abline(coef(lm.tangent3)) 	# to draw horizontal TANGENT for reference point 
                  points(x=sub.tangent3$x, y=sub.tangent3$wave, col='red', pch='.', cex=1.5)
                  } else {
          plot(M$x, x.ref, type='l', pch='.', ylim=c(-1,1), axes=F, xlab=NA, ylab=NA)		# for plotting horizontal line
          } # end of if(use.ref.tangent == TRUE)
          par(new=T)
          plot(M$x, x.end, type='l', pch='.', ylim=c(-1,1),  xlab='data points', ylab='reflection coefficient')		# for plotting horizontal line
           
          } else {
          # use ~horizontal tangent line for the end point (an option is included also for ref point).

          par(mar=c(5, 5, 4, 2) + 0.1)
          plot(M$x, M$wave, type='o', pch='.', ylim=c(-1,1), axes=F, xlab=NA, ylab=NA)
          # plot ~ vertical tangent lines
          abline(coef(lm.tangent1))
          abline(coef(lm.tangent2))
          points(x=sub.tangent1$x, y=sub.tangent1$wave, col='red', pch='.', cex=1.5)        
          points(x=sub.tangent2$x, y=sub.tangent2$wave, col='red', pch='.', cex=1.5)        
          

          par(new=TRUE)

          
          ## choose this to plot horizontal TANGENT lines (i.e. not necessarily horizontal)

          if(use.ref.tangent == TRUE) {
                  abline(coef(lm.tangent3)) 	# to draw horizontal TANGENT for reference point
                  points(x=sub.tangent3$x, y=sub.tangent3$wave, col='red', pch='.', cex=1.5)
                  } else {
          plot(M$x, x.ref, type='l', pch='.', ylim=c(-1,1), axes=F, xlab=NA, ylab=NA)		# for plotting horizontal line
          } # end of if(use.ref.tangent == TRUE)

          abline(coef(lm.tangent4)) 	  # to draw horizontal TANGENT for end point
          points(x=sub.tangent4$x, y=sub.tangent4$wave, col='red', pch='.', cex=1.5)

          par(new=TRUE)


          plot(M$x, x.ref, type='n', pch='.', ylim=c(-1,1), axes=F, xlab=NA, ylab=NA)		# leave it with type='n' (suppress plot)
          par(new=T)
          plot(M$x, x.end, type='n', pch='.', ylim=c(-1,1),  xlab='data points', ylab='reflection coefficient')		# leave it with type='n' (suppress plot)

          } # end of if(use.end.tangent == F)

} # end of if(plot.Lapp == T)



               





                 
      ## ## note: to check waveforms and zoom plots, use zm function ( better to do zm(type='session') or zm(t='s') )








      ## FIND REFERENCE AND END REFLECTION POINTS (y=mx+q)
      int1 <- summary(lm.tangent1)$coefficients[1,1]		# this extracts the Intercept from the linear model
      slope1 <- summary(lm.tangent1)$coefficients[2,1]	# this extracts the Slope

      int2 <- summary(lm.tangent2)$coefficients[1,1]
      slope2 <- summary(lm.tangent2)$coefficients[2,1]

      if(use.end.tangent==F){
      # do not use horizontal TANGENT lines (i.e. use horizontal lines)
            refpoint <- (min.ref.value-int1)/slope1       # reference point

            if(shorted == F){
            endpoint <- (y.min.end-int2)/slope2       # if using the horizontal line  (non-shorted)
            } else {
            endpoint <- (approx.startpoint.value-int2)/slope2   # if using the horizontal line  (shorted)
            }


      } else {
      # use horizontal TANGENT lines (i.e. not necessarily horizontal)
            if(use.ref.tangent == T){
            # if I use horizontal TANGENT line for ref point:
                 int3 <- summary(lm.tangent3)$coefficients[1,1]
                 slope3 <- summary(lm.tangent3)$coefficients[2,1]
                 refpoint <- round((int3-int1)/(slope1-slope3))
                 } else {
                 refpoint <- (min.ref.value-int1)/slope1      # reference point
                 }
  
            int4 <- summary(lm.tangent4)$coefficients[1,1]
            slope4 <- summary(lm.tangent4)$coefficients[2,1]
            
            endpoint <- (int4-int2)/(slope2-slope4)    # if using the horizontal TANGENT line (non-shorted)
            }












               




      ## FIND PERMITTIVITY:
      refpoint <- round(refpoint, digits=6)
      endpoint <- round(endpoint, digits=6)
#      refM <- pto*refpoint
#      startM <- refM+L0
#      startpoint <- startM/pto
#      startpoint <- round(startpoint, digits=6)
#      endM <- endpoint*pto
#      Ka <- ((endM-startM)/Lcal)^2          # I removed <<-
#      Ka <- round(Ka, digits=2)             # round value to 2 digital digits
      pto0 <- (2*L0)/pto                     # offset from refpoint to startpoint (in n. of data points units). It must be multiplied by 2 to account for two way travel.
      startpoint <- refpoint+pto0
      startpoint <- round(startpoint, digits=6)
      tt <- endpoint-startpoint              # in n. of data points
      tt <- tt*pto_ps*10^-12                 # convert to times (s)
      v <- (2*Lcal)/tt                       # velocity in m/s
      Ka <- round((c0/v)^2, digits=2)



      # calculate volumetric water content using Topp et al. 1980:
      vwc.topp <- -5.30e-2 + 2.92e-2 * Ka - 5.5e-4 * Ka^2 + 4.3e-6 * Ka^3
      vwc.topp <- vwc.topp*100   # convert values to percentages
      vwc.topp <- round(vwc.topp, digits=1)       # round value to 1 decimal digits



      ## ## add points and print apparent permittivity/VWC by Topp on the plot
      
      ## ## if critical points and Ka and VWC are not wanted on the plot set plot.points <- F  (or comment unwanted bits below)
      
plot.points <- TRUE
      
if(plot.points == T) {
      
      # add critical points to the graph (useful for example to choose a suitable range of data points for rising tangents) 
#      if(use.end.tangent == F) {
               
      points(x=c(min.ref,
            refpoint, 
            inflection1,
            startpoint, 
            min.end,
            endpoint, 
            inflection2
            ), 
            y = c(M$wave[min.ref],
            M$wave[refpoint],
            M$wave[inflection1],
            M$wave[startpoint], 
            M$wave[min.end],
            M$wave[endpoint],
            M$wave[inflection2]
            ), pch=c(3,1,4,1,3,1,4), col=c('red', 'blue', 'red', 'blue', 'red', 'blue', 'red'))
            
            ## add labels to some points
            text(x=c(refpoint, 
            startpoint, 
            endpoint), 
            y = c(M$wave[refpoint],
            M$wave[startpoint], 
            M$wave[endpoint]), labels=c("refpoint", "startpoint", "endpoint"), pos = 4, offset = 0.5, cex=0.8)
            
            
            
            ## add title with apparent permittivity and VWC calculated by Topp et al., 1980
           # print.Ka <- paste("Apparent permittivity", " = ", as.character(Ka), sep=" ")
           # print.vwc.topp <- paste("Topp VWC (%):", as.character(vwc.topp), sep=" ")
           # mtext(text=print.Ka, side = 3, line = 2, cex=1.2, col='black')
           # mtext(text=print.vwc.topp, side = 3, line = 0.5, cex=0.8, col='black')
           # mtext(text=bquote(~K[a] == .(Ka)), side = 3, line = 2, cex=1.4, col='black')
            mtext(text=bquote(paste(~K[a] == .(Ka), "     ", ~theta[Topp] == .(vwc.topp), "%", sep=' ')), side = 3, line = 0, cex=1.7, col='black')

} # end of if(plot.point == T)


      savePlot(filename = paste(file.name, ".TTA.emf", sep=""), type = c("emf"))    # save the plot (TTA = travel time analysis)


      print(c("refpoint:", round(refpoint, digits=6)), quote=F);
      print(c("startpoint:", round(startpoint, digits=6)), quote=F);
      print(c("endpoint:", round(endpoint, digits=6)), quote=F);
      print(c("Topp VWC (%):", vwc.topp), quote=F)
      print(c("apparent permittivity:", Ka), quote=F) 
      ## ##}    # end of prova


absmin <- min(as.numeric(M$wave))
absmax <- max(as.numeric(M$wave))
differenza <- absmax-absmin
## ##if(differenza>0.2) {try(prova())         # useful in monitoring applications, discards doggy measurements with a shrinked y-axis with respect to normal
## ##} else {
## ## Ka <- NA}                 # I removed <<- 

                 
              
                    #stop()

                            


                       # save the output in an object (file.name, Ka, VWC (Topp))
                       output.Ka[ii,] <- c(file.name, Ka, vwc.topp, refpoint, endpoint, win.Lapp, tt, pto, pto_ps) 
                       #colnames(output.Ka) <- c("file.name", "Ka", "vwc.topp")                              



             


 #graphics.off()










# reset initial settings for next waveform analysis (in case multiple waveforms are analysed at the same time)
interactive.sol <- interactive.sol.ini
interactive.tan1 <- interactive.tan1.ini
interactive.tan2 <- interactive.tan2.ini
interactive.tan3 <- interactive.tan3.ini
interactive.tan4 <- interactive.tan4.ini
probe.name <- probe.name.ini
shorted <- shorted.ini
#tdr.name <- tdr.name.ini
use.end.tangent <- use.end.tangent.ini
use.ref.tangent <- use.ref.tangent.ini








#    # restore original values of n and pto
#    n <- n0
#    pto <- win.Lapp/n     # apparent length (m) of a single point
                                                    

}     # end of for(ii in 1:n.selected.waves)


## convert output to dataframe, change the format of the columns and add temperature and factors contained in the filename:
      output.Ka <- data.frame(output.Ka)
      output.Ka[,1] <- as.character(output.Ka[,1])
      output.Ka[ ,c(2:ncol(output.Ka))] <- data.frame(lapply(output.Ka[ ,c(2:ncol(output.Ka))], as.character), stringsAsFactors=FALSE)
      output.Ka[ ,c(2:ncol(output.Ka))] <- data.frame(lapply(output.Ka[ ,c(2:ncol(output.Ka))], as.numeric), stringsAsFactors=FALSE)
      colnames(output.Ka) <- c("file.name", "Ka", "vwc.topp", "refpoint", "endpoint", "win.Lapp", "win.tt", "pto_m", "pto_ps") 
       
           


                       # export the output to a .txt file
                       write_delim(output.Ka, paste(probe.name, ".Ka.txt", sep=''), delim='\t', append=T, col_names=T)

               

##}
