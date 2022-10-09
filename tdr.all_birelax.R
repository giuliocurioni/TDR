##tdr.all_birelax <- function(...) {
# Giulio Curioni
# 04 July 2018
 
     
## TDR WAVEFORM ANALYSIS FOR PERMITTIVITY, CONDUCTIVITY, INVERSION THROUGH AN FFT
   # the analysis is based on TDR waveforms taken with Campbell Scientific TDR100 and Campbell Scientific 3-wire probes
   # it should work with any probe provided that there is a drop in the reflection coefficient in the probe head
   
     # although not compulsory, waveform should be saved with the following format, keep '-' as the separator EXCEPT temperature (if present). Use '_' to contain temperature, if temperature is measured. See examples below for conductivity measurements used for calibration. The file formats can be different, depending on which software was used to take measurements (e.g. '.mat', '.DAT'):
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
        
## ## IMPORTANT! FOR CORRECT RESULTS EACH MEASUREMENT MUST CONSIST OF BOTH PERMITTIVITY AND CONDUCTIVITY WAVEFORMS (except if doing an inversion). The files can be in different folders but must be IN THE SAME ORDER. The script asks to select a number of waveforms to analyse. First select permittivty measurements, then select conductivity (BEC) measurements (for inversion only one waveform is required). 
# NOTE: if conductivity measurements are not available select twice permittivity measurements but do not look at the conductivity output. 

#check ### before running the script (note: search ### ### for important points for TDR inversion)
           
         source("C:/R.wd/loadall.R")        # load packages and functions
         loadall()
         library(optimization)
         library(GenSA)
         graphics.off()
         options(warn=-1)       # suppress warnings
         options(windowsTimeout = c(1000, 1000)) # see help for windows {grDevices}, try using a longer delay when plotting MANY waveforms
         options(device="windows")
         #options(digits=3)
         

  ### set number of data points
  n0 <- 2048 
  name_standard_format <- TRUE # it is assumed that the file names are saved with the standard format described above

# inversion
### ### 0) is inversion required? 
inversion <- TRUE    # TRUE does the inversion
only.fft.waveform <- FALSE  # use only FFT waveform (TRUE) or select separate waveform for BEC (FALSE)?
                            # note: although it is more convenient to analyse the FFT waveforms only by using the medium distance waveforms, BEC is overestimated typically by 5-10% compared to using long-dist. waveforms.
                            # USE THIS OPTION (i.e. TRUE) ONLY WHEN BEC IS NOT PARTICULARLY IMPORTANT. 
                            # ALSO, do not look at cond_DC optimised values for estimations of DC cond. This is likely overestimated, more than BEC, because the optimisation tends to settle at the upper limit allowed (cond_DC.up)
                            # FFT waveform: taken with longer win.Lapp to see multiple reflections 

## This script imports data from the tables created with waves.to.table.R, which must be run first (comment if not needed). 
  # Run waves.to.table.R first to be able to run this script.
  source("C:/R.wd/scripts/tdr/waves.to.table.R"); waves.to.table()
            
### ### 1) ESSENTIAL! import calibration and measurement parameters. Note: names are case sensitive! IMPORTANT: SELECT ONLY MEASUREMENTS TAKEN WITH THE SAME PROBE!
       probe.name <- "lab1.2"    ## ## choose probe (it must match one of the existing probe.ID in tdr.probe.par.txt)
       mux.level <- "mux0"       ## ## choose the multiplexer level used for the measurements. Allowed options are: mux0, mux1, mux2   
       tdr.name <- "ATU.lab"     ## ## choose the TDR unit used for the measurements. Allowed options are: TDR19200, Dan.PhD, Giulio.PhD, ATU.lab, ATU.red, ATU.blue 
       win.Lapp <- 2.6           ## ## input chosen length parameter (this is important in case something different from calibration is used)  

### 2) is the waveform a shorted waveform in air? YES = TRUE;  NO = FALSE.  Choose accordingly. 
       shorted <- F  
       
## 3) are ~ horizontal tangent lines (not perfectly horizontal lines) necessary (e.g. in dry soil)? Choose accordingly.
# NOTE: this is set to FALSE for shorted measurements.
       use.end.tangent <- F      # for the end point (at the end of the rods)
       use.ref.tangent <- F      # for the reference point (in the probe head) [FALSE for Campbell Scientific probes or any probes with a reflection coeff drop in the head]
      
       if(shorted == T) {use.end.tangent <- FALSE}

### 4) is an interactive solution necessary to better plot the tangents? YES = TRUE;  NO = FALSE.
       interactive.sol <- F

          ## NOTE ON TANGENTS:
           # tangent 1 is the rising tangent near the reference point
           # tangent 2 is the rising tangent near the end point
           # tangent 3 is the horizontal or nearly horizontal tangent line near the reference point
           # tangent 4 is the horizontal or nearly horizontal tangent line near the end point
 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
      
      # set interactive.tan <- FALSE by default (no interactive solution)
      interactive.tan1 <- F
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
      if(3 %in% substr.matrix == T) {interactive.tan3 <- T}
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
tdr.name.ini <- tdr.name
use.end.tangent.ini <- use.end.tangent
use.ref.tangent.ini <- use.ref.tangent


## import probe parameters:
    tdr.probe.par <- read_delim("C:/R.wd/scripts/tdr/tdr.probe.par.txt", delim='\t', col_names = T, na = "NA", skip = 0, n_max = Inf) 
    tdr.probe.par <- mutate_at(tdr.probe.par, c(4:ncol(tdr.probe.par)), as.numeric)
    #tdr.probe.par <- format(tdr.probe.par, digits=8)
    tdr.unit <- tdr.probe.par$tdr.unit
    probe.ID <- tdr.probe.par$probe.ID
    mux <- tdr.probe.par$mux
#    tdr.probe.par <- data.frame(lapply(tdr.probe.par[,4:length(tdr.probe.par)], as.character), stringsAsFactors=FALSE)
#    tdr.probe.par <- data.frame(lapply(tdr.probe.par[,1:length(tdr.probe.par)], as.numeric), stringsAsFactors=FALSE)
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
#    win.Lapp <- tdr.probe.par$win.Lapp[n.probe.ID]          # apparent length (m) of the window (entire TDR plot = length parameter in PCTDR)                                                    
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
    a <- tdr.probe.par$wire.diam[n.probe.ID]          # outer diameter of INNER conductor (m) [i.e. rod diameter]       
    b <- tdr.probe.par$ext.rod.spacing[n.probe.ID]    # inner diameter of OUTER conductor (m) 
    L <- tdr.probe.par$L[n.probe.ID]       # calibrated L after fft (TDR inversion)
    z <- tdr.probe.par$z[n.probe.ID]       # calibrated z after fft (TDR inversion)
    alpha <- tdr.probe.par$alpha[n.probe.ID]          # calibrated alpha after fft (TDR inversion)
    t0 <- tdr.probe.par$t0[n.probe.ID]                # calibrated t0 after fft (TDR inversion) 
    
    # import input functions from FFT calibration in air
    input.fun <- read_delim("C:/R.wd/scripts/tdr/tdr.input.fun.txt", delim = "\t", col_names = F, na = "NA", skip = 0, n_max = Inf)
    input.fun.header <- input.fun[1:3,-1]
    input.fun.values <- input.fun[-1:-3,-1]
    input.fun.values <- data.frame(lapply(input.fun.values[,], as.numeric), stringsAsFactors=FALSE)
    
    input.wave <- input.fun.values[,which(input.fun.header[1,] == tdr.name & input.fun.header[2,] == probe.name & input.fun.header[3,]== mux.level)]
    
    # if input.wave does not exist for the selected probe try importing the air measurements taken with another probe
    ### ### IMPORTANT: this works only for 15cm probes!
    if(length(input.wave) == 0 & mux.level == "mux0") {
                          input.wave <- input.fun.values[,which(input.fun.header[2,] == "lab1.2" & input.fun.header[3,] == "mux0")]
    # also use the calibrated values for the other parameters
    L <- tdr.probe.par$L[which(tdr.probe.par$probe.ID == "lab1.2" & tdr.probe.par$mux == "mux0")]     # calibrated L after fft (TDR inversion)
    z <- tdr.probe.par$z[which(tdr.probe.par$probe.ID == "lab1.2" & tdr.probe.par$mux == "mux0")]     # calibrated z after fft (TDR inversion)
    alpha <- tdr.probe.par$alpha[which(tdr.probe.par$probe.ID == "lab1.2" & tdr.probe.par$mux == "mux0")]  # calibrated alpha after fft (TDR inversion)
    t0 <- tdr.probe.par$t0[which(tdr.probe.par$probe.ID == "lab1.2" & tdr.probe.par$mux == "mux0")]        # calibrated t0 after fft (TDR inversion)
    }
    if(length(input.wave) == 0 & mux.level == "mux2") {
                              input.wave <- input.fun.values[,which(input.fun.header[2,] == "lab1.2" & input.fun.header[3,] == "mux2")]
    # also use the calibrated values for the other parameters
    L <- tdr.probe.par$L[which(tdr.probe.par$probe.ID == "lab1.2" & tdr.probe.par$mux == "mux2")]     # calibrated L after fft (TDR inversion)
    z <- tdr.probe.par$z[which(tdr.probe.par$probe.ID == "lab1.2" & tdr.probe.par$mux == "mux2")]     # calibrated z after fft (TDR inversion)
    alpha <- tdr.probe.par$alpha[which(tdr.probe.par$probe.ID == "lab1.2" & tdr.probe.par$mux == "mux2")]  # calibrated alpha after fft (TDR inversion)
    t0 <- tdr.probe.par$t0[which(tdr.probe.par$probe.ID == "lab1.2" & tdr.probe.par$mux == "mux2")]        # calibrated t0 after fft (TDR inversion)
    }
 
    
    
### ### if using custom parameters (e.g. window length for TDR inversion) manually modify and uncomment the parameters (uncomment only the ones to be modified):
    # note: the start.plot value does not make a difference in the calculation of permittivity so it is not added here
       
#        win.Lapp <- 5.2      
#        start.plot <- 7.6
#        n0 <- 2048 
#        Lcal <- 0.14969	
#        L0 <- 0.04334
#        L <-  0.1498				
#        z <- 	0.3381   
    
    n <- n0
    pto <- win.Lapp/n     # apparent length (m) of a single point

## 5) import waveform(s) - uncomment only the method used. IMPORTANT: SELECT ONLY MEASUREMENTS TAKEN WITH THE SAME PROBE! 
#          path.name <- data.frame(choose.files(caption = "Select file with PERMITTIVITY measurements"))       # select permittivity (Ka) waveforms
#          path.name[,1] <- as.character(path.name[,1])
#          n.waves <- length(path.name[,1])     # number of selected waveforms to be analysed
#              
##              print("Select conductivity waveforms corresponding to permittivity waveforms", quote=F)
##              print("IMPORTANT: to each permittivity waveform MUST correspond a BEC waveform. Ka and BEC waveforms can be in separate folders but MUST be in the same order", quote=F)
#          path.name.bec <- data.frame(choose.files(caption = "Select file with CONDUCTIVITY measurements"))   # select conductivity (bec) waveforms
#          path.name.bec[,1] <- as.character(path.name.bec[,1]) 
#message("Run waves.to.table.R first!")

          # use readr package for faster import:
          perm <- read_delim("waves.Ka.txt", delim = "\t", col_names = F,   
                 locale = locale(date_names = "en", date_format = "%Y-%m-%d", time_format = "%H:%M:%S", decimal_mark = ".", tz = "UTC"),
                 na = "NA", skip = 0, n_max = Inf)
          #perm[,2:n] <- data.frame(lapply(perm[,2:n], as.numeric), stringsAsFactors=FALSE)     
          perm <- data.frame(perm)
          colnames(perm) <- c("file.name",	1:(ncol(perm)-1))        

          cond <- read_delim("waves.BEC.txt", delim = "\t", col_names = F,   
                 locale = locale(date_names = "en", date_format = "%Y-%m-%d", time_format = "%H:%M:%S", decimal_mark = ".", tz = "UTC"),
                 na = "NA", skip = 0, n_max = Inf)
          #cond[,2:n] <- data.frame(lapply(cond[,2:n], as.numeric), stringsAsFactors=FALSE) 
          cond <- data.frame(cond)
          colnames(cond) <- c("file.name",	1:(ncol(cond)-1))        
  
          #number of measurements to analyse
          n.waves <- nrow(perm)
          if(n.waves != nrow(cond)) {stop("Number of rows in permittivity and conductivity dataframes do not match.")}
          n.col <- ncol(perm)

      output.Ka <- matrix(nrow=n.waves, ncol=3)         # create file for the Ka output
      output.BEC <- matrix(nrow=n.waves, ncol=3)        # create file for the BEC output
      output.tdr <- matrix(nrow=n.waves, ncol=29)       # create file for the Jung et al. 2013 method
      temperatura <- c(1:n.waves)
      output.fft <- matrix(nrow=n.waves, ncol=38)       # create file for the TDR inversion output (fft)     

      output.Ka <- data.frame(output.Ka)         # create file for the Ka output
      output.BEC <- data.frame(output.BEC)       # create file for the BEC output
      output.tdr <- data.frame(output.tdr)       # create file for the Jung et al. 2013 method
      output.fft <- data.frame(output.fft)       # create file for the TDR inversion output (fft)
      
      colnames(output.Ka) <- c("file.name", "Ka", "vwc.topp")
      colnames(output.BEC) <- c("file.name.bec", "BEC.mS.m", "RC")
      colnames(output.tdr) <- c("file.name", "Ka", "BEC", "vwc.topp", 
                          "pto", "Lapp", "V1", "V1m", "Vf", 
                          "A.m1.der", "A.m2.der", "A.end.plus.der", "A.pulse.der", 
                          "min.ref.x", "min.ref.y", "refpoint.x", "refpoint.y", "inflection1.x", "inflection1.y", 
                          "max.start.x", "max.start.y", "startpoint.x", "startpoint.y", 
                          "min.end.x", "min.end.y", "endpoint.x", "endpoint.y", "inflection2.x", "inflection2.y")      
      colnames(output.fft) <- c(
                       "file.name",  
                       "BETA.time", "Es.time", "frel1.time", "frel2.time", "Einf.time", "cond_DC.time", "z.time", "Mdisp_time", "r2.time",
                       "BETA.freq", "Es.freq", "frel1.freq", "frel2.freq", "Einf.freq", "cond_DC.freq", "z.freq", "Mdisp_freq", "r2.freq",
                       "f.eff", "Ka_eff_time", "Ka_eff_freq", "Er_eff_time", "Ei_eff_time", "Er_eff_freq", "Ei_eff_freq",
                       "Ka_100_time", "Ka_100_freq", "Er_100_time", "Ei_100_time", "Er_100_freq", "Ei_100_freq",
                       "Ka_1000_time", "Ka_1000_freq", "Er_1000_time", "Ei_1000_time", "Er_1000_freq", "Ei_1000_freq"
                       )      
                        
      for(ii in 1:n.waves) {       # select one or multiple measurements to analyse 

      file.name <- as.character(perm$file.name[ii])
      file.name.bec <- as.character(cond$file.name[ii])
      
              wave <- perm[ii, 2:n.col] 
              wave0 <- t(wave)  # backup original waveform
              wave <- as.numeric(wave0)
              
              bec.wave <- cond[ii, 2:n.col] 
              bec.wave0 <- t(bec.wave)   # backup original waveform
              bec.wave <- as.numeric(bec.wave0)
            
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
      M <- cbind(x, wave, bec.wave)
      M <- data.frame(M)    
      colnames(M) <- c("x", "wave", "bec.wave")

      ## define important parameters
      j <- sqrt(as.complex(-1))   # or simply 0+1i  [better to use j in case I use i in a 'for' loop]
      c0 <- 2.9979e8              # speed of light (m/s)
      E0 <- 8.8542e-12            # absolute dielectric permittivity of free space (F/m)
      mu0 <- 1.2566e-6            # absolute magnetic permeability of free space (H/m)
      mu.abs <- mu0               # note: this assumes that the material is non-magnetic (H/m)
      mu <- mu.abs/mu0            # relative magnetic permeability




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

## TTA: travel time analysis (traditional method for analysing TDR waveforms using tangent lines).  See Heimovaara and Bouten, 1990; Menziani et al., 1996.
# IMPORTANT: some parts are necessary for the frequency domain analysis, DO NOT DELETE
                                                                                                                      
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
if(n <= 2048) {roll.width <- 50} else {roll.width <- 100}
#roll.width <-  100    ### ### only in case of problems, modify this value (use a smaller number for large win.Lapp; e.g. for win.Lapp=2.5 try 100, for win.Lapp=10 use 20)
half.width <- roll.width/2


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
                
        if(shorted == F) {
              #measurement is not shorted.

               ## in theory the second inflection point is the absolute max of the first derivative AFTER the first inflection point:
               der1.max.sub.end  <- der1.max[der1.max > inflection1]
               max.der1.after.inflection1 <- max(der1[der1.max.sub.end])      # this is the global max of der1 after inflection1
               index.abs.max <- which.max(der1[der1.max.sub.end])
               inflection2 <- der1.max.sub.end[index.abs.max]                 # this is the x index corresponding to the second inflection point

               } else {
               #measurement is shorted.
               #second inflection point (descending part of the waveform at the end of the probe (end point), min of the first derivative AFTER inflection1)
               inflection2 <- which.min(der1[inflection1:length(der1)]) + inflection1     # this is the x index corresponding to the second inflection point
               }

## find min of the waveform corresponding to the reference point (head of the probe)

            # take the global minimum among the local minima before inflection1 corresponding to the drop in the probe head
            n.vec.min.ref <- wave.min[wave.min < inflection1]
            n.min.ref <- which.min(wave[n.vec.min.ref])
            min.ref <- n.vec.min.ref[n.min.ref]  # index corresponding to the min in the probe head
            min.ref.value <- M$wave[min.ref]     # min value corresponding to the min in the probe head

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
     if(length(der1.min.after.inflection1) == 0) {der1.min.after.inflection1 <- NA}
     if(is.na(der1.min.after.inflection1)) {
               der1.min.after.inflection1 <- n-10
               } else { #if there isn't a min in der1 after inflection2 try using one of the last points in the waveform (e.g. n-10)
               der1.min.after.inflection1 <- der1.min.after.inflection1[1]  # take the first min of der1 after inflection1
               der1.min.after.inflection1 <- der1.min[der1.min.after.inflection1]   # take the actual index in der1 (not a index of a subset)
               }
                       
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
     upp.lim4 <- (der1.min.after.min.end - min.end) + 5    # upper limit to draw tangent (5 is an arbitrary number of data points to use to the right of der1.min.after.min.end)
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






#-------------------------------------------------------------------------------
                   


             


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
          points(x=sub.tangent1$x, y=sub.tangent1$wave, col='red', pch='.', cex=1.6)        
          points(x=sub.tangent2$x, y=sub.tangent2$wave, col='red', pch='.', cex=1.6)        

          par(new=TRUE)

          ## choose this to plot horizontal TANGENT lines (i.e. not necessarily horizontal)
      
          if(use.ref.tangent == TRUE) {
                  abline(coef(lm.tangent3)) 	# to draw horizontal TANGENT for reference point 
                  points(x=sub.tangent3$x, y=sub.tangent3$wave, col='red', pch='.', cex=1.6)
                  } else {
          plot(M$x, x.ref, type='l', pch='.', ylim=c(-1,1), axes=F, xlab=NA, ylab=NA)		# for plotting horizontal line
          } # end of if(use.ref.tangent == TRUE)
          par(new=T)
          plot(lapp, x.end, type='l', pch='.', ylim=c(-1,1),  xlab='apparent length (m)', ylab='reflection coefficient', font=1, font.lab=1, cex.lab=1.2, cex.axis=1.2)		# for plotting horizontal line
                 
          } else {
          # use ~horizontal tangent line for the end point (an option is included also for ref point).

          par(mar=c(5, 5, 4, 2) + 0.1)
          plot(M$x, M$wave, type='o', pch='.', ylim=c(-1,1), axes=F, xlab=NA, ylab=NA)
          # plot ~ vertical tangent lines
          abline(coef(lm.tangent1))
          abline(coef(lm.tangent2))
          points(x=sub.tangent1$x, y=sub.tangent1$wave, col='red', pch='.', cex=1.6)        
          points(x=sub.tangent2$x, y=sub.tangent2$wave, col='red', pch='.', cex=1.6)
   
          par(new=TRUE)

          ## choose this to plot horizontal TANGENT lines (i.e. not necessarily horizontal)
      
          if(use.ref.tangent == TRUE) {
                  abline(coef(lm.tangent3)) 	# to draw horizontal TANGENT for reference point 
                  points(x=sub.tangent3$x, y=sub.tangent3$wave, col='red', pch='.', cex=1.6)
                  } else {
          plot(M$x, x.ref, type='l', pch='.', ylim=c(-1,1), axes=F, xlab=NA, ylab=NA)		# for plotting horizontal line
          } # end of if(use.ref.tangent == TRUE)
      
          abline(coef(lm.tangent4)) 	  # to draw horizontal TANGENT for end point
          points(x=sub.tangent4$x, y=sub.tangent4$wave, col='red', pch='.', cex=1.6)
          
          par(new=TRUE)

          plot(lapp, x.ref, type='n', pch='.', ylim=c(-1,1), axes=F, xlab=NA, ylab=NA)		# leave it with type='n' (suppress plot)
          par(new=T)
          plot(lapp, x.end, type='n', pch='.', ylim=c(-1,1),  xlab='apparent length (m)', ylab='reflection coefficient', font=1, font.lab=1, cex.lab=1.2, cex.axis=1.2)		# leave it with type='n' (suppress plot)
          
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
          points(x=sub.tangent1$x, y=sub.tangent1$wave, col='red', pch='.', cex=1.6)        
          points(x=sub.tangent2$x, y=sub.tangent2$wave, col='red', pch='.', cex=1.6)        
          
          par(new=TRUE)

          ## choose this to plot horizontal TANGENT lines (i.e. not necessarily horizontal)
      
          if(use.ref.tangent == TRUE) {
                  abline(coef(lm.tangent3)) 	# to draw horizontal TANGENT for reference point 
                  points(x=sub.tangent3$x, y=sub.tangent3$wave, col='red', pch='.', cex=1.6)
                  } else {
          plot(M$x, x.ref, type='l', pch='.', ylim=c(-1,1), axes=F, xlab=NA, ylab=NA)		# for plotting horizontal line
          } # end of if(use.ref.tangent == TRUE)
          par(new=T)
          plot(M$x, x.end, type='l', pch='.', ylim=c(-1,1),  xlab='data points', ylab='reflection coefficient', font=1, font.lab=1, cex.lab=1.2, cex.axis=1.2)		# for plotting horizontal line
           
          } else {
          # use ~horizontal tangent line for the end point (an option is included also for ref point).

          par(mar=c(5, 5, 4, 2) + 0.1)
          plot(M$x, M$wave, type='o', pch='.', ylim=c(-1,1), axes=F, xlab=NA, ylab=NA)
          # plot ~ vertical tangent lines
          abline(coef(lm.tangent1))
          abline(coef(lm.tangent2))
          points(x=sub.tangent1$x, y=sub.tangent1$wave, col='red', pch='.', cex=1.6)        
          points(x=sub.tangent2$x, y=sub.tangent2$wave, col='red', pch='.', cex=1.6)        
          
          par(new=TRUE)

          ## choose this to plot horizontal TANGENT lines (i.e. not necessarily horizontal)

          if(use.ref.tangent == TRUE) {
                  abline(coef(lm.tangent3)) 	# to draw horizontal TANGENT for reference point
                  points(x=sub.tangent3$x, y=sub.tangent3$wave, col='red', pch='.', cex=1.6)
                  } else {
          plot(M$x, x.ref, type='l', pch='.', ylim=c(-1,1), axes=F, xlab=NA, ylab=NA)		# for plotting horizontal line
          } # end of if(use.ref.tangent == TRUE)

          abline(coef(lm.tangent4)) 	  # to draw horizontal TANGENT for end point
          points(x=sub.tangent4$x, y=sub.tangent4$wave, col='red', pch='.', cex=1.6)

          par(new=TRUE)
          
          plot(M$x, x.ref, type='n', pch='.', ylim=c(-1,1), axes=F, xlab=NA, ylab=NA)		# leave it with type='n' (suppress plot)
          par(new=T)
          plot(M$x, x.end, type='n', pch='.', ylim=c(-1,1),  xlab='data points', ylab='reflection coefficient', font=1, font.lab=1, cex.lab=1.2, cex.axis=1.2)		# leave it with type='n' (suppress plot)

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
      refM <- pto*refpoint
      startM <- start.plot+refM+L0
      startpoint <- (startM-start.plot)/pto
      startpoint <- round(startpoint, digits=6)
      endM <- start.plot+endpoint*pto
      Lapp <- endM-startM
      Ka <- (Lapp/Lcal)^2          # I removed <<-
      Ka <- round(Ka, digits=2)             # round value to 2 digital digits
      
      # calculate volumetric water content using Topp et al. 1980:
      vwc.topp <- -5.30e-2 + 2.92e-2 * Ka - 5.5e-4 * Ka^2 + 4.3e-6 * Ka^3
      vwc.topp <- vwc.topp*100   # convert values to percentages
      vwc.topp <- round(vwc.topp, digits=1)       # round value to 1 decimal digits      
      
      # define all points/segments/areas of interest 
      min.ref.x <- min.ref                  # X value of min before refpoint
      min.ref.y <- M$wave[min.ref]          # Y value of min before refpoint
      refpoint.x <- refpoint                # X value of refpoint
      refpoint.y <- M$wave[refpoint]        # Y value of refpoint
      inflection1.x <- inflection1          # X value of inflection1
      inflection1.y <- M$wave[inflection1]  # Y value of inflection2
      max.start.x <- peak.start             # X value of max close to startpoint (bump in the head)
      max.start.y <- M$wave[peak.start]     # Y value of max close to startpoint (bump in the head)
      startpoint.x <- startpoint            # X value of startpoint
      startpoint.y <- M$wave[startpoint]    # Y value of startpoint
      min.end.x <- min.end                  # X value of min before endpoint
      min.end.y <- M$wave[min.end]          # Y value of min before endpoint
      endpoint.x <- endpoint                # X value of endpoint
      endpoint.y <- M$wave[endpoint]        # Y value of endpoint
      inflection2.x <- inflection2          # X value of inflection2
      inflection2.y <- M$wave[inflection2]  # Y value of inflection2
       
      max.min.sub <- M[which(M$x >= max.start.x & M$x <= min.end.x), ]       # select points between bump near the start and min before end point
      start.end.sub <-  M[which(M$x >= startpoint.x & M$x <= endpoint.x), ]  # select points between start and end points
      end.plus.sub <- M[which(M$x > endpoint.x), ]                           # select points after end point
      pulse.sub <- M[which(M$x >= endpoint.x & M$x <= inflection2.x), ]      # select points between end point and inflection2
      
      max.min.sub$wave <- max.min.sub$wave + 1                           # take absolute values (distance from reflection coeff = -1)
      start.end.sub$wave <- start.end.sub$wave + 1                       # take absolute values (distance from reflection coeff = -1)
      end.plus.sub$wave <- end.plus.sub$wave + 1                         # take absolute values (distance from reflection coeff = -1)
      pulse.sub$wave <- pulse.sub$wave + 1                               # take absolute values (distance from reflection coeff = -1)
      
      if(nrow(max.min.sub) != 0){                                        # this is needed in case the max at the start is not present (e.g. in air for short probes)
            der.i <- diff(max.min.sub$wave)                                # derivative
              x.der1 <- 1:(length(max.min.sub$wave)-1)                     # the first derivative has one point less 
              der1.spline <- smooth.spline(x=x.der1, y=der.i)              # smoothing
              der1.max.min <- abs(der1.spline$y)
              } else { der1.max.min <- NA }                                # save smoothed first derivative as a vector 
      der.i <- diff(start.end.sub$wave)
        x.der1 <- 1:(length(start.end.sub$wave)-1)                   # the first derivative has one point less 
        der1.spline <- smooth.spline(x=x.der1, y=der.i)              # smoothing
        der1.start.end <- abs(der1.spline$y)                         # save smoothed first derivative as a vector 
      der.i <- diff(end.plus.sub$wave)
        x.der1 <- 1:(length(end.plus.sub$wave)-1)                    # the first derivative has one point less 
        der1.spline <- smooth.spline(x=x.der1, y=der.i)              # smoothing
        der1.end.plus <- abs(der1.spline$y)                          # save smoothed first derivative as a vector 
      der.i <- diff(pulse.sub$wave)
        x.der1 <- 1:(length(pulse.sub$wave)-1)                       # the first derivative has one point less 
        der1.spline <- smooth.spline(x=x.der1, y=der.i)              # smoothing
        der1.pulse.sub <- abs(der1.spline$y)                         # save smoothed first derivative as a vector 

      A.m1.der <- sum(der1.max.min)                                      # area derivative below waveform between bump near the start and min near the end
      A.m2.der <- sum(der1.start.end)                                    # area derivative below waveform between startpoint and endpoint
      A.end.plus.der <- sum(der1.end.plus)                               # area derivative below waveform after endpoint
      A.pulse.der <- sum(der1.pulse.sub)                                 # area derivative below waveform between endpoint and inflection2
      
      V1 <- abs(max.start.y - min.end.y)
      V1m <- abs(startpoint.y - endpoint.y)
           
      ## import conductivity measurement to find Vf (see jung et al., 2013)
      cond.wave <- M$bec.wave[(n-(long.dist.points-1)):n]     # take the last long.dist.points data points
      RC <- mean(cond.wave)			# measured reflection coefficient at long distance      [not normalised]
      RC <- round(RC, digits=8)
      Vf <- RC + 1               
                     

      ## ## add points and print apparent permittivity/VWC by Topp on the plot
      
      ## ## if critical points and Ka and VWC are not wanted on the plot set plot.points <- F  (or comment unwanted bits below)
      
plot.points <- TRUE
      
if(plot.points == T) {
      
      # add critical points to the graph (useful for example to choose a suitable range of data points for rising tangents) 
#      if(use.end.tangent == F) {
               
      points(x=c(min.ref.x,
            refpoint.x, 
            inflection1.x,
            startpoint.x, 
            min.end.x,
            endpoint.x, 
            inflection2.x
            ), 
            y = c(min.ref.y,
            refpoint.y, 
            inflection1.y,
            startpoint.y, 
            min.end.y,
            endpoint.y, 
            inflection2.y
            ), pch=c(3,1,4,1,3,1,4), col=c('red', 'blue', 'red', 'blue', 'red', 'blue', 'red'), cex=1.2)
            
            ## add labels to some points
            text(x=c(refpoint.x, 
            startpoint.x, 
            endpoint.x), 
            y = c(refpoint.y,
            startpoint.y, 
            endpoint.y), labels=c("refpoint", "startpoint", "endpoint"), pos = 4, offset = 0.5, cex=1.2)
              
            ## add title with apparent permittivity and VWC calculated by Topp et al., 1980
           # print.Ka <- paste("Apparent permittivity", " = ", as.character(Ka), sep=" ")
           # print.vwc.topp <- paste("Topp VWC (%):", as.character(vwc.topp), sep=" ")
           # mtext(text=print.Ka, side = 3, line = 2, cex=1.2, col='black')
           # mtext(text=print.vwc.topp, side = 3, line = 0.5, cex=0.8, col='black')
           # mtext(text=bquote(~K[a] == .(Ka)), side = 3, line = 2, cex=1.4, col='black')
            mtext(text=bquote(paste(~K[a] == .(Ka), "     ", ~theta[Topp] == .(vwc.topp), "%", sep=' ')), side = 3, line = 0, cex=1.7, col='black')

} # end of if(plot.point == T)

      savePlot(filename = paste(file.name, "_TTA.png", sep=""), type = c("png"))    # save the plot (TTA = travel time analysis)

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

                            


                                                     
 
 
 
 
 
 
 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
      
      ## analyse conductivity waveform(s)
      #cond.wave <- M$bec.wave[(n-(long.dist.points-1)):n]     # take the last long.dist.points data points
      #RC <- mean(cond.wave)			# measured reflection coefficient at long distance      # I removed <<-
      #RC <- round(RC, digits=8)
      Rcorr <- ((2*(RC+1))/(Ropen+1))-1	  # corrected reflection coefficient
      RL <- 50*((1+Rcorr)/(1-Rcorr))	    # load resistance       # I removed <<-
      BEC <- Kp/(RL-(L.cable*Rc+R0))	    # bulk electrical conductivity in S/m     # I removed <<-
      BEC <- round(BEC, digits=4)         
      
#      # plot
#      windows(w=12, h=7)
#      par(mar=c(5, 5, 4, 2) + 0.1)
#      plot(x, M$bec.wave, type='l', lwd=1, ylim=c(-1,1), xlab='data points', ylab='reflection coefficient', font=1, font.lab=1, cex.lab=1.2, cex.axis=1.2)
#                  ## add title with BEC and long.dist.RC
##            print.RC <- paste("long dist. reflection coeff:", as.character(RC), sep=" ")
##            print.BEC <- paste("BEC (mS/m):", as.character(BEC*1000), sep=" ")
##            mtext(text=print.RC, side = 3, line = 2, cex=1.2, col='black')
##            mtext(text=print.BEC, side = 3, line = 0.5, cex=1.2, col='black')
#            mtext(text=bquote(paste(BEC == .(BEC*1000), " mS/m", "     ", ~V[f] == .(RC), sep=' ')), side = 3, line = 0, cex=1.7, col='black')
#      
#            print(c("long distance reflection coefficient:", RC), quote=F);
#            print(c("bulk electrical conductivity (mS/m):", BEC*1000), quote=F);
#      
#      savePlot(filename = paste(file.name.bec, "_BEC.png", sep=""), type = c("png"))    # save plot
      
                             # save the output in an object (file.name.bec, BEC, long. dist refl. coeff)
                             output.BEC[ii,] <- c(file.name.bec, BEC*1000, RC)  
      
                             # graphics.off()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #





























                                      









# *** #          
if(inversion == TRUE) {
            
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

## waveform modelling (inversion). See Jones and Or, 2004; Heimovaara, 1994; Heimovaara et al., 1996; Heimovaara, 2001; Huisman et al., 2004; Liu, 2007 (thesis)

## ## use the ORIGINAL waveform (not interpolated in the time domain) 
     
     n <- n0
     pto <- win.Lapp/n         
     x <- 1:n              
#     wave <- smooth(wave0[1:n])   # smooth initial waveform to remove some roughness and stop at n data points (just in case) [note: this smoothing does not do much or nothing at all]
    
     
## ## find sampling frequency, total time, frequency interval, time interval (between points)
            # Note: Liu (2007) used Lapp not 2Lapp. Jones and Or (2004) made a mistake in table 1 (sampling freq is 10 times what they report). Heimovaara (1994) reports Fs values of 149 GHz, comparable with the values calculated here. [note: the procedure here is correct]
                                   
      tot.time <- (2*win.Lapp)/c0       # or:  (2*pto*n)/c0      # total sampling time (s) = (two-way distance)/velocity
      delta.f <- 1/tot.time             # frequency interval (Hz) between points, see Ning Liu (2007) thesis
      Fs <- delta.f * n                 # sampling frequency (Hz)
      delta.t <- tot.time/n             # time interval (s) between points [also, delta.t = 1/Fs, that is why the inverse of the rising time is the effective frequency]     

      if(new.n != n0) {
            # since refpoint and inflection1 are needed, convert them back to the right index on the original waveform (e.g. 2048 data points)
            refpoint <- (refpoint*n0)/new.n   
            inflection1 <- (inflection1*n0)/new.n
            ## create a data frame with X = n. data points, Y = reflection coefficient
            M <- cbind(x, wave, bec.wave)
            M <- data.frame(M)    
            colnames(M) <- c("x", "wave", "bec.wave")   
            }

## normalise waveform to have first point = 0  (Heimovaara et al., 1996) and (?) scale it so that an open measurement in air reaches +1 at long distance (similar to BEC measurements)        
            wave <- wave-wave[1]
            
            
            
#tryCatch(input.wave <- read_delim(paste(probe.name, ".input.wave.txt", sep=""), delim="\t", col_names=F, skip = 0, n_max = Inf), error = function(err) {
#                            stop("Please run tdr.fft.input.fun_cal.R first.")
#                            }
#                      )   
# input.wave <- input.wave[,ncol(input.wave)]    # take the mean of the measured input waves if more than one measurement was taken in air
# alpha <- read_delim(paste(probe.name, ".alpha.txt", sep=""), delim="\t", col_names=F, skip = 0, n_max = Inf)
# alpha <- as.numeric(alpha[nrow(alpha), 2])
# #t0 <- alpha[nrow(t0), 3]
                
 

#windows()

## ## 1) generate simulated input function (use the slope of the rising bit in the head of the probe)

      erf <- function(x) {2 * pnorm(x * sqrt(2)) - 1}   # error function definition, see examples in ?pnorm documentation

#      ## parametric expression of the input function (Heimovaara 2001 and Liu, 2007):
#
#            # follow terminology given in Liu, 2007
#            # X1 <- initial          # refl coeff before step voltage, take first value of the waveform
#            # X2 <- 1                # final value of input waveform (Heimovaara et al 1996)
#            # wave <- wave/X2        # normalise by last point (Heimovaara et al., 1996)   # it makes things worse...
#
#
#      # use 0.5 for Heimovaara 2001 (see Liu 2007)
#      # use Ropen instead of X2 = 1 and mean(wave[1:10]) as X1
      aa <- (mean(wave[1:10])+Ropen)/2   # 0.5              # if using Liu 2007: (X1 + X2)/2          #  Liu didn't divide by 2 
      bb <- (Ropen-mean(wave[1:10]))/2   # 0.5              # if using Liu 2007: (X2 - X1)/2          #  Liu didn't divide by 2 

      # find alpha (Heimovaara 2001, Jones and Or 2004, Liu 2007)
      
     ## recalculate tangent 1 in case n has been changed
     
     if(new.n != n0) {
     ### TANGENT 1: rising tangent 1 (i.e. near reference point in probe head)
     # use an empirical range of data points to find low.lim and upp.lim values:
     # NOTE: this is not necessarily the same as the TTA procedure, but here only initial values are needed before the optimisation so it does not matter too much. 
  
     ## manually modify these limits if tangent for input function is plotted wrong
      low.lim1 <- 2
      upp.lim1 <- 4  
     
     sub.tangent1 <- subset(M, M$x >= inflection1 - low.lim1  &  M$x <= inflection1 + upp.lim1)
     sub.tangent1 <- data.frame(sub.tangent1)
     
                            # run linear model (i.e. tangent)
                            # rising tangent 1 across inflection1
                            tryCatch(lm.tangent1 <- lm(sub.tangent1$wave ~ sub.tangent1$x, data=sub.tangent1), error = function(err) {
                            print("Error in lm.tangent1. Try modifying low.lim1 and upp.lim1 and check inflection1 value.")
                            }
                      )        
     } 
      
     
      tangente <- lm.tangent1
      intercept.tangente <- as.numeric(coef(tangente)[1])
      slope.tangente <- as.numeric(coef(tangente)[2])
      # crossing point between horizontal line at the top of input function and rising tangent (y=mx+q):
      end.input.fun <- (Ropen-intercept.tangente)/slope.tangente     # use Ropen instead of X2 
      end.input.fun <- round(end.input.fun, digits=0)
      start.input.fun <- refpoint
      start.input.fun <- min.ref        # ?
      end.input.fun <- peak.start       # ?
      t.rise.tot <- end.input.fun - start.input.fun

      # take rise time from 10% to 90% of the peak values (see Liu 2007, Robinson et al., 2003a etc.)
      ten.percent <- t.rise.tot*0.1
      start.input.fun <- round(start.input.fun+ten.percent, digits=0)
      end.input.fun <- round(end.input.fun-ten.percent, digits=0)
      t.rise <- end.input.fun - start.input.fun
      #alpha.ini <- 1/t.rise        ## initial value
      # t0 is approx the start of the rising bit
      y.value.t0 <- (Ropen-min.ref.value)/2       # use Ropen instead of X2
      t0.ini <- (y.value.t0-intercept.tangente)/slope.tangente   ## initial value

      xx <- 1:n
      xx <- data.frame(xx)        # this is needed to perform apply

      # create the input function using estimated initial values of alpha and t0
      # see Heimovaara 2001, Jones and Or 2004, Huisman et al, 2004
      #input.fun <- function(x) { (1 + erf(alpha.ini*(x-t0.ini)))/2 }          # see Heimovaara 2001, or simply use aa=bb=0.5
      input.func.ini <- function(x) { (aa + bb * erf(alpha*(x-t0.ini))) }  # see Liu 2007 (thesis) it's basically the same as Heimovaara 2001 if aa=bb=0.5
#
      ## this is the input function calculated with initial parameters alpha.ini and t0.ini
      input.wave.ini <- apply(X=xx, MARGIN=1, FUN=input.func.ini)      # input function waveform

x.minim <- c(refpoint:end.input.fun)            # use only initial portion of waveform for minimisation to find t0
tangente.minim.fun <- slope.tangente*x.minim + intercept.tangente     # tangent to be used for optimisation
               

min.input.fun <- function(r){
# simplex optimisation to find t0, see Heimovaara 2001
      #alpha <- r[1]
      t0 <- r[1]
      zmin <- sqrt( sum( ((aa + bb * erf(alpha*(x.minim-t0))) - tangente.minim.fun)^2 ) ) # quantity to be minimised
      }

  # double check lower and upper limits (note: t0 is actually a number of data points, not time)
	optim.input.fun <- optimise(min.input.fun, lower=startpoint-0.2*startpoint, upper=startpoint+0.2*startpoint) ##, method="L-BFGS-B", lower=c(0,0), upper=c(1,1000))
  if(mux.level == "mux0") {
               t0 <- optim.input.fun$minimum[1] + (startpoint-refpoint) #/2    ## because t0 is optimised against rising tangent in the head it is underestimated. Add the distance from ref point and start point (empirically it seems better to add half of this distance. Look at armonics Re(s11))
               }
  
  if(mux.level == "mux2") {
               t0 <- optim.input.fun$minimum[1] #+ (startpoint-refpoint)/2    ## because t0 is optimised against rising tangent in the head it is underestimated. Add the distance from ref point and start point (empirically it seems better to add half of this distance. Look at armonics Re(s11))
               }
  
  # calculate input.wave based on optimised t0
  input.func <- function(x) { (aa + bb * erf(alpha*(x-t0))) }
  input.wave <- apply(X=xx, MARGIN=1, FUN=input.func)

        ## ## add head reflections to input function (it decreases accuracy in time domain optim...)
#        input.wave[1:startpoint] <- wave[1:startpoint]
#        plot(input.wave, t='l', ylim=c(-1,1))
#        # windows()

             # shift waveform to have initial point = zero (Heimovaara et al., 1996)   
             input.wave <- input.wave - input.wave[1]                


# read.delim(paste(probe.name, ".input.wave.txt", sep=""), header=F)
# input.wave <- input.wave[,ncol(input.wave)]    # take the mean of the measured input waves if more than one measurement was taken in air
# alpha <- read.delim(paste(probe.name, ".alpha.txt", sep=""), header=F)
# alpha <- alpha[nrow(alpha), 2]
# t0 <- alpha[nrow(t0), 3]

 

      ### ### find effective frequency from time rise at the END of the probe (Robinson et al., 2003a)
      tangente2 <- lm.tangent2
      intercept.tangente2 <- as.numeric(coef(tangente2)[1])
      slope.tangente2 <- as.numeric(coef(tangente2)[2])

      endpeak <- der1.min[der1.min > inflection2][1]  # take the first min of first derivative after inflection2
      if(is.na(endpeak)) {endpeak <- n-10}  #if there isn't a min in der1 after inflection2 try using one of the last points in the waveform (e.g. n-10)
     ### END PEAK TANGENT, required to find f.eff
      low.lim.endpeak <- 4
      upp.lim.endpeak <- 4
     
      sub.peak.tangent <- subset(M, M$x >= endpeak - low.lim.endpeak  &  M$x <= endpeak + upp.lim.endpeak)
      sub.peak.tangent <- data.frame(sub.peak.tangent)
               
      # run linear model (i.e. tangent)      
      tryCatch(lm.endpeak.tangent <- lm(sub.peak.tangent$wave ~ sub.peak.tangent$x, data=sub.peak.tangent), error = function(err) {
                            print("Error in lm.endpeak.tangent. Try modifying low.lim.endpeak and upp.lim.endpeak and check endpeak value. It might help to change the value of rollwidth.")
                            }
                      ) 
                       
      endpeak.tan <- lm.endpeak.tangent
      intercept.endpeak.tan <- as.numeric(coef(endpeak.tan)[1])
      slope.endpeak.tan <- as.numeric(coef(endpeak.tan)[2])
      
      # crossing point between tangent2 and horizontal line at the peak after inflection2 (y=mx+q -> x = (q2-q1)/(m1-m2) ) [NOTE: this might not exist in very attenuated waveforms!]:
      endpeak <- (intercept.endpeak.tan - intercept.tangente2) / (slope.tangente2 - slope.endpeak.tan)     ## ##
      t.rise.tot2 <- (endpeak  - endpoint)*delta.t    #use endpoint as the start of the rising bit

      # take rise time from 10% to 90% of the peak values (see Liu 2007, Robinson et al., 2003a etc.)
      ten.percent2 <- t.rise.tot2*0.1
      start2 <- round(endpoint+ten.percent2, digits=0)
      end2 <- round(endpeak-ten.percent2, digits=0)
      t.rise2 <- end2 - start2
      alpha2 <- 1/t.rise2        


 ## find the actual rise time in seconds (Robinson et al., 2003a) and the effective frequency AT THE END OF THE PROBE
    # note: f.eff represents the frequency where the MAJORITY of the power is contained  (Robinson et al., 2003a), the max frequency content seems to be between 2 and 3 times higher than f.eff
 rise.time <- delta.t/alpha2      # rise time in seconds [this is correct]
 f.eff <- 0.35/rise.time          # effective frequency (Hz) at the end of the probe (Robinson et al., 2003a)          
 
 if(f.eff < 50e6) { f.eff <- 200e6 }   ## ## in case f.eff is too small or even negative there is an error in drawing tangents. 
                                       # take the assumption that f.eff = 200 MHz (if this happens it is likely because the waveform is attenuated). 
                 

## plot original waveform and simulated input function

      x <- 1:n      # convert it back to numeric
         
#      # plot waveform, tangent and simulated input function
#      windows(w=12, h=7)
#      plot(x, wave, type='l', lty=1, xlim=c(0,n), ylim=c(-1.2,1.2), xlab=NA, ylab=NA)     # measured waveform (refl coeff)
#      #par(new=T); abline(coef(tangente), xlim=c(0,n), ylim=c(-1.2,1.2), xlab=NA, ylab=NA)             # tangent
#      #par(new=T); plot(x, input.wave.ini, type='l', lty=3, xlim=c(0,n), ylim=c(-1.2,1.2), xlab=NA, ylab=NA)   # simulated input function based on initial parameters
#      par(new=T); plot(x, input.wave, type='l', lty=2, lwd=2, xlim=c(0,n), ylim=c(-1.2,1.2), xlab="index", ylab="reflection coefficient")       
#      legend("bottomright",
#            c(paste("f.eff = ", round(f.eff*10^-6, digits=0), " MHz", sep=''),
#            "measured waveform",
#            "input waveform initial alpha, t0",
#            "input waveform optimised alpha, t0"
#            ),
#            lty=c(0, 1,3,2),
#            lwd=c(0, 1,1,2)
#            #text.col=c("black", "blue", "red")
#            )            
#      points(x=c(start2,end2), y=c(wave[start2], wave[end2]), pch='X', col='blue')  
#      savePlot(filename = paste(file.name, "_measured and input waveforms.png", sep=''), type = c("png"))
          
                    



## ## create new waveforms to cut off cable and head (take only part corresponding to rods)
# it seems that it doesn't increase the accuracy
#wave1 <- wave[startpoint:n]
#l.wave1 <- length(wave1)
#n.add.points <- n - l.wave1
#add.points <- rep(wave[n], n.add.points)     # add a number of points at the end equal to the last point to reach 2048
#wave <- c(wave1, add.points)
#
#input.wave1 <- input.wave[startpoint:n]
#add.points.input <- rep(input.wave[n], n.add.points)  # add a number of points at the end equal to the last point to reach 2048
#input.wave <- c(input.wave1, add.points.input)

## ## normalisation: shift waveform to have initial point = zero (Heimovaara et al., 1996)  [note: it doesn't improve things]
#             # input.initial <- input.wave[1]
#             # input.wave <- input.wave - input.initial       # apparently it works better without normalisation!
#             # input.wave <- input.wave/X2     # normalise by last point (Heimovaara et al., 1996)

#plot(wave1, type='l', ylim=c(-1,1))
#par(new=T)
#plot(input.wave1, type='l', ylim=c(-1,1), col='red')






# windows()






## ## 2) FFT


    Nicolson <- T  # use Nicolson's method, see Ning Liu thesis (2007) and Jones and Or (2004)
                   # if FALSE, use backward difference (i.e. derivative), see Heimovaara (1994) 
                    
        Nic.wave <- matrix(ncol=1, nrow=n)
        Nic.input.wave <- matrix(ncol=1, nrow=n)
        ramp.wave <- wave[n]*(x/n)
        ramp.input.wave <- input.wave[n]*(x/n)         
        Nic.wave <- wave - ramp.wave                     # subtract ramp
        Nic.input.wave <- input.wave - ramp.input.wave   # subtract ramp
        Nic.wave[1] <- 0
        Nic.input.wave[1] <- 0

              
              ## ## zero pad (i.e. add zeros) at the end of the waveform to increase the resolution in the freq domain. 
              # zero padding interpolates in the other domain.
              
#              n.zeros <- 2048       n+n.zeros must be a multiple of 2^x, e.g. choose 2048 (do not exagerate with the zeros!)
#              Nic.wave <- c(Nic.wave, rep(0, n.zeros))     
#              Nic.input.wave <- c(Nic.input.wave, rep(0, n.zeros))   
#             
#              n <- n+n.zeros           
#              x <- 1:n              
              
        if(Nicolson == F){
                    Nic.wave <- c(diff(wave), 0)                     # names are a little confusing but keep as it is
                    Nic.input.wave <- c(diff(input.wave), 0)         # names are a little confusing but keep as it is 
                    }
              
        # fast fourier transform
        fft.wave <- fft(Nic.wave)
        fft.input.wave <- fft(Nic.input.wave)

        ## ## find frequency content
            # this is correct. It has been double checked.
            # Note: Liu used Lapp not 2Lapp. Jones and Or made a mistake in table 1 (sampling freq is 10 times what they report).
              fn <- Fs/2                  # Nyquist frequency (max freq in the fft = Fs/2). It must be higher than the freq bandwith = 1/(2*pi*(t.rise/delta.t)) see Liu thesis
              n.fn <- n/2                   

        f.min <- Fs/n                         # minimum frequency (see Jones and Or 2004) [this is equal to delta.f!]
        f.max <- fn*2                         # maximum freq (note, the actual max freq content is the Nyquist freq fn)

                                               # check http://stackoverflow.com/questions/4364823/how-to-get-frequency-from-fft-result
                                               # https://stackoverflow.com/questions/41394984/matlab-frequency-bin-of-the-positive-and-negative-frequency
                                               # this f comprises both negative and positive freq 
        #f <- seq(1e-9, f.max-delta.f, by=delta.f)   # which is equivalent to f <- (bins*Fs)/n     
        fpos <- seq(1e-9, fn-delta.f, by=delta.f)    # https://stackoverflow.com/questions/41394984/matlab-frequency-bin-of-the-positive-and-negative-frequency
        fneg <- seq(-fn, -delta.f, by=delta.f)
        f <- c(fpos, fneg)

        #bins <- 0:(n-1)
        #f <- (bins*Fs)/n        # see https://www.youtube.com/watch?v=2IkdNsGQgEM towards the end. Very clear.

        f.min.MHz <- round(f.min * 10^-6, digits=2)
        fn.MHz <- round(fn * 10^-6, digits=2)
        print(paste(c("f.min:", f.min.MHz, "MHz"), sep=" "), quote=F)
        print(paste(c("fn:", fn.MHz, "MHz"), sep=" "), quote=F)


        R <- fft.wave
        V0 <- fft.input.wave


## ## ## ## calculate s11 (see Heimovaara, 1994)
        s11 <- R/V0
        
        # find maximum frequency based on magnitude of Re(s11)
        first_below1 <- which(Re(s11) < -1)[1]
        first_above1 <- which(Re(s11) > 1)[1]
        max.useful.index <- min(first_below1,first_above1)
        max.useful.freq <- 2e9 ### ### f[max.useful.index]              # above this frequency the measured data is noise
         
        ##s11[which(f>max.useful.freq | f<max.useful.freq*-1)] <- 0+0i #s11[which(f>max.useful.freq)[1]]     ## ## this significantly increases the accuracy! particularly if applied to the simulated s11. It is OK to remove these frequencies because the actual TDR bandwidth is smaller than max.useful.freq, above it is noise
                       

        mag.s11 <- Mod(s11)      # magnitude (this is generally what is plotted)
        phase.s11 <- Arg(s11)    # phase
        rl.s11 <- 20*log10(mag.s11)     # return loss


#        # double check formulae (see Thomas et al., 2008)
#        s11m <- (10^(rl.s11/20))*cos(phase.s11) + j*(10^(rl.s11/20))*sin(phase.s11)             # this should be the same as s11 above
#        plot(f[1:n.fn], Re(s11[1:n.fn]), type='l', xlim=c(f.min, f.max))
#        par(new=T)
#        plot(f[1:n.fn], Re(s11m[1:n.fn]), col="red", type='l', xlim=c(f.min, f.max))     # plots should match







        ## check that the ifft works (note, the ifft won't match the original waveform if smoothing was applied):

        if(FALSE) {
              i.V0 <- fft(V0, inverse=T)/length(V0)
              windows()
              plot(Re(i.V0), type='l', col='red')    # this should be the same as the plot of Nic.input.wave or input.der.wave (waveform used for fft)
              title("inverted V0")
              windows()
              #par(new=T)
              if(Nicolson == T) {
                    plot(Nic.input.wave, type='l')
                    title("input waveform (Nicolson's method)")
                    } 

              windows()

              i.R <- fft(R, inverse=T)/length(R)
              plot(Re(i.R), type='l', col='red')     # this should be the same as the plot of Nic.wave or der.wave (waveform used for fft)
              title("inverted R")
              windows()
              #par(new=T)
              if(Nicolson == T) {
                    plot(Nic.wave, type='l')
                    title("original waveform (Nicolson's method)")
                    } 
        } # end of if(FALSE){}





         ## ## add here the bits to take the TSC (triple short calibration)





## ## ## ## simulate s11 (s11.sim)
                      
### CHOOSE BANDWIDTH FOR OPTIMISATION IN THE FREQ DOMAIN
## use s11 only up to a certain frequency, modify freq value as wanted. Below this frequency the values are reasonable, above it's mostly noise
max.freq <- f[max.useful.index] #f.eff *2.5 ### ### take the first point where Re(s11) is either >1 or <-1. Above it is mainly noise.

        ## ## SIMULATE E.probe ONLY UP TO max.freq... not correct! keep it but do not use it.
#        ff <- seq(f.min, max.freq, by=(max.freq-f.min)/(n0-1)) #seq(f.min, max.freq, by=(max.freq-f.min)/(n0-1))        # frequency to be used to calculate E.probe: calculate E.probe only up to max.freq! note: max.freq is approximately f.eff but is set by the user           
                   
         indice.optim <- which(f > max.freq-1) 
         n.end <- indice.optim[1]   # index   (note: this is equal to 2048 if max.freq = f[n0]
         f.end <- f[n.end]          # frequency                 



#      ## only for calibration in water:
#      # water (use cole-cole model to simulate Er and Ei, see Robinson et al., 2003a, Huisman et al., 2003a)
#            BETA <- 0               #  if BETA = 0, cole-cole model is debye
#            temp <- as.numeric(readline("input water temperature:  "))           ### ###INPUT TEMPERATURE, relevant if measuring in wet soil or water!
#            Es <- 78.54 * ( 1-4.579e-3*(temp-25) + 1.19e-5*((temp-25)^2) - 2.8e-8*((temp-25)^3) )  # Weast (1986) static permittivity of free water at temperature = temp
#            Einf <- 4.22            # infinite permittivity of free water (Robinson et al., 2003)
#            frel <- 17.4e9          # relaxation frequency of free water
#            cond_DC <- 0.01          # use 0.0001 for distilled water (better not zero)
#
#            colecole <- function(f) { Einf + ( (Es-Einf) / (1 + (j*(f/frel)^(1-BETA))) ) - j*(cond_DC/(2*pi*f*E0)) }
#
#            E.probe <- colecole(f)
#
##      # from Cassidy in Jol 2008 (page 58):
##      A <- 2*pi*f*sqrt( (mu.abs*Er*E0/2) * (sqrt(1+(Ei/Er)^2) - 1 ) )
##      B <- 2*pi*f*sqrt( (mu.abs*Er*E0/2) * (sqrt(1+(Ei/Er)^2) + 1 ) )
##      wave.num <- complex(real=A, imaginary=B)      #this is gamma (complex wave number or propagation constant)






      ### ###
      ## ## ## ## set initial values for the parameters
      BETA.ini <- 0
      BETA <- 0        ## ## IMPORTANT: USING NEGATIVE FREQUENCIES CAUSES THE OPTIMISATION TO FAIL IF BETA IS DIFFERENT FROM ZERO. DO NOT OPTIMISE BETA!
      Es1.ini <- Ka             # for water cal I can use Es, but with soil it is better to use Ka
      Es2.ini <- Ka
      Es.ini <- Ka 
      #frel1.ini <- 17.4e9       # 17.4e9 for water       
      frel2.ini <- 4e9           # (bound?) water
      Einf.ini <- 4.22 #Ka 
      Einf1.ini <- 4.22
      Einf2.ini <- 4.22
      Einf.ini <- 4.22
      cond_DC.ini <- BEC          # use BEC value, use 0.0001 for distilled water (not zero to avoid errors in cond_DC.low and cond_DC.up)
      deltaK1.ini <- Es.ini - Einf.ini
      deltaK2.ini <- Es.ini - Einf.ini                    
#      z.ini <- Zc/Zp     
#      L.ini <- Lcal          # take the value from the air-water cal for time-domain measurements as the initial value
#      Zc.ini <- Zc                
#      E.handle.ini <- 1.8    # from Liu (2007)
#      L.handle.ini <- L0     # take the value from the air-water cal for time-domain measurements as the initial value   
      alpha.ini <- alpha      # take the optimised value as initial value for next optimisation 
      t0.ini <- startpoint #t0            # take the optimised value as initial value for next optimisation 
#      # necessary if optimising also for aa and bb:
#      aa.ini <- aa
#      bb.ini <- bb
      z.ini <- z ##+ 0i
      zi.ini <- 0

# measured waveform to be optimised against:
if(Nicolson == F) {
            i.wave <- cumsum(Nic.wave)
            } else { i.wave <- Nic.wave }

## ## try many frel2 values using the other initial parameters
# note: this procedure is necessary because the optimisation (later) does not manage to optimise frel, although the final results are strongly dependent on the initial values used.
# frel2 is not much affected because there is no beta coefficient to account for the spreading of the relaxation.
initial <- c(seq(10e6, 1e9, by=10e6))      
frel1.ini <- cbind(initial, rep(NA, length(initial)))
frel1.ini <- data.frame(frel1.ini)
colnames(frel1.ini) <- c("frel1.ini", "rss")
         
for(rr in 1:length(initial)) {
      
      E.probe.ini <- Einf.ini + ( (Es.ini-Einf.ini) / (1 + (j*(f/frel1.ini[rr,1])^(1-BETA.ini))) ) + ( (Es.ini-Einf.ini) / (1 + (j*(f/frel2.ini))) ) - j*(cond_DC.ini/(2*pi*f*E0))        

  ro <- (1-z.ini*sqrt(E.probe.ini)) / (1+z.ini*sqrt(E.probe.ini))
  gammaL <- (j*2*pi*f*L*sqrt(E.probe.ini))/c0
  ## simulated s11 (Clarkson et al., 1977)
  s11.sim.ini <- (ro+exp(-2*gammaL)) / (1+ro*exp(-2*gammaL))               # open-circuit (probe)
  #s11.sim.ini <- (ro-exp(-2*gammaL)) / (1-ro*exp(-2*gammaL))              # short-circuit (cell)
  #s11.sim.ini <- ro * (1-exp(-2*gammaL)) / (1-(ro^2)*exp(-2*gammaL) )     # matched-circuit            
  ##s11.sim.ini[which(which(f>max.useful.freq | f<max.useful.freq*-1)] <- s11.sim.ini[which(which(f>max.useful.freq | f<max.useful.freq*-1)[1]]     ## ## this seems to increase the accuracy! particularly if applied to the simulated s11. It is OK to remove these frequencies because the actual TDR bandwidth is smaller than max.useful.freq, above it is noise
  
  
  ## convert to time-domain 
      R.ini <- s11.sim.ini*V0
      # R.ini[1] <- 0+0i
      
      i.s11.sim.ini <- fft(R.ini, inverse=T)/length(R.ini)
      if(Nicolson == F) {
                  i.wave.sim.ini <- cumsum(Re(i.s11.sim.ini))
                  } else {
                  i.wave.sim.ini <- Re(i.s11.sim.ini)
                  # normalise waveform to have first point = 0 
                  i.wave.sim.ini <- i.wave.sim.ini- i.wave.sim.ini[1]
                  }            
      
      # calculate R squared (see Todeschini's book, pag 124, 125)
      tss <- sum( ( i.wave - mean(i.wave) )^2 )     # tot sum of squares
      rss <- sum( ( i.wave.sim.ini - i.wave )^2 )     # residual sum of squares
      r2 <- round(1 - (rss/tss), digits=4)
      
      frel1.ini[rr,2] <- rss
      }
       
      frel1.ini <- frel1.ini$frel1.ini[which.min(frel1.ini$rss)]
 
      ## ## recalculate using best frel1 initial value 
            
     E.probe.ini <- Einf.ini + ( (Es.ini-Einf.ini) / (1 + (j*(f/frel1.ini)^(1-BETA.ini))) ) + ( (Es.ini-Einf.ini) / (1 + (j*(f/frel2.ini))) ) - j*(cond_DC.ini/(2*pi*f*E0))            

      # identify initial Er and Ei values corresponding to TDR effective frequency (i.e. f.eff)
      n.E.ini <- which(f > f.eff)[1]              
      Er.ini <- Re(E.probe.ini[n.E.ini])
      Ei.ini <- Im(E.probe.ini[n.E.ini])

      ## ## uncomment if not available in tdr.probe.par.txt
#      a <- 0.0015    # 0.00159 # 0.00318    # outer diameter of INNER conductor (m) [i.e. rod diameter]
#      b <- 0.0135    # 0.015  # 0.0305      # inner diameter of OUTER conductor (m) 

      # see Cataldo et al. book, pag 39
      LL <- Lcal
      f.up.circ <- c0/(pi*((a+b)/2)*sqrt(1*Er.ini))       # circumferential resonant (cutoff) freq. TEM assumption is valid below this freq.
      f.up.long <- c0/(2*LL*sqrt(Er.ini))                 # longitudinal resonant (cutoff) freq. TEM assumption is valid below this freq. (most important!)

#      # from Cataldo et al (2011), book, pag 13 and 38, impedance of a perfect coax transmission line: ( (60/sqrt(E.probe)) * log(b/a) ) +i*0
#      Zc <- 50 +j*0.0000         # characteristic impedance of the coax transmission line before the probe

#      ## ## method 1 (3 rod probe), see Cataldo book, pag 38
#      Zp <- (1/(4*pi)) * sqrt(mu0/E0)*log( (1-(a/b)^4) / (2*(a/b)^3) ) + j*0.0001

      ## ## method 2
#      Zp <- ( (60/sqrt(1+j*0)) * log(b/a) ) +j*0        # reference characteristic impedance (dieletric is air), see Heimovaara 1994, Liu 2007

      ## ## method 3 - Shuai et al., 2009:
#      k <- ( ((b-a)/2) + a ) / a         # ratio between spacing between two rods (from centre to centre) and rod diameter
#      num.Zp <- (4*pi*E0)/(log(((4*k^2)-1)/(4*k-1))+2*log(2*k-1))
#      den.Zp <- ((3*mu0)/(4*pi))*(log(2*k-1)+(1/3)*log((2*k+1)/(4*k-1)))
#      Zp.shuai <- 1/sqrt(num.Zp/den.Zp)       ## ## Shuai et al. 2009 say Zp=sqrt(num.Zp/den.Zp) but that seems wrong!

## plot model waveform with initial parameters (e.g. z.ini, L.ini)
  ro <- (1-z.ini*sqrt(E.probe.ini)) / (1+z.ini*sqrt(E.probe.ini))
  gammaL <- (j*2*pi*f*L*sqrt(E.probe.ini))/c0
  ## simulated s11 (Clarkson et al., 1977)
  s11.sim.ini <- (ro+exp(-2*gammaL)) / (1+ro*exp(-2*gammaL))               # open-circuit (probe)
  #s11.sim.ini <- (ro-exp(-2*gammaL)) / (1-ro*exp(-2*gammaL))              # short-circuit (cell)
  #s11.sim.ini <- ro * (1-exp(-2*gammaL)) / (1-(ro^2)*exp(-2*gammaL) )     # matched-circuit            
  ##s11.sim.ini[which(f>max.useful.freq | f<max.useful.freq*-1)] <- s11.sim.ini[which(f>max.useful.freq | f<max.useful.freq*-1)[1]]     ## ## this seems to increase the accuracy! particularly if applied to the simulated s11. It is OK to remove these frequencies because the actual TDR bandwidth is smaller than max.useful.freq, above it is noise

  
  ## convert to time-domain 
      R.ini <- s11.sim.ini*V0
      # R.ini[1] <- 0+0i
      
      i.s11.sim.ini <- fft(R.ini, inverse=T)/length(R.ini)
      if(Nicolson == F) {
                  i.wave.sim.ini <- cumsum(Re(i.s11.sim.ini))
                  } else {
                  i.wave.sim.ini <- Re(i.s11.sim.ini)
                  # normalise waveform to have first point = 0 
                  i.wave.sim.ini <- i.wave.sim.ini- i.wave.sim.ini[1]
                  }            
      
      # calculate R squared (see Todeschini's book, pag 124, 125)
      tss <- sum( ( i.wave - mean(i.wave) )^2 )     # tot sum of squares
      rss <- sum( ( i.wave.sim.ini - i.wave )^2 )     # residual sum of squares
      r2 <- round(1 - (rss/tss), digits=4)

#       ## Chen and Or (2006)
#       E.probe.ini <- Einf.ini + ( (Es.ini-Einf.ini) / (1 + (j*(f/frel.ini)^(1-BETA.ini))) ) - j*(cond_DC.ini/(2*pi*f*E0))
#       Er.ini <- Re(E.probe.ini)
#       Ei.ini <- Im(E.probe.ini)    
#       Zp.ini <- Zp/sqrt(E.probe.ini)                # note: this is the impedance of the probe NOT the characteristic impedance of the probe
#                                                     # ( 60/sqrt(E.probe.ini) ) * (b/a)
#       Zh.ini <- Zp/sqrt(E.handle.ini)               # note: this is the impedance of the handle NOT the characteristic impedance of the handle
#                                                     # ( 60/sqrt(E.handle.ini) ) * (b/a)               
#       ro.ini <- (Zp.ini-Zh.ini)/(Zp.ini+Zh.ini)
#       alpha.f.ini <- ((2*pi*f)/c0) * sqrt( (Er.ini/2) * (sqrt(1 + ((Ei.ini/Er.ini)^2)) - 1) )
#       beta.f.ini <- ((2*pi*f)/c0) * sqrt( (Er.ini/2) * (sqrt(1 + ((Ei.ini/Er.ini)^2)) + 1) )
#       H.ini <- exp(-alpha.f.ini*2*L.ini + j*beta.f.ini*2*L.ini)
#       s11.A3.ini <- (ro.ini+H.ini) / (1+ro.ini*H.ini)
#       ro.h.ini <- (Zh.ini-Zc.ini) / (Zh.ini+Zc.ini)
#       H1.ini <- exp( (j*2*pi*f*sqrt(E.handle.ini)*2*L.handle.ini) / c0 )
#      
#       s11.sim.ini <- (ro.h.ini + s11.A3.ini*H1.ini) / (1 + ro.h.ini*s11.A3.ini*H1.ini)
#      
#            R.time <- s11.sim.ini*V0
#            i.s11.sim.ini <- fft(R.time, inverse=T)/length(R.time)
#            i.wave.sim.ini <- Re(i.s11.sim.ini)              



windows(w=12, h=7)
plot(i.wave.sim.ini[1:n0], t='l', lty=2, xlim=c(0, n0), ylim=c(-1.2,1.2))
par(new=T)
plot(i.wave, t='l', xlim=c(0, n0), ylim=c(-1.2,1.2))
title("measured and modelled with initial values")
legend("topright",
             legend=c(paste("BETA.ini = ", BETA.ini, sep=''),
#             paste("deltaK1.ini = ", deltaK1.ini, sep=''),
#             paste("deltaK2.ini = ", deltaK2.ini, sep=''),
             paste("Es.ini = ", Es.ini, sep=''),
             paste("Einf.ini = ", Einf.ini, sep=''),
             paste("frel1.ini = ", frel1.ini, sep=''),
             paste("frel2.ini = ", frel2.ini, sep=''),
             paste("cond_DC.ini = ", cond_DC.ini, sep='')#,
#             paste("L.ini = ", L.ini, sep=''),
#             paste("z.ini = ", z.ini, sep='')
             ), bty='n')
savePlot(filename = paste(file.name, "_waveforms initial values.png", sep=''), type = c("png"))

             

                

                          #stop()


 


### ### set lower and upper limits

#Zc.low <- Zc.ini - 10
#Zc.up <- Zc.ini + 10
#E.handle.low <- E.handle.ini - 0.5
#E.handle.up <- E.handle.ini + 0.5
#L.low <- L.ini-0.001                  # allow +- 1mm from L.ini (Weerts et al., 2001)
#L.up <- L.ini+0.001                   # allow +- 1mm from L.ini (Weerts et al., 2001)
## allow 0.01 variation of z from z.ini. This range seems to produce better results. Zp formula (Cataldo book) is theoretical for 3-rod prob so it should be quite accurate. A variation of +- 0.01 corresponds to a variation of Zc of +- 2 ohm.
z.low <- 0.05  #0.05 according to Weerts et al., 2001
z.up <-  0.50  #0.50 according to Weerts et al., 2001
#zi.low <-  -1e-9     
#zi.up <-   1e-9 #e-8 #0.10 
#E.probe.low <- E.probe.ini - 2        # allow an error of 2 units of permittivity from TTA
#E.probe.up <- E.probe.ini + 2         # allow an error of 2 units of permittivity from TTA
#Er.low <- Er.ini - 1
#Er.up <- Er.ini + 1
#Ei.low <- Ei.ini - 1
#Ei.up <- Ei.ini + 1 
#Es1.low <- Ka -1                          # Es should always be higher than Ka  (allow for 1 unit error on Ka)
#Es1.up <- 80                              # Es in soil should be smaller than Ka in water 
#Es2.low <- Ka -1
#Es2.up <- 80
Es.low <- Ka -2                            # sometimes Ka can be overestimated, especially when using medium-length waveforms and with multiplexers. Therefore allow for some variation.
Es.up <- 80
Einf.low <- 4.20
Einf.up <- 4.24 
#Einf1.low <- 4.22                         
#Einf.up1 <- Ka -1 
#Einf2.low <- 4.22                         # Einf value in water (Robinson et al. 2003)
#Einf2.up <- Ka -1 
frel1.low <- 100e3           #  polarisation in TDR bandwidth
frel1.up <-  1e9             #  polarisation in TDR bandwidth
frel2.low <- 1e9             #  upper freq polarisation                
frel2.up <- 17.4e9           #  upper freq polarisation 
BETA.low <- 0
BETA.up <- 1
cond_DC.low <- cond_DC.ini - cond_DC.ini*0.20     # allow up to 10 % error on BEC measurements
cond_DC.up <- cond_DC.ini + cond_DC.ini*0.20      # allow up to 10 % error on BEC measurements. Because the simulation is done on a medium win.Lapp (but not very long), the simulated BEC is always overestimated (5-10%) compared with the BEC value obtained using very long distance waveforms (i.e. 500 m apparent length). Simulations tend to optimise at cond_DC.up. The value of cond_DC is overestimated, therefore if DC cond is important use BEC measured with long-dist. waveforms. 
alpha.low <- alpha - alpha*0.10          # allow up to 10 % error
alpha.up <- alpha + alpha*0.10           # allow up to 10 % error
#t0.low <- t0 - t0*0.10                   # allow up to 10 % error
#t0.up <- t0 + t0*0.10                    # allow up to 10 % error
#aa.low <- 0
#aa.up <- 1
#bb.low <- 0
#bb.up <- 1
deltaK1.low <- 0
deltaK1.up <- Es.up-Einf.low
deltaK2.low <- 0
deltaK2.up <- Es.up-Einf.low








               
                     





## ## ## ## OPTIMISATION IN THE FREQUENCY DOMAIN



min.fun.freq <- function(r){
## simplex, optimisation in the frequency domain 

#      BETA <- r[1]
      Es <- r[1]
      Einf <- Einf.ini #r[3]
      frel1 <- r[2]
      frel2 <- r[3]            
      cond_DC <- r[4] 
#      alpha <- r[5]
#      t0 <- r[7]
#      #z <- r[7]
#
#        # calculate input.wave based on calibrated alpha and t0                     
#        xx <- data.frame(1:n)
#        input.func <- function(xx) { (aa + bb * erf(alpha*(xx-t0))) }
#        input.wave <- apply(X=xx, MARGIN=1, FUN=input.func)
#        Nic.input.wave <- matrix(ncol=1, nrow=n)
#        ramp.input.wave <- input.wave[n]*(x/n)         
#        Nic.input.wave <- input.wave - ramp.input.wave   # subtract ramp
#        Nic.input.wave[1] <- 0
#        # fast fourier transform
#        fft.input.wave <- fft(Nic.input.wave)
#        V0 <- fft.input.wave
#        # recalculate measured s11 using optimised V0 for a correct comparison of sum of squares between s11 and s11.sim
#        s11 <- R/V0   
      
      #E.probe <- Einf + ( (Es-Einf) / (1 + (j*(f/frel)^(1-BETA))) ) - j*(cond_DC/(2*pi*f*E0))
      E.probe <- Einf + ( (Es-Einf) / (1 + (j*(f/frel1)^(1-BETA))) ) + ( (Es-Einf) / (1 + (j*(f/frel2))) ) - j*(cond_DC/(2*pi*f*E0))         
      
      #z <- complex(real=z, imaginary=zi)
      #z <- z.ini 
      ro <- (1-z*sqrt(E.probe)) / (1+z*sqrt(E.probe))
      gammaL <- (j*2*pi*f*L*sqrt(E.probe))/c0                  

      ## simulated s11 (Clarkson et al., 1977)
      s11.sim <- (ro+exp(-2*gammaL)) / (1+ro*exp(-2*gammaL))                    # open-circuit (probe)
      #s11.sim.ini <- (ro-exp(-2*gammaL)) / (1-ro*exp(-2*gammaL))               # short-circuit (cell)
      #s11.sim.ini <- ro * (1-exp(-2*gammaL)) / (1-(ro^2)*exp(-2*gammaL) )      # matched-circuit
      ##s11.sim[which(f>max.useful.freq | f<max.useful.freq*-1)] <- s11.sim[which(f>max.useful.freq | f<max.useful.freq*-1)[1]]     ## ## this significantly increases the accuracy! particularly if applied to the simulated s11. It is OK to remove these frequencies because the actual TDR bandwidth is smaller than max.useful.freq, above it is noise


      zmin.re <- sum( ( Re(s11[1:n.end]) - Re(s11.sim[1:n.end]) )^2 )  +  sum( ( Im(s11[1:n.end]) - Im(s11.sim[1:n.end]) )^2 )  	
       
      ### ### minimise also the distance between Ka measured during TTA and Ka calculated with Cole-Cole at the effective frequency:
      E.probe.feff <- Einf + ( (Es-Einf) / (1 + (j*(f.eff/frel1)^(1-BETA))) ) + ( (Es-Einf) / (1 + (j*(f.eff/frel2))) ) - j*(cond_DC/(2*pi*f.eff*E0))                    
      Er.feff <- Re(E.probe.feff)
      Ei.feff <- Im(E.probe.feff) 
      min.feff <- (Ka - ((Er.feff * mu) / 2) * (sqrt(1 + (Ei.feff/Er.feff)^2) + 1))^2
      
      zmin <- zmin.re + min.feff/4        # final value to minimise (reduce a bit the weight of min.feff)
      
      }

       
optim.freq.domain <- GenSA(control=list(smooth=F, max.time=3), fn = min.fun.freq, par=c(
                 # optim_sa(fun=min.fun.freq, maximization = FALSE, trace = F, start=c(         
                 #GenSA(control=list(smooth=F, max.time=5), fn = min.fun.freq, par=c(
                 #optim(method = "L-BFGS-B", fn = min.fun.freq, control=list(maxit=1000), par=c( 
                 Es.ini, 
                 frel1.ini, 
                 frel2.ini, 
                 cond_DC.ini)
                 , 
                 lower=c(
                 Es.low, 
                 frel1.low, 
                 frel2.low, 
                 cond_DC.low)
                 , 
                 upper=c(
                 Es.up, 
                 frel1.up, 
                 frel2.up, 
                 cond_DC.up)
                 )   

  BETA.freq <- 0 #round(optim.freq.domain$par[1], digits=4)
  Es.freq <- round(optim.freq.domain$par[1], digits=4)
  Einf.freq <- Einf.ini   
  frel1.freq <- round(optim.freq.domain$par[2], digits=0)
  frel2.freq <- round(optim.freq.domain$par[3], digits=0)     
  cond_DC.freq <- round(optim.freq.domain$par[4], digits=4)   
  z.freq <- Re(z.ini) 
               
  # recalculate s11.sim based on optimised parameters in the freq domain

#        # calculate input.wave based on calibrated alpha and t0                     
#        xx <- data.frame(1:n)
#        input.func <- function(xx) { (aa + bb * erf(alpha*(xx-t0.freq))) }
#        input.wave <- apply(X=xx, MARGIN=1, FUN=input.func)
#        Nic.input.wave <- matrix(ncol=1, nrow=n)
#        ramp.input.wave <- input.wave[n]*(x/n)         
#        Nic.input.wave <- input.wave - ramp.input.wave   # subtract ramp
#        Nic.input.wave[1] <- 0
#        # fast fourier transform
#        fft.input.wave <- fft(Nic.input.wave)
#        V0 <- fft.input.wave
#        # recalculate measured s11 using optimised V0 for a correct comparison of sum of squares between s11 and s11.sim
#        s11 <- R/V0   
    
  #E.probe.freq <- Einf + ( (Es.freq-Einf) / (1 + (j*(f/frel.freq)^(1-BETA.freq))) ) - j*(cond_DC.freq/(2*pi*f*E0))
  E.probe.freq <- Einf.freq + ( (Es.freq-Einf.freq) / (1 + (j*(f/frel1.freq)^(1-BETA.freq))) ) + ( (Es.freq-Einf.freq) / (1 + (j*(f/frel2.freq))) ) - j*(cond_DC.freq/(2*pi*f*E0))
  z <- z.freq
   
  #complex(real=z.freq, imaginary=zi.freq)        
  ro <- (1-z*sqrt(E.probe.freq)) / (1+z*sqrt(E.probe.freq))
  gammaL <- (j*2*pi*f*L*sqrt(E.probe.freq))/c0                              

  ## simulated s11 (Clarkson et al., 1977)
            s11.sim.freq <- (ro+exp(-2*gammaL)) / (1+ro*exp(-2*gammaL))               # open-circuit (probe)
            #s11.sim.freq <- (ro-exp(-2*gammaL)) / (1-ro*exp(-2*gammaL))              # short-circuit (cell)
            #s11.sim.freq <- ro * (1-exp(-2*gammaL)) / (1-(ro^2)*exp(-2*gammaL) )     # matched-circuit           
            ##s11.sim.freq[which(f>max.useful.freq | f<max.useful.freq*-1)] <- s11.sim.freq[which(f>max.useful.freq | f<max.useful.freq*-1)[1]]     ## ## this significantly increases the accuracy! particularly if applied to the simulated s11. It is OK to remove these frequencies because the actual TDR bandwidth is smaller than max.useful.freq, above it is noise
            
#            mag.s11.sim.freq <- Mod(s11.sim.freq)      # magnitude simulated s11 
#            phase.s11.sim.freq <- Arg(s11.sim.freq)    # phase simulated s11 
#            rl.s11.sim.freq <- 20*log10(mag.s11.sim.freq)     # return loss simulated s11
      
      ## convert to time-domain 
      R.freq <- s11.sim.freq*V0
      ## R.freq[1] <- 0+0i
      i.s11.sim.freq <- fft(R.freq, inverse=T)/length(R.freq)
      if(Nicolson == F) {
            i.wave.sim.freq <- cumsum(Re(i.s11.sim.freq))
            } else {
            i.wave.sim.freq <- Re(i.s11.sim.freq)
            # normalise waveform to have first point = 0         
            i.wave.sim.freq <- i.wave.sim.freq-i.wave.sim.freq[1]
            }    
       
      # calculate R squared (see Todeschini's book, pag 124, 125)
      tss.freq <- sum( ( i.wave - mean(i.wave) )^2 )     # tot sum of squares
      rss.freq <- sum( ( i.wave.sim.freq - i.wave )^2 )     # residual sum of squares
      r2.freq <- round(1 - (rss.freq/tss.freq), digits=4)





               #stop()















           

















                      

## ## ## ## OPTIMISATION IN THE TIME DOMAIN



min.fun.time <- function(r){
## simplex, optimisation in the time domain 


#      BETA <- r[1]
      Es <- r[1]
      Einf <- Einf.ini #r[3]
      frel1 <- r[2]
      frel2 <- r[3]            
      cond_DC <- r[4] 
#      alpha <- r[5]
#      t0 <- r[7]
#      #z <- r[7]
#
#        # calculate input.wave based on calibrated alpha and t0                     
#        xx <- data.frame(1:n)
#        input.func <- function(xx) { (aa + bb * erf(alpha*(xx-t0))) }
#        input.wave <- apply(X=xx, MARGIN=1, FUN=input.func)
#        Nic.input.wave <- matrix(ncol=1, nrow=n)
#        ramp.input.wave <- input.wave[n]*(x/n)         
#        Nic.input.wave <- input.wave - ramp.input.wave   # subtract ramp
#        Nic.input.wave[1] <- 0
#        # fast fourier transform
#        fft.input.wave <- fft(Nic.input.wave)
#        V0 <- fft.input.wave

      #E.probe <- Einf + ( (Es-Einf) / (1 + (j*(f/frel)^(1-BETA))) ) - j*(cond_DC/(2*pi*f*E0))
      E.probe <<- Einf + ( (Es-Einf) / (1 + (j*(f/frel1)^(1-BETA))) ) + ( (Es-Einf) / (1 + (j*(f/frel2))) ) - j*(cond_DC/(2*pi*f*E0))         
      E.probe[1] <<- complex(real=0, imaginary=0)      
#      z <- complex(real=z, imaginary=zi)                
      #z <- z.ini
      ro <- (1-z*sqrt(E.probe)) / (1+z*sqrt(E.probe))
      gammaL <- (j*2*pi*f*L*sqrt(E.probe))/c0              

      ## simulated s11 (Clarkson et al., 1977)
      s11.sim <- (ro+exp(-2*gammaL)) / (1+ro*exp(-2*gammaL))    # open-circuit (probe)    
      ##s11.sim[which(f>max.useful.freq | f<max.useful.freq*-1)] <- s11.sim[which(f>max.useful.freq | f<max.useful.freq*-1)[1]]     ## ## this significantly increases the accuracy! particularly if applied to the simulated s11. It is OK to remove these frequencies because the actual TDR bandwidth is smaller than max.useful.freq, above it is noise
      
      ## convert to time-domain
      R.time <- s11.sim*V0
      ## R.time[1] <- 0+0i
      i.s11.sim <- fft(R.time, inverse=T)/length(R.time)
      if(Nicolson == F) {
            i.wave.sim <- cumsum(Re(i.s11.sim[1:n]))
            } else {
            i.wave.sim <- Re(i.s11.sim[1:n])
            # normalise waveform to have first point = 0        
            i.wave.sim <- i.wave.sim-i.wave.sim[1]  
            }        
 
      zmin <- sum( ( i.wave[round(startpoint, digits=0):n] - i.wave.sim[round(startpoint, digits=0):n] )^2 )   # do not use the part of the waveform before startpoint for optimisation          
            
      ### ### minimise also the distance between Ka measured during TTA and Ka calculated with Cole-Cole at the effective frequency:
      E.probe.feff <- Einf + ( (Es-Einf) / (1 + (j*(f.eff/frel1)^(1-BETA))) ) + ( (Es-Einf) / (1 + (j*(f.eff/frel2))) ) - j*(cond_DC/(2*pi*f.eff*E0))                    
      Er.feff <- Re(E.probe.feff)
      Ei.feff <- Im(E.probe.feff) 
      min.feff <- (Ka - ((Er.feff * mu) / 2) * (sqrt(1 + (Ei.feff/Er.feff)^2) + 1))^2
      
      zmin <- zmin + min.feff/4        # final value to minimise (reduce a bit the weight of min.feff)
      
      }
                           
optim.time.domain <- GenSA(control=list(smooth=F, max.time=3), fn = min.fun.time, par=c( 
                 # optim_sa(fun=min.fun.freq, maximization = FALSE, trace = F, start=c(         
                 #GenSA(control=list(smooth=F, max.time=5), fn = min.fun.freq, par=c(
                 #optim(method = "L-BFGS-B", fn = min.fun.freq, control=list(maxit=1000), par=c(
                 Es.ini, 
                 frel1.ini, 
                 frel2.ini, 
                 cond_DC.ini)
                 , 
                 lower=c(
                 Es.low, 
                 frel1.low, 
                 frel2.low, 
                 cond_DC.low)
                 , 
                 upper=c(
                 Es.up, 
                 frel1.up, 
                 frel2.up, 
                 cond_DC.up)
                 )   

  BETA.time <- 0 #round(optim.time.domain$par[1], digits=4)
  Es.time <- round(optim.time.domain$par[1], digits=4)
#  deltaK1.time <- round(optim.time.domain$par[2], digits=2)
#  deltaK2.time <- round(optim.time.domain$par[3], digits=2)
  Einf.time <- Einf.ini #round(optim.time.domain$par[3], digits=2)  
  frel1.time <- round(optim.time.domain$par[2], digits=0)
  frel2.time <- round(optim.time.domain$par[3], digits=0)     
  cond_DC.time <- round(optim.time.domain$par[4], digits=4)   
  z.time <- Re(z.ini)
                 
  # recalculate s11.sim based on optimised parameters in the time domain

#        # calculate input.wave based on calibrated alpha and t0                     
#        xx <- data.frame(1:n)
#        input.func <- function(xx) { (aa + bb * erf(alpha*(xx-t0.time))) }
#        input.wave <- apply(X=xx, MARGIN=1, FUN=input.func)
#        Nic.input.wave <- matrix(ncol=1, nrow=n)
#        ramp.input.wave <- input.wave[n]*(x/n)         
#        Nic.input.wave <- input.wave - ramp.input.wave   # subtract ramp
#        Nic.input.wave[1] <- 0
#        # fast fourier transform
#        fft.input.wave <- fft(Nic.input.wave)
#        V0 <- fft.input.wave
#        # recalculate measured s11 using optimised V0 for a correct comparison of sum of squares between s11 and s11.sim
#        s11 <- R/V0  
     
  #E.probe.time <- Einf + ( (Es.time-Einf) / (1 + (j*(f/frel.time)^(1-BETA.time))) ) - j*(cond_DC.time/(2*pi*f*E0))
  E.probe.time <- Einf.time + ( (Es.time-Einf.time) / (1 + (j*(f/frel1.time)^(1-BETA.time))) ) + ( (Es.time-Einf.time) / (1 + (j*(f/frel2.time))) ) - j*(cond_DC.time/(2*pi*f*E0))
  z <- z.time
  #complex(real=z.time, imaginary=zi.time)        
  ro <- (1-z*sqrt(E.probe.time)) / (1+z*sqrt(E.probe.time))
  gammaL <- (j*2*pi*f*L*sqrt(E.probe.time))/c0                              

  ## simulated s11 (Clarkson et al., 1977)
            s11.sim.time <- (ro+exp(-2*gammaL)) / (1+ro*exp(-2*gammaL))               # open-circuit (probe)
            #s11.sim.time <- (ro-exp(-2*gammaL)) / (1-ro*exp(-2*gammaL))              # short-circuit (cell)
            #s11.sim.time <- ro * (1-exp(-2*gammaL)) / (1-(ro^2)*exp(-2*gammaL) )     # matched-circuit           
            ##s11.sim.time[which(f>max.useful.time | f<max.useful.time*-1)] <- s11.sim.time[which(f>max.useful.time | f<max.useful.time*-1)[1]]     ## ## this significantly increases the accuracy! particularly if applied to the simulated s11. It is OK to remove these timeuencies because the actual TDR bandwidth is smaller than max.useful.time, above it is noise
            
#            mag.s11.sim.time <- Mod(s11.sim.time)      # magnitude simulated s11 
#            phase.s11.sim.time <- Arg(s11.sim.time)    # phase simulated s11 
#            rl.s11.sim.time <- 20*log10(mag.s11.sim.time)     # return loss simulated s11
      
      ## convert to time-domain 
      R.time <- s11.sim.time*V0
      ## R.time[1] <- 0+0i
      i.s11.sim.time <- fft(R.time, inverse=T)/length(R.time)
      if(Nicolson == F) {
            i.wave.sim.time <- cumsum(Re(i.s11.sim.time))
            } else {
            i.wave.sim.time <- Re(i.s11.sim.time)
            # normalise waveform to have first point = 0         
            i.wave.sim.time <- i.wave.sim.time-i.wave.sim.time[1]
            }    
       
      # calculate R squared (see Todeschini's book, pag 124, 125)
      tss.time <- sum( ( i.wave - mean(i.wave) )^2 )     # tot sum of squares
      rss.time <- sum( ( i.wave.sim.time - i.wave )^2 )     # residual sum of squares
      r2.time <- round(1 - (rss.time/tss.time), digits=4)


      ## calculate magnitude of dispersion
          Er_100_freq <- Re(E.probe.freq[which(f > 100e6)[1]])
          Er_1000_freq <- Re(E.probe.freq[which(f > 1e9)[1]])
          Ei_100_freq <- Im(E.probe.freq[which(f > 100e6)[1]])
          Ei_1000_freq <- Im(E.probe.freq[which(f > 1e9)[1]])
          Ka_100_freq <- ((Er_100_freq * mu) / 2) * (sqrt(1 + (Ei_100_freq/Er_100_freq)^2) + 1)
          Ka_100_freq <- round(Ka_100_freq, digits=2)
          Ka_1000_freq <- ((Er_1000_freq * mu) / 2) * (sqrt(1 + (Ei_1000_freq/Er_1000_freq)^2) + 1)
          Ka_1000_freq <- round(Ka_1000_freq, digits=2)
          
          Er_100_time <- Re(E.probe.time[which(f > 100e6)[1]])
          Er_1000_time <- Re(E.probe.time[which(f > 1e9)[1]])
          Ei_100_time <- Im(E.probe.time[which(f > 100e6)[1]])
          Ei_1000_time <- Im(E.probe.time[which(f > 1e9)[1]])
          Ka_100_time <- ((Er_100_time * mu) / 2) * (sqrt(1 + (Ei_100_time/Er_100_time)^2) + 1)
          Ka_100_time <- round(Ka_100_time, digits=2)
          Ka_1000_time <- ((Er_1000_time * mu) / 2) * (sqrt(1 + (Ei_1000_time/Er_1000_time)^2) + 1)
          Ka_1000_time <- round(Ka_1000_time, digits=2)

       # note! if BETA = 0 Er Mdisp = 0, but Ka Mdisp should still be  useful
       Mdisp_Er_freq <- round(abs(Er_100_freq - Er_1000_freq), digits=2)        # magnitude of dispersion of real permittivity 100MHz-1GHz
       Mdisp_Er_time <- round(abs(Er_100_time - Er_1000_time), digits=2)        # magnitude of dispersion of real permittivity 100MHz-1GHz
       Mdisp_freq <- round(abs(Ka_100_freq - Ka_1000_freq), digits=2)           # magnitude of dispersion of apparent permittivity 100MHz-1GHz
       Mdisp_time <- round(abs(Ka_100_time - Ka_1000_time), digits=2)           # magnitude of dispersion of apparent permittivity 100MHz-1GHz

       print(c("Er Mdisp_freq (100MHz-1GHz):", Mdisp_Er_freq), quote=F)
       print(c("Er Mdisp_time (100MHz-1GHz):", Mdisp_Er_time), quote=F)         
       print(c("Ka Mdisp_freq (100MHz-1GHz):", Mdisp_freq), quote=F)
       print(c("Ka Mdisp_time (100MHz-1GHz):", Mdisp_time), quote=F)












## ## ## ## plots

windows(w=12, h=7)
plot(f[1:n.end], Re(s11[1:n.end]), type='o', pch=1, xlim=c(0, f.end), ylim=c(-1,1), xlab=NA, ylab=NA )
par(new=T)
plot(f[1:n.end], Re(s11.sim.freq[1:n.end]), type='o', pch=4, lty=2, xlim=c(0, f.end), ylim=c(-1,1), xlab="frequency (Hz)", ylab="Re(s11)" )
      legend("topright",
            c("measured","simulated"),
            lty=c(1,2),
            bg="white"
            )
title("freq domain optimisation")
savePlot(filename = paste(file.name, "_Re(s11).png", sep=''), type = c("png"))
dev.off()

windows(w=12, h=7)

## measured and simulated waveforms
      plot(i.wave, type='l', ylim=c(-1.2,1.2), ylab=NA)            
      par(new=T)
      plot(i.wave.sim.freq[1:n0], type='l', lty=2, ylim=c(-1.2,1.2), ylab=NA)
      legend("bottomright",
         c("measured", "simulated"),
         lty=c(1,2),
         bg="white"
         )
      title("freq domain optimisation")
      legend("bottomleft",
             legend=c(
             paste("BETA.freq = ", BETA.freq, sep=''),
             paste("Es.freq = ", Es.freq, sep=''),
             paste("frel1.freq = ", frel1.freq, sep=''),
             paste("frel2.freq = ", frel2.freq, sep=''),
             paste("Einf.freq = ", Einf.freq, sep=''),
             paste("cond_DC.freq = ", cond_DC.freq, sep=''),
             paste("z.freq = ", z.freq, sep='')
), bty='n')
      text( x=n0/2, y=1.1, labels=paste("RSS = ", round(rss.freq, digits=4), sep="") )      
      par(new=T); plot(input.wave[1:n0], type='l', col='black', ylim=c(-1.2,1.2), ylab="reflection coefficient")
      savePlot(filename = paste(file.name, "_wave.sim.freq.png", sep=''), type = c("png"))
      dev.off()
      
      windows(w=12, h=7)

      plot(i.wave[1:n0], type='l', ylim=c(-1.2,1.2), ylab=NA, xlim=c(0,n0))              
      par(new=T)
      plot(i.wave.sim.time[1:n0], type='l', lty=2, ylim=c(-1.2,1.2), ylab=NA, xlim=c(0,n0))   
      legend("bottomright",
         c("measured", "simulated"),
         lty=c(1,2),
         bg="white"
         )
      title("time domain optimisation")
      legend("bottomleft",
             legend=c(
             paste("BETA.time = ", BETA.time, sep=''),
             paste("Es.time = ", Es.time, sep=''),
             paste("frel1.time = ", frel1.time, sep=''),
             paste("frel2.time = ", frel2.time, sep=''),
             paste("Einf.time = ", Einf.time, sep=''),
             paste("cond_DC.time = ", cond_DC.time, sep=''),
             paste("z.time = ", z.time, sep='')
), bty='n')
      text( x=n0/2, y=1.1, labels=paste("RSS = ", round(rss.time, digits=4), sep="") )          
      par(new=T); plot(input.wave[1:n0], type='l', col='black', ylim=c(-1.2,1.2), ylab="reflection coefficient", xlim=c(0,n0))   
      savePlot(filename = paste(file.name, "_wave.sim.time.png", sep=''), type = c("png"))
      dev.off()

                       # save objects with shorter names for quicker analysis
                       E.t <- E.probe.time
                       Er.t <- Re(E.t)
                       Ei.t <- Im(E.t)
                       E.f <- E.probe.freq
                       Er.f <- Re(E.f)
                       Ei.f <- Im(E.f)

                       E.t.eff <- E.t[which(f > f.eff)[1]]
                       E.f.eff <- E.f[which(f > f.eff)[1]]
                       Er_eff_time <- Re(E.t.eff)
                       Ei_eff_time <- Im(E.t.eff)
                       Er_eff_freq <- Re(E.f.eff)
                       Ei_eff_freq <- Im(E.f.eff)
                      
                       Ea.f <- ((Er.f * mu) / 2) * (sqrt(1 + (Ei.f/Er.f)^2) + 1)
                       Ka_eff_freq <- round(Ea.f[which(f > f.eff)[1]], digits=2)    # print apparent permittivity calculated from Er and Ei at the effective frequency
                       Ea.t <- ((Er.t * mu) / 2) * (sqrt(1 + (Ei.t/Er.t)^2) + 1)
                       Ka_eff_time <- round(Ea.t[which(f > f.eff)[1]], digits=2)    # print apparent permittivity calculated from Er and Ei at the effective frequency
                             
windows(h=7, w=12)
plot(f[1:(n/2)], Ea.t[1:(n/2)], xlim=c(0,1e9), ylim=c(0,50), t='o'); 
par(new=T); 
plot(f[1:(n/2)], Ea.f[1:(n/2)], xlim=c(0,1e9), ylim=c(0,50), t='o', col='red')
legend("topright", legend=c("Ka (time optim)", "Ka (freq optim)", paste("Ka_eff_time = ", Ka_eff_time), paste("Ka_eff_freq = ", Ka_eff_freq), paste("Ka = ", Ka)), text.col=c("black", "red", "black", "red", "blue"))                             
savePlot(filename = paste(file.name, "_Ka_complex.png", sep=''), type = c("png"))

    #  windows(w=12, h=7)
     dev.off()

print(c("f.eff:", round(f.eff*10^-6, digits=0), "MHz"), quote=F)             # effective frequency (MHz)
print(c("f.up.circ:", round(f.up.circ*10^-6, digits=0), "MHz"), quote=F)     # circumferential resonance freq. (MHz)
print(c("f.up.long:", round(f.up.long*10^-6, digits=0), "MHz"), quote=F)     # longitudinal resonance freq. (MHz), probably more important!


##}

} # *** #  end of if(inversion == TRUE)
 

                       # save the output in an object (file.name, Ka, VWC (Topp))
                       output.Ka[ii,] <- c(file.name, Ka, vwc.topp)                             
                       output.tdr[ii,] <- c(file.name, Ka, BEC, vwc.topp, 
                                               pto, Lapp, V1, V1m, Vf, 
                                               A.m1.der, A.m2.der, A.end.plus.der, A.pulse.der,      
                                               min.ref.x, min.ref.y, refpoint.x, refpoint.y, inflection1.x, inflection1.y, 
                                               max.start.x, max.start.y, startpoint.x, startpoint.y, 
                                               min.end.x, min.end.y, endpoint.x, endpoint.y, inflection2.x, inflection2.y)

                       if (inversion == TRUE) {
                       # save the fft output in an object
                       output.fft[ii,] <- c(
                                           file.name, 
                                           BETA.time, Es.time, frel1.time, frel2.time, Einf.time, cond_DC.time, z.time, Mdisp_time, r2.time, 
                                           BETA.freq, Es.freq, frel1.freq, frel2.freq, Einf.freq, cond_DC.freq, z.freq, Mdisp_freq, r2.freq,
                                           f.eff, Ka_eff_time, Ka_eff_freq, Er_eff_time, Ei_eff_time, Er_eff_freq, Ei_eff_freq, 
                                           Ka_100_time, Ka_100_freq, Er_100_time, Ei_100_time, Er_100_freq, Ei_100_freq,
                                           Ka_1000_time, Ka_1000_freq, Er_1000_time, Ei_1000_time, Er_1000_freq, Ei_1000_freq
                                           )                
                       }              

                 
      ## add temperature
      file.name.col <- output.tdr[ii,1]
              
              ## split file.name by '-' or '.' (with the exception of temperature that is contained by two '_')
              #Splits <- unlist(strsplit(file.name.col[1], "[.]"))  # it assumes that each waveform is named in a standard way
              Splits <- unlist(strsplit(file.name.col, "-"))  # it assumes that each waveform is named in a standard way
              firstNA <- which(is.na(Splits))
              
              if(length(firstNA) == 0) { 
                                 temperatura[ii] <- as.numeric(substring(Splits[length(Splits)], 1, (nchar(Splits[length(Splits)])-1))) } else{
                                 temperatura[ii] <- as.numeric(substring(Splits[firstNA-1], 1, (nchar(Splits[firstNA-1])-1)))
                                 }
              
#              # add extra columns with additional information from file name
#              extracol$probe.name[ii] <- Splits[1]
#              extracol$mux[ii] <- Splits[which(Splits == "mux0" | Splits == "mux1" | Splits == "mux2")]
#              extracol$repetition[ii] <- Splits[pmatch("rep", Splits, nomatch = NA_integer_, duplicates.ok = FALSE)]       # partial match!
#              extracol$sample_rep[ii] <- unlist(strsplit(extracol$repetition[ii], "[.]"))
             
 

graphics.off()


# reset initial settings for next waveform analysis (in case multiple waveforms are analysed at the same time)
interactive.sol <- interactive.sol.ini
interactive.tan1 <- interactive.tan1.ini
interactive.tan2 <- interactive.tan2.ini
interactive.tan3 <- interactive.tan3.ini
interactive.tan4 <- interactive.tan4.ini
probe.name <- probe.name.ini
shorted <- shorted.ini
tdr.name <- tdr.name.ini
use.end.tangent <- use.end.tangent.ini
use.ref.tangent <- use.ref.tangent.ini

    # restore original values of n and pto
    n <- n0
    pto <- win.Lapp/n     # apparent length (m) of a single point
                                                    

#}     # end of for(ii in 1:n.waves)

 

                       # save the output in an object (file.name, Ka, VWC (Topp))
                       output.Ka[ii,] <- c(file.name, Ka, vwc.topp)                             
                       output.tdr[ii,] <- c(file.name, Ka, BEC, vwc.topp, 
                                               pto, Lapp, V1, V1m, Vf, 
                                               A.m1.der, A.m2.der, A.end.plus.der, A.pulse.der,      
                                               min.ref.x, min.ref.y, refpoint.x, refpoint.y, inflection1.x, inflection1.y, 
                                               max.start.x, max.start.y, startpoint.x, startpoint.y, 
                                               min.end.x, min.end.y, endpoint.x, endpoint.y, inflection2.x, inflection2.y)

                       if (inversion == TRUE) {
                       # save the fft output in an object
                       output.fft[ii,] <- c(
                                           file.name, 
#                                          Ka, BEC, vwc.topp, 
#                                          pto, Lapp, V1, V1m, Vf, 
#                                          A.m1.der, A.m2.der, A.end.plus.der, A.pulse.der,      
#                                          min.ref.x, min.ref.y, refpoint.x, refpoint.y, inflection1.x, inflection1.y, 
#                                          max.start.x, max.start.y, startpoint.x, startpoint.y, 
#                                          min.end.x, min.end.y, endpoint.x, endpoint.y, inflection2.x, inflection2.y,
                                           BETA.time, Es.time, frel1.time, frel2.time, Einf.time, cond_DC.time, z.time, Mdisp_time, r2.time, 
                                           BETA.freq, Es.freq, frel1.freq, frel2.freq, Einf.freq, cond_DC.freq, z.freq, Mdisp_freq, r2.freq,
                                           f.eff, Ka_eff_time, Ka_eff_freq, Er_eff_time, Ei_eff_time, Er_eff_freq, Ei_eff_freq,
                                           Ka_100_time, Ka_100_freq, Er_100_time, Ei_100_time, Er_100_freq, Ei_100_freq,
                                           Ka_1000_time, Ka_1000_freq, Er_1000_time, Ei_1000_time, Er_1000_freq, Ei_1000_freq
                                           )                
                       }                     

graphics.off()


# reset initial settings for next waveform analysis (in case multiple waveforms are analysed at the same time)
interactive.sol <- interactive.sol.ini
interactive.tan1 <- interactive.tan1.ini
interactive.tan2 <- interactive.tan2.ini
interactive.tan3 <- interactive.tan3.ini
interactive.tan4 <- interactive.tan4.ini
probe.name <- probe.name.ini
shorted <- shorted.ini
tdr.name <- tdr.name.ini
use.end.tangent <- use.end.tangent.ini
use.ref.tangent <- use.ref.tangent.ini

    # restore original values of n and pto
    n <- n0
    pto <- win.Lapp/n     # apparent length (m) of a single point
                                                    

}     # end of for(ii in 1:n.waves)

            
            
                  output.tdr <- output.tdr %>%
                                mutate_at(colnames(output.tdr)[-1], as.numeric) #%>%
#                                add_column(temp = temperatura)
#                  extracol <- extracol[colSums(!is.na(extracol)) > 0]    # remove columns with only NAs          
#                  output.tdr <- cbind(output.tdr, extracol, split.name) 
                    
                  output.Ka <- output.Ka %>%
                               mutate_at(colnames(output.Ka)[-1], as.numeric)
                               output.Ka[,1] <- as.character(output.Ka[,1])
                    
                  output.BEC <- output.BEC %>%
                               mutate_at(colnames(output.BEC)[-1], as.numeric)
                               output.BEC[,1] <- as.character(output.BEC[,1])
           

                      # export the output to a .txt file
#                       write.table(output.Ka, paste(probe.name, ".Ka.txt", sep=''), sep='\t', row.names=F, quote=F, append=T)
#                       write.table(output.tdr, paste(probe.name, ".output.tdr.txt", sep=""), sep='\t', row.names=F, quote=F, append=T)
#                       write.table(output.BEC, paste(probe.name, ".BEC.txt", sep=''), sep='\t', row.names=F, col.names=c("file.name.bec", "BEC.mS.m", "RC"), quote=F, append=T)
                       write_delim(output.Ka, paste(probe.name, ".Ka.txt", sep=''), delim='\t', append = T, col_names=T)
                       write_delim(output.tdr, paste(probe.name, ".output.tdr.txt", sep=''), delim='\t', append = T, col_names=T)
                       write_delim(output.BEC, paste(probe.name, ".BEC.txt", sep=''), delim='\t', append = T, col_names=T)

# * #                       
if(inversion == TRUE) {
                      output.fft <- data.frame(output.fft)
                      output.fft <- output.fft %>%
                                    mutate_at(colnames(output.fft)[-1], as.numeric) 
                      output.fft <- cbind(output.tdr, output.fft[,-1])
                        
                       #write.table(output.fft, paste(probe.name, ".fft.txt", sep=''), sep='\t', row.names=F, quote=F, append=T)
                       write_delim(output.fft, paste(probe.name, ".fft.txt", sep=''), delim='\t', append = T, col_names=T)                      
                       
options(warn=0)     # turn warnings back on
                       
message(paste(c("at the effective frequency of ", round(f.eff/10^6, digits=0), " MHz, E.f = ", round(E.f.eff, digits=2)), sep=''))
message(paste(c("at the effective frequency of ", round(f.eff/10^6, digits=0), " MHz, E.t = ", round(E.t.eff, digits=2)), sep=''))
message(paste(c("at the effective frequency of ", round(f.eff/10^6, digits=0), " MHz, Ea.f = ", round(Ka_eff_freq, digits=2)), sep=''))
message(paste(c("at the effective frequency of ", round(f.eff/10^6, digits=0), " MHz, Ea.t = ", round(Ka_eff_time, digits=2)), sep=''))

} # * #                   

##}
