##tdr.cal.perm <- function(...) {
# Giulio Curioni
# 15 July 2019


# perform AIR-WATER or ACETONE-WATER calibration for permittivity

#check ### before running the script (note: search ### ### for important points for TDR inversion)
        
         source("C:\\R.wd\\loadall.R")        # load packages and functions
         loadall()
         graphics.off()
         options(warn=-1)       # suppress warnings


message("\n\nStep1 - Take TDR waveforms in distilled water and air (possibly shorted), or in acetone, using an appropriate window length so that the end reflection is clearly visible in water.\n        Avoid using long waveforms. ALSO, record the temperature of water and acetone if used.")
message("Step2 - Analyse waveforms using tdr.perm.tta.R or tdr.perm.tta.acclima.R. Probes must be analysed separately and the script must be run separately for each medium,\n        particularly if shorted measurements in air have been taken.")
message("Step3 - Manually rename results of analysis (i.e. probe.name.Ka.txt) to be able to distinguish between water and second reference medium.")
answer1 <- readline("\nAre steps 1 to 3 completed? press ENTER to continue, press any other key to cancel.\n")                           
if(answer1 != "") {stop("Complete steps 1 & 2 before proceeding.")}
waterT <- readline("\nInput water temperature during calibration.\n"); waterT <- as.numeric(waterT)
medium2T <- readline("\nIf second reference medium was acetone input its temperature during calibration. Otherwise press ENTER. Do NOT input temperature if medium2 is air.\n")
medium2T  <- as.numeric(medium2T ) 

# import results of analysis
path.name <- data.frame(choose.files(caption = "Select file with analysed WATER measurements. It should be in R.wd (i.e. probe.name.Ka.txt)"))
path.name <- as.character(path.name[1,1])
water <- read_delim(path.name, delim="\t", col_names = TRUE)

path.name <- data.frame(choose.files(caption = "Select file with analysed AIR (possibly shorted) or ACETONE measurements. It should be in R.wd (i.e. probe.name.Ka.txt)"))
path.name <- as.character(path.name[1,1])
medium2 <- read_delim(path.name, delim="\t", col_names = TRUE)    # this in general is air

# take averages
water$tt <- water$endpoint - water$refpoint
water_tt <- mean(water$tt, na.rm=T)     # in case of repeat measurements
medium2$tt <- medium2$endpoint - medium2$refpoint
medium2_tt <- mean(medium2$tt, na.rm=T)     # in case of repeat measurements

# convert data points to time units
pto_ps <- water$pto_ps[1]
 
 
if(length(pto_ps) != 0) {


# Acclima, use actual time (s)  ------------------------------------------------

   
if(water$pto_ps[1] != medium2$pto_ps[1]) {stop("Check that the step between two points is consistent for all measurements.\nDo not change win.Lapp or traveltime! if this happens repeat the measurements.")}
water_tt <- water_tt * pto_ps * 10^-12  # time in s
medium2_tt <- medium2_tt * pto_ps * 10^-12  # time in s
# reference permittivity values
water_perm <- 78.54 * ( 1-4.579e-3*(waterT-25) + 1.19e-5*((waterT-25)^2) - 2.8e-8*((waterT-25)^3) )  # Weast (1986) static permittivity of free water at temperature = waterT
# if reference temperature is not available it is assumed that the reference second medium is air, otherwise it is acetone
if(is.na(medium2T)) {medium2_perm <- 1.00059} else {medium2_perm <- 23.2047 - 0.1017*medium2T}
# calculated calibrated L0 (offset from refpoint to startpoint in time units) and calibrated Lcal (length of probe in metres)
c0 <- 2.9979e8              # speed of light (m/s)
Lcal_tt <- (medium2_tt - water_tt)/(sqrt(medium2_perm)-sqrt(water_perm))      
Lcal <- (Lcal_tt*c0)/2      # divide by 2 to account for two way travel
L0 <- ((water_tt*c0)/2) - (Lcal*sqrt(water_perm))         # careful with units, I can't sum times with lengths. Convert water travel times to apparent length.


} else {


# Campbell Scientific (use apparent length, m) ---------------------------------


winLapp <- readline("\nInput the window length (winLapp) used during calibration and press ENTER. The script works only if 2048 data points were taken, if a different number was used change manually the script.\n")
winLapp  <- as.numeric(winLapp) 
nn <- 2048      ### ### change manually if a different number of data points was used
pto_ps <- winLapp/nn
water_tt <- water_tt * pto_ps # time in m
medium2_tt <- medium2_tt * pto_ps # time in m
# reference permittivity values
water_perm <- 78.54 * ( 1-4.579e-3*(waterT-25) + 1.19e-5*((waterT-25)^2) - 2.8e-8*((waterT-25)^3) )  # Weast (1986) static permittivity of free water at temperature = waterT
# if reference temperature is not available it is assumed that the reference second medium is air, otherwise it is acetone
if(is.na(medium2T)) {medium2_perm <- 1.00059} else {medium2_perm <- 23.2047 - 0.1017*medium2T}
# calculated calibrated L0 (offset from refpoint to startpoint in time units) and calibrated Lcal (length of probe in metres)
c0 <- 2.9979e8              # speed of light (m/s)
Lcal_tt <- (medium2_tt - water_tt)/(sqrt(medium2_perm)-sqrt(water_perm))      
Lcal <- Lcal_tt
L0 <- water_tt - (Lcal*sqrt(water_perm))         # careful with units, I can't sum times with lengths. 

Lcal <- Lcal_tt   
L0 <- water_tt - (Lcal*sqrt(water_perm))         # careful with units, I can't sum times with lengths. 


}


Lcal <- round(Lcal, digits=5)
L0 <- round(L0, digits=5)

message(paste("Lcal and L0 are, respectively (m):", as.character(Lcal), as.character(L0), sep=" "))

message("\nCopy the above values to tdr.probe.par.txt in C:/R.wd/scripts/tdr")
         
