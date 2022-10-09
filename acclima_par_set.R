acclima_par_set <- function(winLapp) {
# Giulio Curioni
# 17 AUgust 2018

# convert units used by Acclima software to the units used by Campbell Scientific (e.g. time to apparent length etc.)
# also make sure that output has 2048 data points

# check ### before running the script

### input: desired winLapp in metres        
         
# function to round to the nearest 5 (max resolution (time step) is 5 picoseconds
mround <- function(x,base){
          base*round(x/base)
          } 

# space = velocity * time
N <- 2048 ### change if necessary (number of data points)
          ### according to Acclima email, the max time is 20000 ps. In fact it can be longer, but given the max time resolution: 20000/5 = 4000 (max data points is 4000)
c0 <- 2.9979e8 # speed of light, m/s
tt <- winLapp/c0
tt <- tt * 10^12 # convert s to ps
time_step <- tt/N # time step between two consecutive points to have N data points

# round to the nearest 5:
time_step <- mround(x=time_step, base=5)
tt <- time_step * N  # recalculate total time so that number of points is N
winLapp_actual <- (tt*10^-12)*c0     # recalculate winLapp so that number of points is N and total time is tt

print(c("time (ps):", tt), quote=F)
print(c("time step (ps):", time_step), quote=F)
print(c("corresponding winLapp (m):", winLapp_actual), quote=F)


### alternatively (this allows far more control on the actual winLapp)

#choose tt
#keep time step = 5
#calculate N e corresponding winLapp (as long as N is more than 2048 it's fine!)

tt0 <- 25000    ### change time (ps)
time_step0 <- 10   # keep it constant
N0 <- tt0/time_step0
winLapp_actual0 <- (tt0*10^-12)*c0
message(" ")
print(c("time (ps):", tt0), quote=F)
print(c("time step (ps):", time_step0), quote=F)
print(c("corresponding winLapp (m):", winLapp_actual0), quote=F)
print(c("n. data points:", N0), quote=F)

if(N0 >=4000) {print("n. data points over the limit of 4000!")}
}
         