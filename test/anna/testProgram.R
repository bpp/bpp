#####
#Testing 3 sequences from 1 population, simple.ctl
times <- read.csv("~/bpp/test/anna/testSimple.txt", header =FALSE)
print("simple.ctl: 3 sequences from 1 population")
firstCoal <- times[seq(from = 1, to = dim(times)[1] -1, by = 2), ]
secCoal <- times[seq(from = 2, to = dim(times)[1], by = 2), ]

# need to condition on time before third lineage is added
# Look at first time and second time
twoLineages <- which(firstCoal[,3] < .05 ) # Indices of runs where coalescent event occurs before last sample is added 
beforeThird <- firstCoal[twoLineages, 3] 
# Mean time of first coalescence event conditional on the event being before the
# third sample is added 
#mean(beforeThird) - .02 # should be truncated exponential 
theta <- 2/.04
b <- .05-.02
print(paste("Theoretical: ", 1/theta - b * (exp(theta * b) - 1)^ -1, ", Experimental: ", mean(beforeThird) - .02))
# Calculating truncated exponential
b <- .05-.02 # Time between second and third samples 
theta <-1/.02

# These should be very close 
#1/theta  -b *(exp(theta*b) -1)^-1 # truncated exponential mean
# Difference between expected and realized time of the first coalescent event conditional on it occurs before the last sampling event
abs (1/theta  -b *(exp(theta*b) -1)^-1 - (mean(beforeThird) - .02)) / (1/theta  -b *(exp(theta*b) -1)^-1 )

# Second coalescent time conditional on the first coalescent event occurring before the last sample was added
time2 <- secCoal[twoLineages ,3] 
#mean(time2) -.05 # should be .020, mean of exponential distribution with rate .020
print(paste("Theoretical: ", .02, ", Experimental: ", mean(time2) -0.05))
abs(mean(time2) -.05 - .02) /(.02)

# If coalescence before third lineage is added, 
# look at first and second coalescence time
threeLineages <- which(firstCoal[,3] > .05 )
time13 <- firstCoal[threeLineages,3]
theo <- 1/(3 * 2/.04)
print(paste("Theoretical: ", theo, ", Experimental: ", (mean(time13) - .05)))
abs(theo  - (mean(time13) - .05))  / theo# Waiting time for coalescent event with three lineages, given there are three lineages at time .05

time23 <- secCoal[threeLineages, 3]
print(paste("Theoretical: ", .02, ", Experimental: ",mean(time23 -time13) ))
abs(mean(time23 -time13) - .02) / .02# Waiting time for coalescent event with two lineages



### Checking formulas with simulation 
# testExp <- rexp(n =1000000, rate =1/.02)
# mean(testExp)
# 
# truncExp <- which(testExp <= .03)
# truncUpper <- which(testExp > .03)
# freqL <- length(truncExp)/length(testExp)
# mean(testExp[truncExp])
# 
# (freqL * mean(testExp[truncExp]))  + ((1 - freqL) * mean(testExp[-truncExp]))
#####

# Testing 2 populations, 2 sequences per population 
# simple2.ctl
times <- read.csv("~/bpp/test/anna/testSimple2.txt", header =FALSE)
print("simple2.ctl: 2 populations, 2 sequences per population ")
firstCoal <- times[seq(from =1, to = dim(times)[1] -2, by = 3), ]
secCoal <- times[seq(from =2, to = dim(times)[1]- 1, by = 3), ]
thirdCoal <- times[seq(from =3, to = dim(times)[1], by = 3), ]

# 0 and 1 coalesce in tip population
rowNum <- intersect(intersect(which(times[ ,1] == 0), which(times[,2] == 1)),
          which(times[,3] < .07))
#mean(times[rowNum, 3])- 0.05 # subtract time of last sample

b <- .07-.05 # speciation time minus time of last sample
theta <-2/.04
print(paste("Theoretical: ", (1/theta  -b *(exp(theta*b) -1)^-1), ", Experimental: ", (mean(times[rowNum, 3])- 0.05)))
abs((1/theta  -b *(exp(theta*b) -1)^-1) - (mean(times[rowNum, 3])- 0.05)) / (1/theta  -b *(exp(theta*b) -1)^-1)

#length(rowNum)/ length(firstCoal[,1])
#PDF
#1- exp(-.02 *theta) # probably coalesce before speciation time .02 = .07-0.05
print(paste("Theoretical: ", (1- exp(-.02 *theta)), ", Experimental: ", length(rowNum)/ length(firstCoal[,1]) ))
abs(length(rowNum)/ length(firstCoal[,1]) - (1- exp(-.02 *theta) )) / (1- exp(-.02 *theta) )

# 2 and 3 coalesce in tip population 
rowNum <- intersect(intersect(which(times[ ,1] == 2), which(times[,2] == 3)),
                    which(times[,3] < .07))
#mean(times[rowNum, 3])- 0.06 # subtract time of last sample

b <- .07-.06 # speciation time minus time of last sample
theta <-2/.03 # k chose 2 is where the 2 is from
print(paste("Theoretical: ", 1/theta  -b *(exp(theta*b) -1)^-1, ", Experimental: ", mean(times[rowNum, 3])- 0.06 ))
abs(1/theta  -b *(exp(theta*b) -1)^-1 - (mean(times[rowNum, 3])- 0.06)) / (1/theta  -b *(exp(theta*b) -1)^-1 )

#length(rowNum)/ length(firstCoal[,1])
#PDF
#1- exp(-.01 *theta) # probably coalesce before speciation time .02 = .07-0.06
print(paste("Theoretical: ",(1- exp(-.01 *theta)), ", Experimental: ", length(rowNum)/ length(firstCoal[,1]) ))
abs(length(rowNum)/ length(firstCoal[,1]) - (1- exp(-.01 *theta)))/ ((1- exp(-.01 *theta)))

# Coalesce times in ancestral pop given by pairs coalesce 
# in daughter pops (only one time)
rowNum <- intersect(which(firstCoal[, 3] < .07), which(secCoal[ , 3] < .07))
#mean(thirdCoal[rowNum, 3]) - 0.07
theta <- 1/(2 / .035)
print(paste("Theoretical: ", theta, ", Experimental: ", (mean(thirdCoal[rowNum, 3]) - 0.07)))
abs((mean(thirdCoal[rowNum, 3]) - 0.07)- theta) /theta

# Coalescent time in ancestral pop given only 0 and 1 coalesce 
# in daughter pop (2 times)
# or 
# Coalescent time in ancestral pop given only 2 and 3 coalesce 
# in daughter pop (2 times)
rowNum <- intersect(which(firstCoal[, 3] < .07), which(secCoal[ , 3] > .07))
# mean(secCoal[rowNum, 3]) - 0.07
theta <- 1/(2 *3 / .035)
print(paste("Theoretical: ", theta, ", Experimental: ", (mean(secCoal[rowNum, 3]) - 0.07)))
abs(theta- (mean(secCoal[rowNum, 3]) - 0.07)) /theta

#mean(thirdCoal[rowNum, 3] - secCoal[rowNum, 3]) 
theta <- 1/(2 / .035)
print(paste("Theoretical: ", theta, ", Experimental: ", (mean(thirdCoal[rowNum, 3] - secCoal[rowNum, 3]))))
abs(theta- (mean(thirdCoal[rowNum, 3] - secCoal[rowNum, 3]))) / theta

# coalescent times in ancestral pop given no coalescent events
# in ancestral pop
rowNum <- which(firstCoal[,3] > .07)
firstEvent <- firstCoal[rowNum, ]
#mean(firstEvent[,3]) -.07

theta <- 1/(4 *3 / .035)
print(paste("Theoretical: ", theta, ", Experimental: ", (mean(firstEvent[,3]) -.07)))
abs(theta- (mean(firstEvent[,3]) -.07)) /theta

secondEvent <- secCoal[rowNum, 3] -firstCoal[rowNum, 3]
#mean(secondEvent)
theta <- 1/(3 *2 / .035)
print(paste("Theoretical: ", theta, ", Experimental: ", mean(secondEvent)))
abs(theta- mean(secondEvent)) /theta

thirdEvent <- thirdCoal[rowNum, 3] -secCoal[rowNum, 3]
#mean(thirdEvent)
theta <- 1/(2 / .035)
print(paste("Theoretical: ", theta, ", Experimental: ", mean(thirdEvent)))
abs(theta- mean(thirdEvent))/theta

# Should I check the proportion of times we fall into each case?

# Testing 2 populations, 2 sequences per population 
# simple2.ctl
times <- read.csv("~/bpp/test/anna/testSimple2phase.txt", header =FALSE)
print("simple2phase.ctl: 2 populations, 1 sequences per population, unphased ")
print("This is percent error")
firstCoal <- times[seq(from =1, to = dim(times)[1] -2, by = 3), ]
secCoal <- times[seq(from =2, to = dim(times)[1]- 1, by = 3), ]
thirdCoal <- times[seq(from =3, to = dim(times)[1], by = 3), ]

# 0 and 1 coalesce in tip population
rowNum <- intersect(intersect(which(times[ ,1] == 0), which(times[,2] == 1)),
                    which(times[,3] < .07))
#mean(times[rowNum, 3])- 0.05 # subtract time of last sample

b <- .07-.05 # speciation time minus time of last sample
theta <-2/.04
abs((1/theta  -b *(exp(theta*b) -1)^-1) - (mean(times[rowNum, 3])- 0.05)) / (1/theta  -b *(exp(theta*b) -1)^-1)

#length(rowNum)/ length(firstCoal[,1])
#PDF
#1- exp(-.02 *theta) # probably coalesce before speciation time .02 = .07-0.05
abs(length(rowNum)/ length(firstCoal[,1]) - (1- exp(-.02 *theta) )) / (1- exp(-.02 *theta) )

# 2 and 3 coalesce in tip population 
rowNum <- intersect(intersect(which(times[ ,1] == 2), which(times[,2] == 3)),
                    which(times[,3] < .07))
#mean(times[rowNum, 3])- 0.06 # subtract time of last sample

b <- .07-.06 # speciation time minus time of last sample
theta <-2/.03 # k chose 2 is where the 2 is from
abs(1/theta  -b *(exp(theta*b) -1)^-1 - (mean(times[rowNum, 3])- 0.06)) / (1/theta  -b *(exp(theta*b) -1)^-1 )

#length(rowNum)/ length(firstCoal[,1])
#PDF
#1- exp(-.01 *theta) # probably coalesce before speciation time .02 = .07-0.06
abs(length(rowNum)/ length(firstCoal[,1]) - (1- exp(-.01 *theta)))/ ((1- exp(-.01 *theta)))

# Coalesce times in ancestral pop given by pairs coalesce 
# in daughter pops (only one time)
rowNum <- intersect(which(firstCoal[, 3] < .07), which(secCoal[ , 3] < .07))
#mean(thirdCoal[rowNum, 3]) - 0.07
theta <- 1/(2 / .035)
abs((mean(thirdCoal[rowNum, 3]) - 0.07)- theta) /theta

# Coalescent time in ancestral pop given only 0 and 1 coalesce 
# in daughter pop (2 times)
# or 
# Coalescent time in ancestral pop given only 2 and 3 coalesce 
# in daughter pop (2 times)
rowNum <- intersect(which(firstCoal[, 3] < .07), which(secCoal[ , 3] > .07))
# mean(secCoal[rowNum, 3]) - 0.07
theta <- 1/(2 *3 / .035)
abs(theta- (mean(secCoal[rowNum, 3]) - 0.07)) /theta

#mean(thirdCoal[rowNum, 3] - secCoal[rowNum, 3]) 
theta <- 1/(2 / .035)
abs(theta- (mean(thirdCoal[rowNum, 3] - secCoal[rowNum, 3]))) / theta

# coalescent times in ancestral pop given no coalescent events
# in ancestral pop
rowNum <- which(firstCoal[,3] > .07)
firstEvent <- firstCoal[rowNum, ]
#mean(firstEvent[,3]) -.07

theta <- 1/(4 *3 / .035)
abs(theta- (mean(firstEvent[,3]) -.07)) /theta

secondEvent <- secCoal[rowNum, 3] -firstCoal[rowNum, 3]
#mean(secondEvent)
theta <- 1/(3 *2 / .035)
abs(theta- mean(secondEvent)) /theta

thirdEvent <- thirdCoal[rowNum, 3] -secCoal[rowNum, 3]
#mean(thirdEvent)
theta <- 1/(2 / .035)
abs(theta- mean(thirdEvent))/theta



### Testing ancestral sampling
# Testing 2 populations, 1 sequence in one pop, 2 in the other, 
# one in the ancestral pop
# simple3.ctl
times <- read.csv("~/bpp/test/anna/testSimple3.txt", header =FALSE)
print("simple3.ctl: 2 populations, 1 sequence in one pop, 2 in the other, one in the ancestral pop")
firstCoal <- times[seq(from =1, to = dim(times)[1] -2, by = 3), ]
secCoal <- times[seq(from =2, to = dim(times)[1]- 1, by = 3), ]
thirdCoal <- times[seq(from =3, to = dim(times)[1], by = 3), ]

#####
# Coalescent event in tip population
rowNum <- which(firstCoal[,3] < .07)
#mean(firstCoal[rowNum, 3]) -0.05 # time of second sample 

b <- .07-.05 # speciation time minus time of last sample
theta <-2/.04 # k chose 2 is where the 2 is from
#1/theta  -b *(exp(theta*b) -1)^-1
print(paste("Theoretical: ", (1/theta  -b *(exp(theta*b) -1)^-1) , ", Experimental: ", (mean(firstCoal[rowNum, 3]) -0.05)))
abs((1/theta  -b *(exp(theta*b) -1)^-1) - (mean(firstCoal[rowNum, 3]) -0.05 )) / (1/theta  -b *(exp(theta*b) -1)^-1)

# Exactly 1 coalescent event before final sample added, conditional 
# on coalescent event tip population
rowNum2 <- which(secCoal[,3] < .074)
rowNumBoth <- intersect(rowNum, rowNum2)
#mean(secCoal[rowNumBoth, 3]) - .07

b <- .074-.07 # last sample time minus speciation time 
theta <-2/.035 # k chose 2 is where the 2 is from
#1/theta  -b *(exp(theta*b) -1)^-1
print(paste("Theoretical: ", (1/theta  -b *(exp(theta*b) -1)^-1) , ", Experimental: ", (mean(secCoal[rowNumBoth, 3]) - .07)))
abs((1/theta  -b *(exp(theta*b) -1)^-1)- (mean(secCoal[rowNumBoth, 3]) - .07)) / (1/theta  -b *(exp(theta*b) -1)^-1)

# No coalescent events before final sample added
rowNum  <- which(firstCoal[,3] < .07)
rowNum2 <- which(secCoal[,3] > 0.074)
rowNumBoth <- intersect(rowNum, rowNum2)
#mean(secCoal[rowNumBoth, 3]) - .074
#1/(3* 2 / .035)
print(paste("Theoretical: ", (1/(3* 2 / .035)) , ", Experimental: ", (mean(secCoal[rowNumBoth, 3]) - .074)))
abs((1/(3* 2 / .035)) - (mean(secCoal[rowNumBoth, 3]) - .074)) / (1/(3* 2 / .035))

#####
# No coalescent event in tip population
# First coalescent event is before final sample added
rowNum <- which(firstCoal[,3] > .07)
rowNum2 <- which(firstCoal[,3] < 0.074)
rowNumBoth <- intersect(rowNum, rowNum2)

#mean(firstCoal[rowNumBoth, 3]) - .07

b <- .074-.07 # last sample time minus speciation time 
theta <- 3*2/.035 # k chose 2 is where the 2 is from
#1/theta  -b *(exp(theta*b) -1)^-1
print(paste("Theoretical: ", (1/theta  -b *(exp(theta*b) -1)^-1) , ", Experimental: ", (mean(firstCoal[rowNumBoth, 3]) - .07)))
abs((1/theta  -b *(exp(theta*b) -1)^-1) - (mean(firstCoal[rowNumBoth, 3]) - .07)) / (1/theta  -b *(exp(theta*b) -1)^-1)

# start here again
# # 2nd coalescent event before final sample added
# rowNum <- which(firstCoal[,3] > .07)
# rowNum2 <- which(secCoal[,3] < 0.74)
# rowNumBoth <- intersect(rowNum, rowNum2)

# No coalescent events before final sample added
rowNum <- which(firstCoal[,3] > 0.074)
#mean(firstCoal[rowNum, 3]) - .074
#1/(4 *3 / .035)
print(paste("Theoretical: ", (1/(4 *3 / .035)) , ", Experimental: ", (mean(firstCoal[rowNum, 3]) - .074)))
abs((1/(4 *3 / .035)) -  (mean(firstCoal[rowNum, 3]) - .074)) / (1/(4 *3 / .035))

print(paste("Theoretical: ", (1/(2 *3 / .035)) , ", Experimental: ", mean(secCoal[rowNum, 3]- firstCoal[rowNum, 3])))
abs((1/(2 *3 / .035))-mean(secCoal[rowNum, 3]- firstCoal[rowNum, 3])) / (1/(2 *3 / .035))
#1/(2 *3 / .035)

#mean(thirdCoal[rowNum, 3]- secCoal[rowNum, 3]) 
#1/(2 / .035)
print(paste("Theoretical: ", (1/(2 / .035)) , ", Experimental: ", (mean(thirdCoal[rowNum, 3]- secCoal[rowNum, 3]))))
abs((1/(2 / .035)) - (mean(thirdCoal[rowNum, 3]- secCoal[rowNum, 3]))) / (1/(2 / .035))
## ANNA: check the proportions for ancestral sampling

#####
times <- read.csv("~/bpp/test/anna/testSimple4.txt", header =FALSE)
print("simple4.ctl: See control file for info")
firstCoal <- times[seq(from =1, to = dim(times)[1] -2, by = 4), ]
secCoal <- times[seq(from =2, to = dim(times)[1]- 1, by = 4), ]
thirdCoal <- times[seq(from =3, to = dim(times)[1], by = 4), ]
fourthCoal <- times[seq(from =4, to = dim(times)[1], by = 4), ]


# Coalesce before first speciation time 
rowNum <- which(firstCoal[,3] < .07)
#mean(firstCoal[rowNum, 3]) -0.05 # time of second sample 

b <- .07-.05 # speciation time minus time of last sample
theta <-2/.04 # k chose 2 is where the 2 is from
#1/theta  -b *(exp(theta*b) -1)^-1
print(paste("Theoretical: ", (1/theta  -b *(exp(theta*b) -1)^-1) , ", Experimental: ", (mean(firstCoal[rowNum, 3]) -0.05 )))
abs((1/theta  -b *(exp(theta*b) -1)^-1) - (mean(firstCoal[rowNum, 3]) -0.05 )) / (1/theta  -b *(exp(theta*b) -1)^-1)

# First coalesce between first speciation and sample added
rowNum <- intersect(which(firstCoal[,3] > .07), which(firstCoal[,3] < .074))
#mean(firstCoal[rowNum, 3]) -0.05 # time of second sample 

b <- .074-.07 # speciation time minus time of last sample
theta <-2 * 3/.035 # k chose 2 is where the 2 is from
#1/theta  -b *(exp(theta*b) -1)^-1
print(paste("Theoretical: ", (1/theta  -b *(exp(theta*b) -1)^-1) , ", Experimental: ", (mean(firstCoal[rowNum, 3]) -0.07 )))
abs((1/theta  -b *(exp(theta*b) -1)^-1) - (mean(firstCoal[rowNum, 3]) -0.07 )) / (1/theta  -b *(exp(theta*b) -1)^-1)


# First coalesce between sample added and second speciation 
rowNum <- intersect(which(firstCoal[,3] > .074), which(firstCoal[,3] < .08))
b <- .08-.074 # speciation time minus time of last sample
theta <- 3 * 4 /.035 # k chose 2 is where the 2 is from
#1/theta  -b *(exp(theta*b) -1)^-1
print(paste("Theoretical: ", (1/theta  -b *(exp(theta*b) -1)^-1) , ", Experimental: ",(mean(firstCoal[rowNum, 3]) -0.074 )))
abs((1/theta  -b *(exp(theta*b) -1)^-1) - (mean(firstCoal[rowNum, 3]) -0.074 )) / (1/theta  -b *(exp(theta*b) -1)^-1)

# First coalesce older than all speciation events
rowNum <- which(firstCoal[,3] > .08)
theta <- 5 * 4 /.05 # k chose 2 is where the 2 is from
#1/theta  -b *(exp(theta*b) -1)^-1
print(paste("Theoretical: ", (1/theta) , ", Experimental: ",(mean(firstCoal[rowNum, 3]) -0.08 )))
abs((1/theta) - (mean(firstCoal[rowNum, 3]) -0.08 )) / (1/theta)


# Coalesce before first speciation time, next coalescent time is in AB before next sequence is added
b <- .074- .07
theta <- 2 / .035
rowNum <- which(firstCoal[,3] < .07)
rowNum2 <- which(secCoal[,3] < .074)
rowsSecond <- intersect(rowNum2,rowNum) 
print(paste("Theoretical: ", (1/theta  -b *(exp(theta*b) -1)^-1) , ", Experimental: ", (mean(secCoal[rowsSecond,3]) - .07) ))
abs((1/theta  -b *(exp(theta*b) -1)^-1) - (mean(secCoal[rowsSecond,3]) - .07) ) / (1/theta  -b *(exp(theta*b) -1)^-1)

# Second coalescent event is in between sample added and ABC  speciation
b <- .08- .074
theta <- 2 * 3 / .035
rowNum2 <- intersect(intersect(which(secCoal[,3] > .074), which(secCoal[,3] < .08)), which(firstCoal[,3]<.074))
print(paste("Theoretical: ", (1/theta  -b *(exp(theta*b) -1)^-1) , ", Experimental: ", (mean(secCoal[rowNum2,3]) - .074) ))
abs((1/theta  -b *(exp(theta*b) -1)^-1) - (mean(secCoal[rowNum2,3]) - .074) ) / (1/theta  -b *(exp(theta*b) -1)^-1)

# Second coalescent event older than ABC 
rowNum <- intersect(which(secCoal[,3] >.08), which(firstCoal[,3] < .08))
theta <- 1/(3 *4 / .05)
print(paste("Theoretical: ", theta , ", Experimental: ", (mean(secCoal[rowNum, 3]) - .08) ))
abs(theta - (mean(secCoal[rowNum, 3]) - .08)) / (theta)

# Second coalescent event older than ABC, finding waiting time until third
rowNum <- which(secCoal[,3] >.08)
theta <- 1/(3 *2 / .05)
print(paste("Theoretical: ", theta , ", Experimental: ", (mean(thirdCoal[rowNum, 3] - secCoal[rowNum, 3]) )))
abs(theta - (mean(thirdCoal[rowNum, 3] - secCoal[rowNum, 3]))) / (theta)

# If 3rd coalescent time is before ABC
rowNum <- which(thirdCoal[,3] <.08)
#mean(fourthCoal[rowNum, 3]) - .08
theta <- 2  / .05
print(paste("Theoretical: ", 1/theta , ", Experimental: ", (mean(fourthCoal[rowNum, 3]) - .08) ))
abs(1/theta - (mean(fourthCoal[rowNum, 3]) - .08)) / (1/theta)

# If 3rd coalescent after ABC and 2nd before 
# Third coal
rowNum <- intersect(which(thirdCoal[,3] >.08), which(secCoal[,3] < .08))
theta <- 1/(3 *2 / .05)
print(paste("Theoretical: ", theta , ", Experimental: ", (mean(thirdCoal[rowNum, 3]) - .08) ))
abs(theta - (mean(thirdCoal[rowNum, 3]) - .08)) / (theta)

# If 3rd coalescent after ABC
# Fourth coal
rowNum <- which(thirdCoal[,3] >.08)
theta <- 2 / .05
print(paste("Theoretical: ", 1/theta , ", Experimental: ", (mean(fourthCoal[rowNum, 3] - thirdCoal[rowNum, 3])) ))
abs(1/theta - (mean(fourthCoal[rowNum, 3] - thirdCoal[rowNum, 3]))) / (1/theta)


# 
times <- read.csv("~/bpp/test/anna/testIM.txt", header =FALSE)
print("testIM.ctl, This is percent error")
firstCoal <- times[seq(from =1, to = dim(times)[1] -1, by = 2), ]
secCoal <- times[seq(from =2, to = dim(times)[1], by = 2), ]
mono <- intersect(which(firstCoal[, 1] != 2), which(firstCoal[, 2] != 2))
empirical <- length(mono)/ dim(firstCoal)[1]

M <- 0.3*4
# See Takahata and Slatkin 1990 for formula
theoretical <- (1 + 7*M/6 + M^2/3) / (1 + 5*M/2+M^2 )

A <- sqrt(1 + 4 * M ^2)
abs(empirical - theoretical) / theoretical
