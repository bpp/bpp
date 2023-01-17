print("a1 c1\n")
times <- read.csv("~/bpp/test/anna/testInference/gtrees/test3/coal2_a1_c1", header=FALSE)
theta <-  .006/2
theta 
mean(times[,1]) + .001 -.016 -theta

print("a1 a3")
times <- read.csv("~/bpp/test/anna/testInference/gtrees/test3/coal2_a1_a3", header=FALSE)
theta <-  .007/2
b <- .008-.0022
coalA <- times[which(times[,1] < .008 - .001),]
coalAB <- times[intersect(which(times[,1] > .008 - .001), which(times[,1] < .016 - .001)),]
coalABC <- times[which(times[,1] > .016 - .001),]
theta <-  2/.007
1/theta  -b *(exp(theta*b) -1)^-1 
mean(coalA) - .0022+.001 - (1/theta  -b *(exp(theta*b) -1)^-1 )

theta <- 2/.008
b <- .016-.008
#mean(coalAB)-.008+.001
1/theta  -b *(exp(theta*b) -1)^-1 
mean(coalAB)-.008+.001 - (1/theta  -b *(exp(theta*b) -1)^-1 )

# Check?
theta <- .006/2
theta
mean(coalABC) - (.016-.001 + theta)

print("a1 a2")
times <- read.csv("~/bpp/test/anna/testInference/gtrees/test3/coal2_a1_a2", header=FALSE)
theta <-  .007/2
b <- .008- .0015
coalA <- times[which(times[,1] < .008 - .001),]
coalAB <- times[intersect(which(times[,1] > .008 - .001), which(times[,1] < .016 - .001)),]
coalABC <- times[which(times[,1] > .016 - .001),]
#length(coalA)/length(times[,1])
theta <- 1/theta
#mean(coalA)-(.0015-.001)
1/theta  -b *(exp(theta*b) -1)^-1 
mean(coalA)-(.0015-.001) - (1/theta  -b *(exp(theta*b) -1)^-1)

theta <- 2/.008
b <- .016-.008
#mean(coalAB)-.008+.001
1/theta  -b *(exp(theta*b) -1)^-1 
mean(coalAB)-.008+.001 - (1/theta  -b *(exp(theta*b) -1)^-1) 

theta <- .006/2
theta
mean(coalABC) - (.016-.001 + theta)


print("a1 c5")
times <- read.csv("~/bpp/test/anna/testInference/gtrees/test3/coal2_a1_c5", header=FALSE)
theta <-  .006/2
theta 
#mean(times[,1]) - (.016- .001)
-(theta - ( mean(times[,1]) - (.016- .001)))

print("a3 b2")
times <- read.csv("~/bpp/test/anna/testInference/gtrees/test3/coal2_a3_b2", header=FALSE)
theta <-  .008/2
b <- .016-.008
#theta 
coalAB <- times[which(times[,1] < .016 - .0022),]
coalABC <- times[which(times[,1] > .016 - .0022),]

theta <- 1/theta
1/theta  -b *(exp(theta*b) -1)^-1 
mean(coalAB) - (.008- .0022) - (1/theta  -b *(exp(theta*b) -1)^-1 )

theta <- .006/2
#.016-.0022 + theta
theta
mean(coalABC) - (.016-.0022 + theta)


print("b1 b2")
times <- read.csv("~/bpp/test/anna/testInference/gtrees/test3/coal3_b1_b2", header=FALSE)
theta <-  .004/2
b <- .008- .0023
coalB <- times[which(times[,1] < .008 - .001),]
coalAB <- times[intersect(which(times[,1] > .008 - .001), which(times[,1] < .016 - .001)),]
coalABC <- times[which(times[,1] > .016 - .001),]
#length(coalB)/length(times[,1])
theta <- 1/theta
1/theta  -b *(exp(theta*b) -1)^-1 
mean(coalB)-(.0023-.001) - (1/theta  -b *(exp(theta*b) -1)^-1 )

theta <- 2/.008
b <- .016-.008
#mean(coalAB)-.008+.001
1/theta  -b *(exp(theta*b) -1)^-1 
mean(coalAB)-.008+.001 - (1/theta  -b *(exp(theta*b) -1)^-1 )

theta <- .006/2
#.016-.001 + theta
theta
#mean(coalABC)
-(.016-.001 + theta - mean(coalABC))


print("b1 b3")
times <- read.csv("~/bpp/test/anna/testInference/gtrees/test3/coal3_b1_b3", header=FALSE)
theta <-  .004/2
b <- .008- .003
coalB <- times[which(times[,1] < .008 - .001),]
coalAB <- times[intersect(which(times[,1] > .008 - .001), which(times[,1] < .016 - .001)),]
coalABC <- times[which(times[,1] > .016 - .001),]
#length(coalB)/length(times[,1])
theta <- 1/theta
1/theta  -b *(exp(theta*b) -1)^-1 
mean(coalB)-(.003-.001) - (1/theta  -b *(exp(theta*b) -1)^-1 )

theta <- 2/.008
b <- .016-.008
#mean(coalAB)-.008+.001
1/theta  -b *(exp(theta*b) -1)^-1 
mean(coalAB)-.008+.001 - (1/theta  -b *(exp(theta*b) -1)^-1 )

theta <- .006/2
theta
mean(coalABC)- (.016-.001 + theta)


print("a1 ab1")
times <- read.csv("~/bpp/test/anna/testInference/gtrees/test3/coal4_a1_ab1", header=FALSE)
theta <-  .008/2
b <- .016-.009
coalAB <- times[which(times[,1] < .016 - .001),]
coalABC <- times[which(times[,1] > .016 - .001),]
theta <- 1/theta
1/theta  -b *(exp(theta*b) -1)^-1 
mean(coalAB)-(.009-.001) - (1/theta  -b *(exp(theta*b) -1)^-1 )

theta <- .006/2
theta
mean(coalABC)- (.016-.001 + theta)

print("a4 ab2")
times <- read.csv("~/bpp/test/anna/testInference/gtrees/test3/coal4_a4_ab2", header=FALSE)
theta <-  .008/2
b <- .016-.0095
coalAB <- times[which(times[,1] < .016 - .0025),]
coalABC <- times[which(times[,1] > .016 - .0025),]
theta <- 1/theta
1/theta  -b *(exp(theta*b) -1)^-1 
mean(coalAB)-(.0095-.0025) - (1/theta  -b *(exp(theta*b) -1)^-1 )

theta <- .006/2
theta
mean(coalABC)- (.016-.0025 + theta)

print("b3 ab3")
times <- read.csv("~/bpp/test/anna/testInference/gtrees/test3/coal4_b3_ab3", header=FALSE)
theta <-  .008/2
b <- .016-.012
coalAB <- times[which(times[,1] < .016 - .003),]
coalABC <- times[which(times[,1] > .016 - .003),]
theta <- 1/theta
1/theta  -b *(exp(theta*b) -1)^-1 
mean(coalAB)-(.012-.003) - (1/theta  -b *(exp(theta*b) -1)^-1 )

theta <- .006/2
theta
mean(coalABC)- (.016-.003 + theta)
