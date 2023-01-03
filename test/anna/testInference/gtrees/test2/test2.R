
# With simple2.ctl
times <- read.csv("~/bpp/test/anna/testInference/gtrees/test2/coal2_a1_c1", header=FALSE)
theta <-  .006/2
theta 
mean(times[,1]) + .001 -.016 -theta

times <- read.csv("~/bpp/test/anna/testInference/gtrees/test2/coal2_a1_a3", header=FALSE)
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
mean(coalAB)-.008+.001
1/theta  -b *(exp(theta*b) -1)^-1 
mean(coalAB)-.008+.001 - (1/theta  -b *(exp(theta*b) -1)^-1 )

theta <- .006/2
.016-.001 + theta
mean(coalABC)


times <- read.csv("~/bpp/test/anna/testInference/gtrees/test2/coal2_a1_a2", header=FALSE)
theta <-  .007/2
b <- .008- .0015
coalA <- times[which(times[,1] < .008 - .001),]
coalAB <- times[intersect(which(times[,1] > .008 - .001), which(times[,1] < .016 - .001)),]
coalABC <- times[which(times[,1] > .016 - .001),]
length(coalA)/length(times[,1])
theta <- 1/theta
mean(coalA)-(.0015-.001)
1/theta  -b *(exp(theta*b) -1)^-1 
mean(coalA)-(.0015-.001) - (1/theta  -b *(exp(theta*b) -1)^-1)

theta <- 2/.008
b <- .016-.008
mean(coalAB)-.008+.001
1/theta  -b *(exp(theta*b) -1)^-1 
mean(coalAB)-.008+.001 - (1/theta  -b *(exp(theta*b) -1)^-1) 

theta <- .006/2
.016-.001 + theta
mean(coalABC)


times <- read.csv("~/bpp/test/anna/testInference/gtrees/test2/coal2_a1_c5", header=FALSE)
theta <-  .006/2
theta 
mean(times[,1]) - (.016- .001)
theta - ( mean(times[,1]) - (.016- .001))

times <- read.csv("~/bpp/test/anna/testInference/gtrees/test2/coal2_a3_b2", header=FALSE)
theta <-  .008/2
b <- .016-.008
theta 
coalAB <- times[which(times[,1] < .016 - .0022),]
coalABC <- times[which(times[,1] > .016 - .0022),]

theta <- 1/theta
1/theta  -b *(exp(theta*b) -1)^-1 
mean(coalAB) - (.008- .0022) - (1/theta  -b *(exp(theta*b) -1)^-1 )

theta <- .006/2
.016-.0022 + theta
mean(coalABC)


times <- read.csv("~/bpp/test/anna/testInference/gtrees/test2/coal3_b1_b2", header=FALSE)
theta <-  .004/2
b <- .008- .0023
coalB <- times[which(times[,1] < .008 - .001),]
coalAB <- times[intersect(which(times[,1] > .008 - .001), which(times[,1] < .016 - .001)),]
coalABC <- times[which(times[,1] > .016 - .001),]
length(coalB)/length(times[,1])
theta <- 1/theta
mean(coalB)-(.0023-.001) - (1/theta  -b *(exp(theta*b) -1)^-1 )
1/theta  -b *(exp(theta*b) -1)^-1 

theta <- 2/.008
b <- .016-.008
mean(coalAB)-.008+.001
1/theta  -b *(exp(theta*b) -1)^-1 
mean(coalAB)-.008+.001 - (1/theta  -b *(exp(theta*b) -1)^-1 )

theta <- .006/2
.016-.001 + theta
mean(coalABC)
.016-.001 + theta - mean(coalABC)


times <- read.csv("~/bpp/test/anna/testInference/gtrees/test2/coal3_b1_b3", header=FALSE)
theta <-  .004/2
b <- .008- .003
coalB <- times[which(times[,1] < .008 - .001),]
coalAB <- times[intersect(which(times[,1] > .008 - .001), which(times[,1] < .016 - .001)),]
coalABC <- times[which(times[,1] > .016 - .001),]
length(coalB)/length(times[,1])
theta <- 1/theta
mean(coalB)-(.003-.001) - (1/theta  -b *(exp(theta*b) -1)^-1 )
1/theta  -b *(exp(theta*b) -1)^-1 

theta <- 2/.008
b <- .016-.008
mean(coalAB)-.008+.001
mean(coalAB)-.008+.001 - (1/theta  -b *(exp(theta*b) -1)^-1 )
1/theta  -b *(exp(theta*b) -1)^-1 

theta <- .006/2
.016-.001 + theta
mean(coalABC)- (.016-.001 + theta)
mean(coalABC)
