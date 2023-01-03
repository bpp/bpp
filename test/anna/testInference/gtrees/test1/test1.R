# With simple.ctl
times <- read.csv("~/bpp/test/anna/testInference/gtrees/test1/coal_A^a1_C^c1", header=FALSE)
theta <-  .006/2
theta 
mean(times[,1]) + .001 -.2

times <- read.csv("~/bpp/test/anna/testInference/gtrees/test1/coal_A^a1_A^a3", header=FALSE)
theta <-  .004/2
theta 
mean(times[,1]) - (.0022- .001)

times <- read.csv("~/bpp/test/anna/testInference/gtrees/test1/coal_A^a1_A^a2", header=FALSE)
theta <-  .004/2
theta 
mean(times[,1]) - (.0015- .001)

times <- read.csv("~/bpp/test/anna/testInference/gtrees/test1/coal_A^a1_C^c5", header=FALSE)
theta <-  .006/2
theta 
mean(times[,1]) - (.2- .001)

times <- read.csv("~/bpp/test/anna/testInference/gtrees/test1/coal_A^a3_B^b2", header=FALSE)
theta <-  .0035/2
theta 
mean(times[,1]) - (.1- .0022)
