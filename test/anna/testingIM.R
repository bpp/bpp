# Transition probability matrix for migration
prob.mat<-function(alpha,beta,t){
  exp.val<-exp(-(alpha+beta)*t)
  probabilities<-1/(alpha+beta)*matrix(c(beta+alpha*exp.val, beta-beta*exp.val, alpha-alpha*exp.val, alpha+beta*exp.val),nrow=2)
  return(probabilities)
}

# Multiplying the migration rates because this the the rate a lineage moves between populations
# (See coalescent calculations)
(P.mat<- prob.mat(alpha=0.2 /.004 *4 , beta=0.2 /.004 *4 , t=.005))
#####
# Check with no ancestral population (tau is old)
times <- read.csv("~/bpp/test/anna/mig/testIM1.txt", header =FALSE)
mean(times[,3]) -.005 # Subtract the sample time of sequence sampled further back in time

n <- 2 # Number of populations
M <- .2*4 # Multiply the rate used in bpp by 4 to account for the different definition of M
theta <- .004
# Expected coalescent time for two sequences sampled from different populations = theta +theta/(2M), where M = 4Nm
#(theta + theta  /(2 * M))

P.mat[1,2] * (theta) + P.mat[1,1] * (theta + theta  /(2 * M))

#####
# Check with less extreme tau
# ANNA FIGURE THIS OUT
times <- read.csv("~/bpp/test/anna/mig/testIM2.txt", header =FALSE)
theta_tip <- .04
theta_anc <- .03
mean(times[,3]) 
mean(times[,3]) /theta_tip#-.005 # Subtract the sample time of sequence sampled further back in time
n <- 2 # Number of populations
M <- .2*4

D <- (n*M + n - 1)^2 - (4*(n-1)*M)
lambda_1 <- (n * M + n - 1 - sqrt(D)) / (2 * (n-1))
lambda_2 <- (n * M + n - 1 + sqrt(D)) / (2 * (n-1))


a <- theta_anc/ theta_tip #1 # N_anc/N_tip
tau <- 1.4 / theta_tip #.2 * 2 /theta_tip
mean(times[,3] *2 /theta_tip)

E_T0 <-(n + ((n-1)/sqrt(D)) * ((lambda_2 -1)*(a-(1/lambda_1)) * exp(-lambda_1 * tau) + 
                        (1-lambda_1) * (a - (1/lambda_2)) * exp (-lambda_2 *tau))) 
E_T0 
E_T0 * theta_tip/2
                        

E_T1 <- (n + ((n-1) / M) + (n-1)/sqrt(D) * (lambda_2 *(a-(1/lambda_1)) * exp(-lambda_1 * tau) - 
                               lambda_1 *(a - 1/lambda_2) * exp (-lambda_2 *tau)))
E_T1
E_T1 * theta_tip/2

(P.mat<- prob.mat(alpha=0.2 /.04 *4 , beta=0.2 /.04 *4 , t=.005))
P.mat[1,2] * (E_T0 * theta_tip/2) + P.mat[1,1] * (E_T1 * theta_tip/2) +.005
mean(times[,3]) 

#####

times <- read.csv("~/bpp/test/anna/mig/testIM2.txt", header =FALSE)
#mean(times[,3]) 
#mean(times[,3]) /theta_tip-.005 # Subtract the sample time of sequence sampled further back in time
n <- 2 # Number of populations
M <- .2*4

D <- (n*M + n - 1)^2 - (4*(n-1)*M)
lambda_1 <- (n * M + n - 1 - sqrt(D)) / (2 * (n-1))
lambda_2 <- (n * M + n - 1 + sqrt(D)) / (2 * (n-1))
theta_tip <- .04
theta_anc <- .03

a <- theta_anc/ theta_tip #1 # N_anc/N_tip
tau <- (1.4 -.005) * 2 / theta_tip #.2 * 2 /theta_tip
mean(times[,3] *2 /theta_tip)

E_T0 <-(n + ((n-1)/sqrt(D)) * ((lambda_2 -1)*(a-(1/lambda_1)) * exp(-lambda_1 * tau) + 
                                 (1-lambda_1) * (a - (1/lambda_2)) * exp (-lambda_2 *tau))) 
E_T0 
E_T0 * theta_tip/2


E_T1 <- (n + ((n-1) / M) + (n-1)/sqrt(D) * (lambda_2 *(a-(1/lambda_1)) * exp(-lambda_1 * tau) - 
                                              lambda_1 *(a - 1/lambda_2) * exp (-lambda_2 *tau)))
E_T1
E_T1 * theta_tip/2


(P.mat<- prob.mat(alpha=0.2 /.04 *4 , beta=0.2 /.04 *4 , t=.005))
P.mat[1,2] * (E_T0 * theta_tip/2) + P.mat[1,1] * (E_T1 * theta_tip/2) +.005
mean(times[,3]) 


#####
times <- read.csv("~/bpp/test/anna/mig/testIM3.txt", header =FALSE)
#mean(times[,3]) 
#mean(times[,3]) /theta_tip-.005 # Subtract the sample time of sequence sampled further back in time
n <- 2 # Number of populations
M <- .2*4

D <- (n*M + n - 1)^2 - (4*(n-1)*M)
lambda_1 <- (n * M + n - 1 - sqrt(D)) / (2 * (n-1))
lambda_2 <- (n * M + n - 1 + sqrt(D)) / (2 * (n-1))
theta_tip <- .04
theta_anc <- .03

a <- theta_anc/ theta_tip #1 # N_anc/N_tip
tau <- (.1 -.005)* 2/ theta_tip #.2 * 2 /theta_tip
mean(times[,3] *2 /theta_tip)

E_T0 <-(n + ((n-1)/sqrt(D)) * ((lambda_2 -1)*(a-(1/lambda_1)) * exp(-lambda_1 * tau) + 
                                 (1-lambda_1) * (a - (1/lambda_2)) * exp (-lambda_2 *tau))) 
E_T0 
E_T0 * theta_tip/2


E_T1 <- (n + ((n-1) / M) + (n-1)/sqrt(D) * (lambda_2 *(a-(1/lambda_1)) * exp(-lambda_1 * tau) - 
                                              lambda_1 *(a - 1/lambda_2) * exp (-lambda_2 *tau)))
E_T1
E_T1 * theta_tip/2

(P.mat<- prob.mat(alpha=0.2 /.04 *4 , beta=0.2 /.04 *4 , t=.005))
P.mat[1,2] * (E_T0 * theta_tip/2) + P.mat[1,1] * (E_T1 * theta_tip/2) +.005
mean(times[,3]) 

#####
# Now test with larger number of sequences but considering just two sequences.
times <- read.csv("~/bpp/test/anna/mig/IM4/coal_a1_b2", header =FALSE, sep = '')
n <- 2 # Number of populations
M <- .2*4

D <- (n*M + n - 1)^2 - (4*(n-1)*M)
lambda_1 <- (n * M + n - 1 - sqrt(D)) / (2 * (n-1))
lambda_2 <- (n * M + n - 1 + sqrt(D)) / (2 * (n-1))
theta_tip <- .04
theta_anc <- .03

a <- theta_anc/ theta_tip 
tau <- (.05 -.006)* 2/ theta_tip 
mean(times[,2] *2 /theta_tip)

E_T0 <-(n + ((n-1)/sqrt(D)) * ((lambda_2 -1)*(a-(1/lambda_1)) * exp(-lambda_1 * tau) + 
                                 (1-lambda_1) * (a - (1/lambda_2)) * exp (-lambda_2 *tau))) 
E_T0 
E_T0 * theta_tip/2


E_T1 <- (n + ((n-1) / M) + (n-1)/sqrt(D) * (lambda_2 *(a-(1/lambda_1)) * exp(-lambda_1 * tau) - 
                                              lambda_1 *(a - 1/lambda_2) * exp (-lambda_2 *tau)))
E_T1
E_T1 * theta_tip/2

(P.mat<- prob.mat(alpha=0.2 /.04 *4 , beta=0.2 /.04 *4 , t=.006))
P.mat[1,2] * (E_T0 * theta_tip/2) + P.mat[1,1] * (E_T1 * theta_tip/2) + .006  
mean(times[,2]) 


#####
# Now test with larger number of sequences but considering just two sequences.
times <- read.csv("~/bpp/test/anna/mig/IM4/coal_a2_b2", header =FALSE, sep = '')
n <- 2 # Number of populations
M <- .2*4

D <- (n*M + n - 1)^2 - (4*(n-1)*M)
lambda_1 <- (n * M + n - 1 - sqrt(D)) / (2 * (n-1))
lambda_2 <- (n * M + n - 1 + sqrt(D)) / (2 * (n-1))
theta_tip <- .04
theta_anc <- .03

a <- theta_anc/ theta_tip 
tau <- (.05 -.006)* 2/ theta_tip 
mean(times[,2] *2 /theta_tip)

E_T0 <-(n + ((n-1)/sqrt(D)) * ((lambda_2 -1)*(a-(1/lambda_1)) * exp(-lambda_1 * tau) + 
                                 (1-lambda_1) * (a - (1/lambda_2)) * exp (-lambda_2 *tau))) 
E_T0 
E_T0 * theta_tip/2


E_T1 <- (n + ((n-1) / M) + (n-1)/sqrt(D) * (lambda_2 *(a-(1/lambda_1)) * exp(-lambda_1 * tau) - 
                                              lambda_1 *(a - 1/lambda_2) * exp (-lambda_2 *tau)))
E_T1
E_T1 * theta_tip/2

(P.mat<- prob.mat(alpha=0.2 /.04 *4 , beta=0.2 /.04 *4 , t=.006 - .001))
P.mat[1,2] * (E_T0 * theta_tip/2) + P.mat[1,1] * (E_T1 * theta_tip/2) + .006 - .001 
mean(times[,2]) 

