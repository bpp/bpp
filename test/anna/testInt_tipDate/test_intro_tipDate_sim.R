setwd("~/bpp/test/anna/testInt_tipDate")

# All are with a two tip (A and B) tree, with bidirectional introgression
##### 
# Two samples from A
times <- read.table("ex_1/coal_time", header = FALSE)[,1]

# Option 1- coalesce before hybridization
t_2 <- 0.03
h_time <- 0.035
theta <- 2/0.01
rows <- which(times < h_time)
# Mean waiting time conditional on occuring before h_time
mean(times[rows]) - t_2
b <- h_time - t_2

# Mean of truncated exponential is  1/theta  -b *(exp(theta*b) -1)^-1
ex <- 1/theta - b * (exp(theta * b) -1 )^-1
ex

# Probability of coalescing before hybridization 
mean(times < h_time) -
(1 - exp(-theta * (h_time-t_2)))

# Option 2- coalesce after hybridization before root
t_root <- 0.045
phi <- 0.1
p_both_stay <- (1-phi)^2
p_both_leave <- phi ^2
p_one_leave <- 2 * (1-phi) * phi

b <- t_root - h_time
ex <- 1/theta - b * (exp(theta * b) -1 )^-1
ex
rows <- intersect(which(times < t_root), which(times > h_time))
mean(times[rows]) - h_time

# probability of coalesce is 
# sum of probability in same pop * prob coalesce, 
# conditional on greater than h_time
p_same <- p_both_stay + p_both_leave
# Prob coalesce given time is older than hyrbidization time
row_htime <- which(times > h_time)
mean(times[row_htime] < t_root)
p_same * (1 - exp(-theta * (t_root-h_time)))


# Option 3- coalesce after root
rows <- which (times > t_root)
1/theta
mean(times[rows]) - t_root

#####
# Two samples from B
times <- read.table("ex_2/coal_time", header = FALSE)[,1]

# Option 1- coalesce before hybridization
t_2 <- 0.025
h_time <- 0.035
theta <- 2/0.01
rows <- which(times < h_time)
# Mean waiting time conditional on occuring before h_time
mean(times[rows]) - t_2
b <- h_time - t_2

# Mean of truncated exponential is  1/theta  -b *(exp(theta*b) -1)^-1
ex <- 1/theta - b * (exp(theta * b) -1 )^-1
ex

# Probability of coalescing before hybridization 
mean(times < h_time) -
  (1 - exp(-theta * (h_time-t_2)))

# Option 2- coalesce after hybridization before root
t_root <- 0.045
phi <- 0.3
p_both_stay <- (1-phi)^2
p_both_leave <- phi ^2
p_one_leave <- 2 * (1-phi) * phi

b <- t_root - h_time
ex <- 1/theta - b * (exp(theta * b) -1 )^-1
ex
rows <- intersect(which(times < t_root), which(times > h_time))
mean(times[rows]) - h_time

# probability of coalesce is 
# sum of probability in same pop * prob coalesce, 
# conditional on greater than h_time
p_same <- p_both_stay + p_both_leave
# Prob coalesce given time is older than hyrbidization time
row_htime <- which(times > h_time)
mean(times[row_htime] < t_root)
p_same * (1 - exp(-theta * (t_root-h_time)))


# Option 3- coalesce after root
rows <- which (times > t_root)
1/theta
mean(times[rows]) - t_root


#####
# Two samples from B, unequal theta after introgression
times <- read.table("ex_3/coal_time", header = FALSE)[,1]
pop_X <- read.table("ex_3/debug_X", header = FALSE)[,1]
pop_Y <- read.table("ex_3/debug_Y", header = FALSE)[,1]
# Option 1- coalesce after hybridization before root
h_time <- 0.035
theta_1 <- 2/0.01
theta_2 <- 2/0.035


t_root <- 0.045
phi <- 0.3

p_both_stay <- (1-phi) * (1-phi)
p_both_leave <- phi * phi
p_diff <- 2 * (1-phi) * phi
p_same <- p_both_stay + p_both_leave

# Prob coalesce given time is older than hyrbidization time
rows <- which(times > h_time)
mean(times[rows] < t_root)
b <- t_root-h_time
p_both_leave * (1 - exp(-theta_1 * b)) + p_both_stay * (1 - exp(-theta_2 * b)) 


# Note that the mean time of coalescence given coalescence occurs 
# before the root (and after introgression) can not be calculated simply 
# taking the probability of both lineages being in X times mean in X + prob 
# in Y * mean in Y / prob in same pop
# This is conditioning on coalescence changes the probability of being in each 
# of the populations from phi
# Thus we need to separate based on which lineages are in which pop- can get this
# from the debugging info

b <- t_root - h_time
ex_1 <- 1/theta_1 - b * (exp(theta_1 * b) -1 )^-1
ex_2 <- 1/theta_2 - b * (exp(theta_2 * b) -1 )^-1
ex_1
ex_2
rows <- intersect(which(times < t_root), which(times > h_time))

x_row <- which(pop_X ==2)
y_row <- which(pop_Y ==2)


mean(times[intersect(x_row, rows)]) -h_time
mean(times[intersect(y_row, rows)]) -h_time


# Option 2- coalesce after root
rows <- which (times > t_root)
theta <- 2/0.01
1/theta
mean(times[rows]) - t_root

#####
# One sample from A, one from B, unequal theta after introgression
times <- read.table("ex_4/coal_time", header = FALSE)[,1]
pop_X <- read.table("ex_4/debug_X", header = FALSE)[,1]
pop_Y <- read.table("ex_4/debug_Y", header = FALSE)[,1]

# Option 1- coalesce after hybridization before root
h_time <- 0.035
t_root <- 0.045
theta_1 <- 2/0.01
theta_2 <- 2/0.035
phi_1 <- 0.1
phi_2 <- 0.3
p_both_stay <- (1-phi_1) * (1-phi_2)
p_both_leave <- phi_1 * phi_2
p_both_a <- (1-phi_1) * phi_2
p_both_b <- phi_1 * (1-phi_2)


# Prob coalesce given time is younger than root
mean(times < t_root)
b <- t_root - h_time
p_both_a * (1 - exp(-theta_1 * b)) + p_both_b* (1 - exp(-theta_2 * b)) 

#
ex_1 <- 1/theta_1 - b * (exp(theta_1 * b) -1 )^-1
ex_2 <- 1/theta_2 - b * (exp(theta_2 * b) -1 )^-1
ex_1
ex_2
rows <- intersect(which(times < t_root), which(times > h_time))

x_row <- which(pop_X ==2)
y_row <- which(pop_Y ==2)

mean(times[intersect(x_row, rows)]) -h_time
mean(times[intersect(y_row, rows)]) -h_time
  
# Option 2- coalesce after root
rows <- which (times > t_root)
theta <- 2/0.01
1/theta
mean(times[rows]) - t_root

#####
# Two samples from A, one from B, unequal theta after introgression
# Coal times from the A samples
times <- read.table("ex_5/coal_time", header = FALSE, sep = ",")
pop_X <- read.table("ex_5/debug_X", header = FALSE)[,1] # Number of lineages in X
pop_Y <- read.table("ex_5/debug_Y", header = FALSE)[,1] # Number of lineages in Y 
firstCoal <- times[seq(from = 1, to = dim(times)[1] -1, by = 2), ]
secCoal <- times[seq(from = 2, to = dim(times)[1], by = 2), ]
#first_event_aa <- which(firstCoal$V1 + firstCoal$V2 == 1) 
#second_event_aa <- which(firstCoal$V1 + firstCoal$V2 != 1) 
#times <- c(firstCoal[first_event_aa, 3], secCoal[second_event_aa, 3])

# Option 1- coalesce before hybridization
t_2 <- 0.03
h_time <- 0.035
theta <- 2/0.01
rows <- which(firstCoal < h_time)
# Mean waiting time conditional on occurring before h_time
mean(firstCoal[rows]) - t_2
b <- h_time - t_2

# Mean of truncated exponential is  1/theta  -b *(exp(theta*b) -1)^-1
ex <- 1/theta - b * (exp(theta * b) -1 )^-1
ex

# Probability of coalescing before hybridization 
mean(firstCoal < h_time) -
  (1 - exp(-theta * (h_time-t_2)))
prob_no_coal <- 1- (1 - exp(-theta * (h_time-t_2)))

# Option 2- coalesce after hybridization before root
t_root <- 0.045
phi_1 <- 0.1
phi_2 <- 0.3
p_both_A_stay <- (1-phi_1)^2
p_both_A_leave <- phi_1 ^2
p_one_A_leave <- 2 * (1-phi_1) * phi_1

# Both in A stay, one in B comes
A_3 <- p_both_A_stay * phi_2 
mean(pop_X[which(firstCoal > h_time)] == 3)
A_3
# Both in A stay, B doesnt or 
A_2 <- p_both_A_stay * (1-phi_2) + p_one_A_leave * phi_2  #Given no coal
mean(pop_X[which(firstCoal > h_time)] == 2)
A_2
A_1 <- p_one_A_leave * (1-phi_2) + p_both_A_leave * phi_2
mean(pop_X[which(firstCoal > h_time)] == 1)
A_1
A_0 <- p_both_A_leave * (1-phi_2)
A_0
mean(pop_X[which(firstCoal > h_time)] == 0)


t <- t_root - h_time
prob <- A_3 * (1 - exp(- 2 * 3/(0.01) * t)) + 
        A_2 * (1 - exp(- 2/0.01 * t)) + 
        A_1 * (1 - exp(- 2/0.035  * t)) + 
        A_0 * (1 - exp(- 2 * 3/(0.035) * t))
 
theta_1 <- 2/0.01
theta_2 <- 2/0.035

#Given no coalescence before hybridization, the probability of coalescing before root 
mean( (firstCoal[which(firstCoal > h_time)] < t_root)) 
prob 

# conditional on coalesence before hybridization, probability of coalescing before root 

# Both in A stay, B doesnt or 
A_2 <- (1-phi_1) * phi_2 
mean(pop_X[which(firstCoal < h_time)] == 2)
A_2
A_1 <- (1-phi_1) * (1-phi_2) + phi_1 * phi_2 
mean(pop_X[which(firstCoal < h_time)] == 1)
A_1
A_0 <- phi_1 * (1-phi_2)
A_0
mean(pop_X[which(firstCoal > h_time)] == 0)

# ANNA: you are going to have to know which lineages move to get distribution of times-
# leave this for now

# Option 3- coalesce after root
rows <- which (firstCoal > t_root)
1/theta/3
(mean(firstCoal[rows]) - t_root)
mean((secCoal[rows]) - firstCoal[rows])

rows <- which ((firstCoal < t_root & secCoal > t_root) == TRUE)
1/theta
mean(secCoal[rows]) - t_root

# # Coal times from A sample1 and B sample
# first_event_aa <- which(firstCoal$V1 + firstCoal$V2 == 2) 
# second_event_aa <- which(firstCoal$V1 + firstCoal$V2 != 2) 
# times <- c(firstCoal[first_event_aa, 3], secCoal[second_event_aa, 3])
# 
# # Option 1- coalesce after hybridization before root
# h_time <- 0.035
# t_root <- 0.045
# theta_1 <- 2/0.01
# theta_2 <- 2/0.035
# phi_1 <- 0.1
# phi_2 <- 0.3
# p_both_stay <- (1-phi_1) * (1-phi_2)
# p_both_leave <- phi_1 * phi_2
# p_both_a <- (1-phi_1) * phi_2
# p_both_b <- phi_1 * (1-phi_2)
# 
# 
# # Prob coalesce given time is younger than root
# mean(times < t_root)
# b <- t_root - h_time
# p_both_a * (1 - exp(-theta_1 * b)) + p_both_b* (1 - exp(-theta_2 * b)) 
# 
# ex_1 <- 1/theta_1 - b * (exp(theta_1 * b) -1 )^-1
# ex_2 <- 1/theta_2 - b * (exp(theta_2 * b) -1 )^-1
# rows <-  which(times < t_root)
# 
# #ANNA: Again you are going to need to know which lineages are where to get the
# # mean coalescence time (or do some conditional math)
# 
# 
# # Option 2- coalesce after root
# rows <- which (times > t_root)
# 1/theta
# mean(times[rows]) - t_root
# 
# # Coal times from A sample2 and B sample
# first_event_aa <- which(firstCoal$V1 + firstCoal$V2 == 3) 
# second_event_aa <- which(firstCoal$V1 + firstCoal$V2 != 3) 
# times <- c(firstCoal[first_event_aa, 3], secCoal[second_event_aa, 3])
# 
# # Option 1- coalesce after hybridization before root
# h_time <- 0.035
# t_root <- 0.045
# theta_1 <- 2/0.01
# theta_2 <- 2/0.035
# phi_1 <- 0.1
# phi_2 <- 0.3
# p_both_stay <- (1-phi_1) * (1-phi_2)
# p_both_leave <- phi_1 * phi_2
# p_both_a <- (1-phi_1) * phi_2
# p_both_b <- phi_1 * (1-phi_2)
# 
# 
# # Prob coalesce given time is younger than root
# mean(times < t_root)
# b <- t_root - h_time
# p_both_a * (1 - exp(-theta_1 * b)) + p_both_b* (1 - exp(-theta_2 * b)) 
# 
# ex_1 <- 1/theta_1 - b * (exp(theta_1 * b) -1 )^-1
# ex_2 <- 1/theta_2 - b * (exp(theta_2 * b) -1 )^-1
# rows <-  which(times < t_root)
# mean(times[rows]) - h_time 
# 
# (p_both_a * (ex_1) + p_both_b *(ex_2) ) / (p_both_a + p_both_b)
# 
# 
# # Option 2- coalesce after root
# rows <- which (times > t_root)
# 1/theta
# mean(times[rows]) - t_root
# 
#####
# Two samples from A and 2 from B
# Samples from A and A
times <- read.table("ex_6/A_A", header = FALSE, sep = " ")[,2]
t_a1 <- 0.025
times <- times + t_a1

# Option 1- coalesce before hybridization
t_2 <- 0.03
h_time <- 0.035
theta <- 2/0.008
rows <- which(times < h_time)
# Mean waiting time conditional on occuring before h_time
mean(times[rows]) - t_2
b <- h_time - t_2

# Mean of truncated exponential is  1/theta  -b *(exp(theta*b) -1)^-1
ex <- 1/theta - b * (exp(theta * b) -1 )^-1
ex

# Probability of coalescing before hybridization 
mean(times < h_time) -
  (1 - exp(-theta * (h_time-t_2)))

# Option 2- coalesce after hybridization before root
t_root <- 0.045
phi <- 0.1
theta_1 <- 2/0.01
theta_2 <- 2/0.035
p_both_stay <- (1-phi)^2
p_both_leave <- phi ^2
p_one_leave <- 2 * (1-phi) * phi

b <- t_root - h_time

# Prob coalesce given time is older than hyrbidization time
row_htime <- which(times > h_time)
mean(times[row_htime] < t_root)
coal <- p_both_stay * (1 - exp(-theta_1 * b)) + p_both_leave * (1 - exp(-theta_2 * b))
coal
no_coal <- p_both_stay * ( exp(-theta_1 * b)) + p_both_leave * ( exp(-theta_2 * b)) + p_one_leave



# Option 3- coalesce after root
rows <- which (times > t_root)
theta_r <- 2/0.005
1/theta_r
mean(times[rows]) - t_root

# Samples from B and B 
times <- read.table("ex_6/B_B", header = FALSE, sep = " ")[,2]
t_b1 <- 0.02
times <- times + t_b1

# Option 1- coalesce before hybridization
t_2 <- 0.027
h_time <- 0.035
theta <- 2/0.015
rows <- which(times < h_time)
# Mean waiting time conditional on occuring before h_time
mean(times[rows]) - t_2
b <- h_time - t_2

# Mean of truncated exponential is  1/theta  -b *(exp(theta*b) -1)^-1
ex <- 1/theta - b * (exp(theta * b) -1 )^-1
ex

# Probability of coalescing before hybridization 
mean(times < h_time) -
  (1 - exp(-theta * (h_time-t_2)))

# Option 2- coalesce after hybridization before root
t_root <- 0.045
phi <- 0.3
theta_2 <- 2/0.01
theta_1 <- 2/0.035
p_both_stay <- (1-phi)^2
p_both_leave <- phi ^2
p_one_leave <- 2 * (1-phi) * phi

b <- t_root - h_time
ex_1 <- 1/theta_1 - b * (exp(theta_1 * b) -1 )^-1
ex_2 <- 1/theta_2 - b * (exp(theta_2 * b) -1 )^-1


# Prob coalesce given time is older than hyrbidization time
row_htime <- which(times > h_time)
mean(times[row_htime] < t_root)
coal <- p_both_stay * (1 - exp(-theta_1 * b)) + p_both_leave * (1 - exp(-theta_2 * b))
coal
no_coal <- p_both_stay * ( exp(-theta_1 * b)) + p_both_leave * ( exp(-theta_2 * b)) + p_one_leave



# Option 3- coalesce after root
rows <- which (times > t_root)
theta_r <- 2/0.005
1/theta_r
mean(times[rows]) - t_root
