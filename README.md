matlab-estimate-flat-multinomial
================================

This is a Matlab function that estimates the number of categories of a multinomial distribution

Estimates the number of categories of a uniform multinomial distribution

ESTIMATEFLATMULTINOMIAL(N_trial, N_different, [N_category_max])
  solves the Bayesian estimation problem of a uniform multinomial
  distribution of unknown dimension (i.e. number of categories, N_category)
  using the number of trials (N_trial) and the number of different
  outcomes (N_different) as input. Assumes a flat prior for N_category.

  The sampling distribution is 
      P(n_1, n_2, ... n_k) = [n! / (n_1! * n_2! * ... * n_k!)] * (1/k)^n
          where                       
          k  = N_category
          n  = sum(n_j) = N_trial
          n_j = the number of outcomes from category j
          sum( sgn(n_j) ) = N_different
  
  The posterior distribution is
      P(N_category = k) = A * k^(-n) * k! / ((k - N_different)!)
          where A is a normalization constant

  [N_category_dist, avg, variance, flag] = ESTIMATEFLATMULTINOMIAL 
      has four ouputs:
  N_category_dist =  N_category_max by 3 matrix:
      col 1: integers from 1 to N_category_max
      col 2: PDF of N_category
      col 3: (1 - CDF) of N_category
  avg = expectation value of N_category (=NaN if non-existent)
  variance = variance of N_category (=NaN if non-existent)
  flag = 0, posterior is normalized (on [1,Inf) )
       = 1, posterior is non-normalizable

  N_category_max is the largest possible value for N_category.
  If specified, the posterior is normalized over [1,N_category_max].

  If N_category_max is unspecified, 
  it is considered to be effectively infinite, and
  the normalization is carreid out over the entire [0,Inf) domain.
  Similarly, avg and variance are calculated over [0,Inf), if they exist.
  The maximal dimension of N_category_dist is determined
  by the criteria that the probability in the tail of the
  posterior ( > N_category_max ) is smaller than 
  TAIL_WEIGHT_CUTOFF, if it's normalizable, and that the 
  probability of N_category = N_category_max is smaller
  than TAIL_VALUE_CUTOFF, if not normalizable.
