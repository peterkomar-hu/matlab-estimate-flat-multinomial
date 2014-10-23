matlab-estimate-flat-multinomial
================================

This is a Matlab function that estimates the number of categories (N_category) of a flat multinomial distribution
  P(n_1, n_2, ... n_k) = [n! / (n_1! * n_2! * ... * n_k!)] * (1/k)^n
from a uniform prior and number of trials (N_trial) and number of different outcomes (N_different)

EstimateFlatMultinomial(N_trial, N_different, [N_category_max])
has two mandatory input arguments:
  N_trial (positive integer)
  N_different (positive integer, <= N_trial)
and one optional input argument:
  N_category_max (positive integer, >= N_different)

[N_category_dist, avg, variance, flag] = EstimateFlatMultinomial 
has four ouputs:
  N_category_dist =  N_category_max by 3 matrix:
    col 1: integers from 1 to N_category_max
    col 2: PDF of N_category
    col 3: (1 - CDF) of N_category
  avg = expectation value of N_category (=NaN if non-existent)
  variance = variance of N_category (=NaN if non-existent)
  flag = 0, posterior is normalized (on [1,Inf) )
       = 1, posterior is non-normalizable
