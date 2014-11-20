function [N_category_dist, avg, variance, flag] = ...
    EstimateFlatMultinomial(N_trial, N_different, varargin)
% Estimates the number of categories of a uniform multinomial distribution
%
% ESTIMATEFLATMULTINOMIAL(N_trial, N_different, [N_category_max])
%   solves the Bayesian estimation problem of a uniform multinomial
%   distribution of unknown dimension (i.e. number of categories, N_category)
%   using the number of trials (N_trial) and the number of different
%   outcomes (N_different) as input. Assumes a flat prior for N_category.
%
%   The sampling distribution is 
%       P(n_1, n_2, ... n_k) = [n! / (n_1! * n_2! * ... * n_k!)] * (1/k)^n
%           where                       
%           k  = N_category
%           n  = sum(n_j) = N_trial
%           n_j = the number of outcomes from category j
%           sum( sgn(n_j) ) = N_different
%   
%   The posterior distribution is
%       P(N_category = k) = A * k^(-n) * k! / ((k - N_different)!)
%           where A is a normalization constant
%
%   [N_category_dist, avg, variance, flag] = ESTIMATEFLATMULTINOMIAL 
%       has four ouputs:
%   N_category_dist =  N_category_max by 3 matrix:
%       col 1: integers from 1 to N_category_max
%       col 2: PDF of N_category
%       col 3: (1 - CDF) of N_category
%   avg = expectation value of N_category (=NaN if non-existent)
%   variance = variance of N_category (=NaN if non-existent)
%   flag = 0, posterior is normalized (on [1,Inf) )
%        = 1, posterior is non-normalizable
%
%   N_category_max is the largest possible value for N_category.
%   If specified, the posterior is normalized over [1,N_category_max].
%
%   If N_category_max is unspecified, 
%   it is considered to be effectively infinite, and
%   the normalization is carreid out over the entire [0,Inf) domain.
%   Similarly, avg and variance are calculated over [0,Inf), if they exist.
%   The maximal dimension of N_category_dist is determined
%   by the criteria that the probability in the tail of the
%   posterior ( > N_category_max ) is smaller than 
%   TAIL_WEIGHT_CUTOFF, if it's normalizable, and that the 
%   probability of N_category = N_category_max is smaller
%   than TAIL_VALUE_CUTOFF, if not normalizable.

TAIL_VALUE_CUTOFF = 1.0e-5;
TAIL_WEIGHT_CUTOFF = 1.0e-2;

% input check
if mod(N_trial,1) || N_trial < 1
    error('N_trial must be a positive integer');
end
if mod(N_different,1) || N_different < 1
    error('N_different must be a positive integer');
end
if N_trial < N_different
    error('N_trials must be => N_different');
end

if nargin > 2
% CASE 1:  N_category_max is specfified
    option = 1;
    
    % dealing with optional argument
    N_category_max = cell2mat(varargin(1));

    % input check
    if mod(N_category_max,1) || N_category_max < 1
        error('N_category_max must be a positive integer');
    end
    if N_category_max < N_different
        error('N_category_max must be => N_different');
    end

else
% CASE 2:  N_category_max is NOT specfified
    option = 2;
    
    % start with a large enough value
    N_category_max = 2*N_different;
end    

satisfied = 0;
while ~satisfied
% This loop runs only once if N_category_max input argument is specified.
% But it runs many times, to find the smallest N_category_max
% that satisfies the CUTOFF conditions, if N_category_max argument is not
% specified.
    
    % container for considered N_category values
    N_category = transpose(linspace(1, N_category_max, N_category_max));

    % N_category values with non-zero posterior probability 
    % (N_category >= N_different)
    M = N_category(N_different :  N_category_max);

    % container for non-zero PDF values
    pdf = zeros(length(M),1);

    % posterior distribution (not normalized)
    for i = 1:length(M)
        pdf(i) = exp( ...
            (N_trial - N_different) * log(N_different) ... this constant is here to avoid underflow
            -N_trial * log(M(i)) ...
            + gammaln( M(i) + 1 ) ...
            - gammaln( M(i)-N_different + 1 ) ...
        );
    end

    % normalization   
    
    if N_trial > N_different + 1 || option == 1
    % when the posterior is normalizable
        normalizable = 1;
    else
        normalizable = 0;
    end
    
    % the explicit part ( < N_category_max )
    norm = sum(pdf);
    
    if option == 2 && normalizable
        % weight in the tail
        % evaluated approximately (works for large N_category_max)
        weight_tail =  tailmoment(0, N_trial, N_different, N_category_max);
        
        % added to the weight in the pdf
        norm =  norm + weight_tail;
        
        % normalized probability in the tail
        prob_tail = weight_tail / norm;
    end    
    
    pdf = pdf / norm;

    if option == 1 
    % when N_category_max is specified as an input
        satisfied = 1;
    elseif option == 2 && normalizable && prob_tail < TAIL_WEIGHT_CUTOFF
    % when N_category_max is not specified in the input 
    % AND the posterior is normalizable
    % AND the missing tail has small weight
        satisfied = 1;
    elseif option == 2 && ~normalizable && pdf(length(M)) < TAIL_VALUE_CUTOFF
    % when N_category_max is not specified in the input 
    % AND the posterior is NOT normalizable
    % AND the last bin has small probability
        satisfied = 1;
    else
        % increase the range of N_category by a factor of 2
        N_category_max = 2 * N_category_max;   
    end
end


% CDF calculated from PDF
% cdf is the complement cumulative distribution ( 1 - CDF )
cdf = zeros(length(pdf),1);
for i = linspace((length(pdf)-1), 1, (length(pdf)-1))
    cdf(i) = cdf(i+1) + pdf(i+1);
end


% append values for the impossible region [1, N_different-1]
pdf = padarray(pdf, N_different - 1, 0, 'pre');
cdf = padarray(cdf, N_different - 1, 1, 'pre');

% concatenate columns to get the distribution
N_category_dist = horzcat(N_category, pdf, cdf);


% flag indicates if the posterior is normalizable
flag = 1- normalizable;


% average
% explicit contribution (from < N_category_max)
avg = sum(pdf .* N_category);
if option == 2 
    if N_trial > N_different + 2
        % tail contribution (from > N_category_max)
        avg = avg + tailmoment(1, N_trial, N_different, N_category_max) / norm;
    else
        % avg = Inf
        avg = NaN;
    end
end

% variance
% explicit contribution (from < N_category_max)
variance = sum(pdf .* N_category .* N_category);
if option == 2 
    if N_trial > N_different + 3
        % tail contribution (from > N_category_max)
        variance = variance + tailmoment(2, N_trial, N_different, N_category_max) / norm;
    else
        % variance = Inf
        variance = NaN;
    end
end
variance = variance - avg^2;


end


function s = tailmoment(order, N_trial, N_different, N_category_max)
% approximates the contribution to the statistical moment given by order
% hiding in the tail [N_category_max + 1,Inf) for the unnormalized 
% posterior distribution
%
% it calculates 
% sum(N_different^(N_trial-N_differnt)/M^(N_trial-order) * M! / (M-N_different)!,  
%   over M = N_category_max + 1 : Inf )
% approximately by using the asymptotic series expansion around M = Inf.
    N = N_trial;
    K = N_different;
    Mmax = N_category_max;
    r = order;
    
    s = exp(  (N-K) * ( log(K) - log(Mmax) ) ... 
            + (1+r) * log(Mmax) - log(N-K-(1+r)) ...
            - (N-K -(1+r)) / (2*Mmax) ...
            );
end
