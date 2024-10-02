#### Problem Set 2 - Markov chain process
using Distributions, LinearAlgebra


### Part A - Tauchen method ###############################################################
# Input: tauchen function (do not change)
function tauchen(mean, sd, rho, num_states; q=3)

    uncond_sd = sd/sqrt(1-rho^2)
    y = range(-q*uncond_sd, stop = q*uncond_sd, length = num_states)
    d = y[2]-y[1]

    Pi = zeros(num_states,num_states)

    for row = 1:num_states
      # end points
          Pi[row,1] = cdf(Normal(),(y[1] - rho*y[row] + d/2)/sd)
          Pi[row,num_states] = 1 - cdf(Normal(), (y[num_states] - rho*y[row] - d/2)/sd)

      # middle columns
          for col = 2:num_states-1
              Pi[row, col] = (cdf(Normal(),(y[col] - rho*y[row] + d/2) / sd) -
                             cdf(Normal(),(y[col] - rho*y[row] - d/2) / sd))
          end
    end

  yy = y .+ mean # center process around its mean

  Pi = Pi./sum(Pi, dims = 2) # renormalize

  return Pi, yy
end 

# Parameters value
# rho = ;
# sigma = ;
# N = ;
# N_iter = 2000;

## Applying the tauchen method (function below) to get the Markov chain process



## Create a function to get the stable distribution
#function get_invdist(*guess*,*something important*)
#    for iter in 1:N_iter
#        *fixed point equation*
#        *error computation*
#        if *tolerance criteria*
#            return *output*
#        elseif iter == N_iter
#            error("No solution found after $iter iterations")
#        end
#        *guess update*
#    end
#end


## Find the stable distribution for the process generated above 
# *initial guess* 
# *output = get_invdist(*initial guess*,*something important*)


## What is the mean income? 
# *answer* 

## What is the share of households with income y_5?
# *answer* 





### Part B - Recessions ###############################################################
P = [0.971 0.029 0.000 
    0.145 0.778 0.077
    0.000 0.508 0.492] 


## Assume we are at normal growth, what is the proba of going in recession in t+1? 
# * initial distribution*
# * proba computation *

## In t+6? 
# * proba computation *

## Compute the stable distribution
# *initial guess* 
# *output = get_invdist(*initial guess*,*something important*)










