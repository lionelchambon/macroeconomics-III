#### Problem Set 2 - Markov chain process
using Distributions, LinearAlgebra


### Part A - Tauchen method ###############################################################

### Q1

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

### Q2

# Parameters value

rho = 0.8
sigma = 0.1225
num_states = 5;
N_iter = 2000;

## Applying the tauchen method (function below) to get the Markov chain process

 Pi, yy = tauchen(0, sigma, rho, num_states)  
 println(Pi) ## Generating the transition matrix

 println(collect(yy)) ## Generating the possible realizations of income
inc = exp.(yy) ## Recovering standardized values

### Q3

## Create a function to get the stable distribution

function get_invdist(first_guess, P)

    iter = 1000
    tol=1e-8
    π = (first_guess / sum(first_guess))'

    for iter in 1:iter
        π_new = π * P
        
        if norm(π_new - π) < tol
            return π_new / sum(π_new)  
        end
    
        π = π_new
    end

end


first_guess = [0, 0, 1, 0, 0]

invdist = get_invdist(first_guess, Pi)

println(invdist)

### Q4

## What is the mean income? 

expected_y = sum(invdist .* inc')
println("The expected income is $expected_y.")

## What is the share of households with income y_5?

share_y_5 = invdist[5]*100
println("The share of households with income y_5 is $share_y_5 %.")


### Part B - Recessions ###############################################################
P = [0.971 0.029 0.000 
    0.145 0.778 0.077
    0.000 0.508 0.492] 

### Q1

## Assume we are at normal growth, what is the proba of going in recession in t+1? 

# We can see straight away that if the probability of remaining in normal growth given that 
# we are in the state of normal growth is 0.971, the probability of being in any kind of  recession
# in t+1 is 1-0.971 = 0.029.

## In t+6? 

P6 = P^6 
println(P6)

# First, we raise the stochastic matrix to the power of 6. Then, by using the same computation as before,
# we find the probability of being in a recession in 6 months to be ~0.12%.

### Q2

first_guess_2 = [1, 0, 0]

invdist_2 = get_invdist(first_guess_2, P)
println(invdist_2)
println(invdist_2[3])

# The long-run probability of a severe recession is around 2.5%.











