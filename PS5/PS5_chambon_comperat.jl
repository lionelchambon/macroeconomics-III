##############################################################################################################################
## BLOCK 0 ##############################################################################################################################
##############################################################################################################################
using Distributions, LinearAlgebra

##Define Parameters
struct par_model
    beta::Float64              # discount factor
    w::Float64           # wage 
    mu::Float64            # risk aversion from CRRA parameter
    maxits
    tol
end

par = par_model(0.95,1.0,2.0,3000,1e-6)

##Defining the Grid for the Endogenous State Variable: Capital
Na = 300;
b = -0.2;
amax =  60;
agrid = collect(range(b, length = Na, stop = amax));

##Defining the Grid for the Exogenous State Variable: Technology Shock
rho = 0.9             # persistence of the AR(1) process
sigma = 0.2              # standard deviation of the AR(1) process
Ns = 5
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

prob = tauchen(0, sigma, rho, Ns)[1];
logs = tauchen(0, sigma, rho, Ns)[2];
sgrid = exp.(logs);



##############################################################################################################################
## BLOCK 3 ##############################################################################################################################
##############################################################################################################################

## Intermediate, do not touch, before the 25/11 session 
# we just loop over 3 possible interest rates
r_vec = [0.01, 0.02, 0.05]
Vguess = zeros(Ns,Na)
agg_a = zeros(length(r_vec))

for (ir,r) in enumerate(r_vec)

    agg_a[ir] = aiyagari(r, par, agrid, sgrid,prob, Vguess) 

end 



##############################################################################################################################
## BLOCK 2 ##############################################################################################################################
##############################################################################################################################

aiyagari = function(r, par, agrid, sgrid, prob, Vguess)
    Ns = length(sgrid)
    Na = length(agrid)

    # Call the VFI function and get the policy index 
    policy_a_index = VFI(r, agrid, sgrid, Vguess, prob, par)

    ##### 1. Building the transition matrix  ####################
    # Build Q as a 4D array
    Q = zeros(Ns,Na,Ns,Na)

    ### TO FILL ####### 
    # Compute Q 
    ##################

    # Then reshape it if Q was 4D
    Q = reshape(Q, Ns*Na, Ns*Na)

    # Check that the rows sum to 1! 
    sum(Q, dims=2)

    ###### 2. Computing the stable distribution  ###################################
    dist = ones(1, Ns*Na) / (Ns*Na);
    dist = get_stable_dist(dist, Q,par)

    # Check that the distribution vector sums to 1! 
    sum(dist)
    
    # Reshape dist as a Ns x Na dimension (more readable) (don't have to)
    # dist = reshape(dist, Ns, Na)

    ###### 3. Computing the aggregate #############################################
    agg_a = 0
    ### TO FILL ####### 
    # Compute agg_a 
    ##################

    return agg_a
end 




function get_stable_dist(invdist, P,par)
    for iter in 1:par.maxits
        invdist2 = invdist * P
        if maximum(abs, invdist2 .- invdist) < 1e-9
            println("Found solution after $iter iterations")
            return invdist2
        elseif iter == par.maxits
            error("No solution found after $iter iterations")
            return invdist
        end
        err = maximum(abs, invdist2 - invdist)
        invdist = invdist2
    end
end





##############################################################################################################################
## BLOCK 1 ##############################################################################################################################
##############################################################################################################################

VFI = function(r, agrid, sgrid, V0, prob, par)

    Ns = length(sgrid)
    Na = length(agrid)

    U = zeros(Ns,Na,Na)

    for is in 1:Ns                     # Loop Over skills Today
        for ia in 1:Na                 # Loop Over assets Today
            for ia_p in 1:Na           # Loop Over assets Tomorrow
                a = agrid[ia];     # Technology Today
                a_p = agrid[ia_p];     # Capital Today
                s = sgrid[is];    # Capital Tomorrow
                # Solve for Consumption at Each Point
                c = (1+r)*a + s*par.w - a_p
                if c .< 0
                    U[is,ia,ia_p] = -10^6;
                else()
                    U[is,ia,ia_p] = c^(1-par.mu)/(1-par.mu);
                end
            end
        end
    end


    Vnew = copy(V0);  # The new value function I obtain after an iteration
    Vguess = copy(V0);  # the  value function from which I start in each new iteration
    policy_a_index = Array{Int64,2}(undef,Ns,Na);
    tv = zeros(Na)

    ### TO FILL ####### 
    # VFI loop
    ##################



    return policy_a_index

end 