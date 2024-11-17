using Distributions, LinearAlgebra

##Define Parameters
beta = 0.95 # Discount factor
rho = 0.9; # persistence of the AR(1) process
sigma = 0.2; # Standart deviation of the AR(1) process
Ns = 5; # Number of states of the Markov Chain discretization 
b = -0.2 # borrowing limit
r = 0.025 # interest rate 
w = 1  # wage 
mu = 2 # risk aversion. CRRA             

##Defining the Grid for the Endogenous State Variable: Capital
Na = 300;
amax =  60;
agrid = collect(range(b, length = Na, stop = amax));


##Defining the Grid for the Exogenous State Variable: Technology Shock
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



## Build the 3-Dimensional Contemporaneous Utility Grid for the System
# Initialize the 3D Array
U = zeros(Ns,Na,Na)

###### TO FILL ############################################################################################################
 for is in 1:Ns                     # Loop Over skills Today
     for ia in 1:Na                 # Loop Over assets Today
         for ia_p in 1:Na           # Loop Over assets Tomorrow

             c = sgrid[is] * w + (1 +r) * agrid[ia] - agrid[ia_p] 

             if c .< 0 
                 
                  U[is, ia, ia_p] = -1e10

             else()
                 
                  U[is, ia, ia_p] = (c^(1 - mu) - 1) / (1 - mu)

             end
         end
     end
 end
##################################################################################################################


##Value Function Iteration

#Initial Guess of the Value Function
V0 = zeros(Ns,Na);

#Second upgraded initial guess under a = a'. 
###### TO FILL ############################################################################################################

#Using what we discussed in class, but this increased my number of iterations to 136 instead of reducing it...

V0_alt = zeros(Ns, Na)
for is in 1:Ns
    for ia in 1:Na
        V0_alt[is, ia] = (1/(1-beta))*U[is, ia, ia]
    end
end

#Trying another guess:

V0_bis = zeros(Ns, Na)
for is in 1:Ns
    for ia in 1:Na
        c = sgrid[is] * w + (1 + r) * agrid[ia] - agrid[ia]  
        if c > 0
            V0_bis[is, ia] = (c^(1 - mu) - 1) / (1 - mu)  
        else
            V0_bis[is, ia] = -1e10  
        end
    end
end
Vguess = copy(V0_bis)

##################################################################################################################
#Calculate the Guess of the Expected Value Function

tol = 1e-4;
its = 0;
maxits = 3000; # Define the maximum number of iterations
Vnew = copy(Vguess);  # The new value function I obtain after an iteration
Vguess = copy(V0_alt);  # the  value function from which I start in each new iteration. Replace V0 with V0_alt or V0_bis to compare
policy_a_index = Array{Int64,2}(undef,Ns,Na);
tv = zeros(Na)

###### TO FILL ############################################################################################################
 for iter in 1:maxits
     for is in 1:Ns                     # Loop Over skills Today
         for ia in 1:Na                 # Loop Over assets Today
          for ia_p in 1:Na

          EV = sum(prob[is, :] .* Vguess[:, ia_p])
          tv[ia_p] = U[is, ia, ia_p] + beta * EV

          end

         Vnew[is, ia], policy_a_index[is, ia] = findmax(tv)

     end
    end

     if maximum(abs,Vguess.-Vnew) < tol
         println("Found solution after $iter iterations")
     return nothing
     elseif iter==maxits
         println("No solution found after $iter iterations")
     return nothing
     end
     err = maximum(abs,Vguess.-Vnew)
     global Vguess = copy(Vnew)
     println("#iter = $iter, εᵥ = $err")
 end
##################################################################################################################

# Found solution after 133 iterations using the initial guess.
# Found solution after 136 iterations using improved guess a=a'.
# Found solution after 131 iterations when scrapping the division by 1-beta.


##Policy function for assets
policy_a = Array{Float64,2}(undef,Ns,Na);

for is in 1:Ns
policy_a[is,:] = agrid[policy_a_index[is,:]] 
end


##Policy function for consumption
policy_c = Array{Float64,2}(undef,Ns,Na) 
###### TO FILL ############################################################################################################

for is in 1:Ns
  for ia in 1:Na
      policy_c[is, ia] = sgrid[is] * w + (1 + r) * agrid[ia] - policy_a[is, ia]
  end
end

##################################################################################################################


##Plotting the results
using Plots
vplot = plot(agrid, Vnew[:,:]', label="Value Function", xlabel="Assets (a)", ylabel="Value (V)", title="Value Function V(a, s)")
display(vplot)

aplot = plot(agrid, policy_a[:,:]', label="Policy for Assets", xlabel="Assets (a)", ylabel="Assets Tomorrow (a')", title="Policy Function for Assets")
display(aplot)

cplot = plot(agrid, policy_c[:,:]', label="Policy for Consumption", xlabel="Assets (a)", ylabel="Consumption (c)", title="Policy Function for Consumption")
display(cplot)
