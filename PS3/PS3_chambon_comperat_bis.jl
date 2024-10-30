# Defining parameters
beta = 0.9
r = 0.2
Ny = 50        # Number of Grid Points
Na = 100        # Number of Grid Points
b = 1


################ INITIALIZATION ##############################################################################

##Defining the Grid for the Endogenous State Variable: Capital
ymin = 5; ymax = 10;    # Bounds for Grid
grid_y = collect(range(ymin, ymax, length = Ny)); # Grid for income

amin = -5; amax = 5;    # Bounds for Grid
grid_a = collect(range(amin, amax, length = Na)); # Grid for assets

# No credit constraint 
v_opt = zeros(Ny,Ny);
index_a_opt= Array{Int64,2}(undef,Ny,Ny); # Hint! 
c_1_opt = zeros(Ny,Ny);
c_2_opt = zeros(Ny,Ny);
a_opt = zeros(Ny,Ny);

# With credit constraint
v_opt_bis = zeros(Ny,Ny);
c_1_opt_bis = zeros(Ny,Ny);
c_2_opt_bis = zeros(Ny,Ny);
a_opt_bis = zeros(Ny,Ny);


################ NO CREDIT CONSTRAINT ##############################################################################

# Fill the V array for all possible choice of a 
V = zeros(Ny,Ny,Na);
for ia in 1:Na
    a = grid_a[ia]
    for iy1 in 1:Ny
        y1 = grid_y[iy1]
        for iy2 in 1:Ny
            y2 = grid_y[iy2]
            # compute the consumption given y1, y2, a
            c1 = y1 - a
            c2 = y2 + (1 + r)*a
            if c1 > 0 && c2 > 0
                V[iy1, iy2, ia] = log(c1) + beta*log(c2)
                c_1_opt[iy1, iy2] = c1
                c_2_opt[iy1, iy2] = c2
                a_opt[iy1, iy2] = a

            else
                V[iy1, iy2, ia] = -10^9
            end
        end
    end
end

# Max and recover a_opt
for iy1 in 1:Ny
    for iy2 in 1:Ny
# Take the max over a with the findmax function and store the output
        v_opt[iy1, iy2], index_a_opt[iy1, iy2] = findmax(V[iy1, iy2, :])
        
        a_opt[iy1, iy2] = grid_a[index_a_opt[iy1, iy2]]  

        c_1_opt[iy1, iy2] = grid_y[iy1] - a_opt[iy1, iy2]  
        c_2_opt[iy1, iy2] = grid_y[iy2] + (1 + r) * a_opt[iy1, iy2]  
    end
end


################ WITH CREDIT CONSTRAINT ##############################################################################

a_opt_bis = copy(a_opt)

binding = findall(x -> x < -b, a_opt_bis)
binding_set = Set(binding)

for iy1 in 1:Ny
    for iy2 in 1:Ny

        if (iy1, iy2) in binding_set
 
            a_opt_bis[iy1, iy2] = -b

            c_1_opt_bis[iy1, iy2] = grid_y[iy1] - a_opt_bis[iy1, iy2]  # Optimal c1 under constraint
            c_2_opt_bis[iy1, iy2] = grid_y[iy2] + (1 + r) * a_opt_bis[iy1, iy2]  # Optimal c2 under constraint

            if c_1_opt_bis[iy1, iy2] > 0 && c_2_opt_bis[iy1, iy2] > 0
                v_opt_bis[iy1, iy2] = log(c_1_opt_bis[iy1, iy2]) + beta * log(c_2_opt_bis[iy1, iy2])
            else
                v_opt_bis[iy1, iy2] = -10^9  
            end
        else
            # If constraint does not bind, retain the original values
            c_1_opt_bis[iy1, iy2] = c_1_opt[iy1, iy2]
            c_2_opt_bis[iy1, iy2] = c_2_opt[iy1, iy2]
            v_opt_bis[iy1, iy2] = v_opt[iy1, iy2]
        end
    end
end


################ PLOTS ##############################################################################
using Plots


plot(grid_y,a_opt[:,40],title="Value Function", label="No credit constraint")
plot!(grid_y,a_opt_bis[:,40],title="Value Function", label="With credit constraint")


plot(grid_y,v_opt[:,40],title="Asset Policy Function", label="No credit constraint")
plot!(grid_y,v_opt_bis[:,40],title="Asset Policy Function", label="With credit constraint")