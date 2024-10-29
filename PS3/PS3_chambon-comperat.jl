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
    for iy1 in 1:Ny
        for iy2 in 1:Ny
            # compute the consumption given y1, y2, a
            #
            #if ## check if the consumptions are negative
            #    # fill V
            #else
            #    # fill V
            #end
        end
    end
end

# Max and recover a_opt
for iy1 in 1:Ny
    for iy2 in 1:Ny
        # take the max over a with the findmax function and store the output
        # recover a_opt
        # recover c1_opt
        # recover c2_opt
    end
end


################ WITH CREDIT CONSTRAINT ##############################################################################

a_opt_bis = copy(a_opt);
for iy1 in 1:Ny
    for iy2 in 1:Ny
        #if ## check if the constraint binds 
        #    ## if yes, than ... + careful: -b 
        #end 
        # ## recover c_1
        # ## recover c_2
        # ## recover v 
    end
end


################ PLOTS ##############################################################################
using Plots


plot(grid_y,a_opt[:,40],title="Value Function", label="No credit constraint")
plot!(grid_y,a_opt_bis[:,40],title="Value Function", label="With credit constraint")


plot(grid_y,v_opt[:,40],title="Asset Policy Function", label="No credit constraint")
plot!(grid_y,v_opt_bis[:,40],title="Asset Policy Function", label="With credit constraint")

