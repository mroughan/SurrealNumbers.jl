# set up parameters for all simulations

Ns = [250 500 1000 2000 4000]
Lambdas = [0.5 1 1.5 2 2.5 3.0 3.5]
Alphas = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9]
Gs = [1 5 10]

# Ns = [2000]
# Lambdas = [2.0]
# Alphas = [0.6]
# Gs = [1]

# Ns = [5000]
# Lambdas = [1.5]
# Alphas = [0.2]
# Gs = [5]

m = 50 # number of iterations
Nseeds = 30
# Nseeds = 1
s0 = 0 # start seed number is (s0+1)

parent_dist = "P"
split_dist = "U"
# split_dist = "B"
distributions = "$(parent_dist)-$(split_dist)"

# parameters for determining convergence
m_final_states = 10 # number of final states to use to measure converged state
ϵ = 0.01 # how close to get to the convereged state
