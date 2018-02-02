include("simulate.jl")

data_dir = "/home/jbaker/documents/data/comparisons/"

x = readdlm( "$(data_dir)datasets/multivariate_lognormal/2" )

out_dir = "$(data_dir)/tuning/multivariate_lognormal/2/sghmc/"
mkpath( out_dir )

eta = 2.5e-4
alpha = 1.25

output = sim_sghmc( "multivariate_lognormal", x, eta, alpha )
writedlm( "$(out_dir)$(eta)-$(alpha)", output )
