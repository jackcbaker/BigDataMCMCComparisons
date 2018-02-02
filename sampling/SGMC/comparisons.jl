include("simulate.jl")
include("./models/models.jl")

"""
Make random runs from SGLD and SGHMC for a given model
"""
function comparison( model::statistical_model, dimension::Int64, seed::Int64,
                data_dir::AbstractString, step_const::Float64, eta::Float64, alpha::Float64 )
            
    # Simulate using SGLD
    mkpath( "$(data_dir)/sgld/$seed" )
    out_file = "$(data_dir)sgld/$seed/theta"
    srand( seed )
    output = sim_sgld( model, step_const )
    writedlm( out_file, output )

    # Simulate using SGHMC
    mkpath( "$(data_dir)sghmc/$seed" )
    out_file = "$(data_dir)sghmc/$seed/theta"
    output = sim_sghmc( model, eta, alpha )
    writedlm( out_file, output )
end

function comparison_sgld( model::lda, seed::Int64, data_dir::AbstractString, step_const::Float64 )
            
    # Simulate using SGLD
    mkpath( "$(data_dir)samples/lda/sgld/compare/$seed" )
    out_file = "$(data_dir)samples/lda/sgld/compare/$seed/theta"
    output = sim_sgld( model, 10, step_const )
    writedlm( out_file, output )
end

function comparison_sghmc( model::lda, seed::Int64, data_dir::AbstractString, eta::Float64, alpha::Float64 )
            
    # Simulate using SGLD
    mkpath( "$(data_dir)samples/lda/sghmc/compare/$seed" )
    out_file = "$(data_dir)samples/lda/sghmc/compare/$seed/theta"
    output = sim_sghmc( model, 10, eta, alpha )
    writedlm( out_file, output )
end
