__precompile__(false)
module OptProject

    export build_Q_c, newton, quasi_Newton

    using LinearAlgebra

    include("solvers.jl")
    include("helpers.jl")

end # module
