include("postfunc.jl")

# @everywhere using MPI

@everywhere shotval = 450;
@everywhere maxt    = 218;

println("starting to do things")

@time include("postcalc.jl")

println("Done!")
