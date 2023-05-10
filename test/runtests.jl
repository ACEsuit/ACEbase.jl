using ACEbase
using Test

# the main point of this package is having the interface functions 
# but we also have some FIO functionality and that we should test 
# rigorously to not break this important functionality. 

@testset "ACEbase.jl" begin
    include("test_fio.jl")
end
