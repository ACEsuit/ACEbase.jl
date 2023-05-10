
using Test, ACEbase, SparseArrays, StaticArrays
using ACEbase.Testing 

##

# list of object to test FIO for ... 

objects = [ 
   rand(100), 
   rand(10, 30), 
   rand(ComplexF64, 30, 30), 
   rand(1:10, 3, 3, 3), 
   # sparse
   sprand(100, 100, 0.01), 
   sprand(10, 1000, 0.003), 
   # static arrays 
   (@SVector rand(5)),
   (@MVector rand(Float32, 3)), 
   (@SMatrix [1 2 3; 4 5 6; 7 8 9]), 
   #  complex static arrays are not yet supported
]

for O in objects
   println_slim(@test all(test_fio(O)))
end


