module ACEbase

using Reexport, StaticArrays

include("def.jl")



abstract type AbstractBasis end
abstract type ACEBasis <: AbstractBasis end
abstract type OneParticleBasis{T} <: ACEBasis end
abstract type Discrete1pBasis{T} <: OneParticleBasis{T} end 
abstract type ScalarACEBasis <: ACEBasis end

abstract type AbstractState end

abstract type AbstractConfiguration end

abstract type AbstractContinuousState <: AbstractState end

abstract type AbstractDiscreteState <: AbstractState end

isdiscrete(::AbstractContinuousState) = false
isdiscrete(::AbstractDiscreteState) = true


"""
`function fltype`

Return the output floating point type employed by some object, typically a
calculator or basis.

DEPRECATE THIS!!!
"""
function fltype end

fltype(T::DataType) = T

fltype_intersect(o1, o2) =
   fltype_intersect(fltype(o1), fltype(o2))

fltype_intersect(T1::DataType, T2::DataType) =
   typeof(one(T1) * one(T2))

"""
`function rfltype end`

Return the real floating point type employed by some object,
typically a calculator or basis, this is normally the same as fltype, but
it can be difference e.g. `rfltype = real ∘ flype
"""
rfltype(args...) = real(fltype(args...))


function gradtype end

function valtype end 

valtype(basis::ACEBasis, cfg::AbstractConfiguration) =
      valtype(basis, zero(eltype(cfg)))


# """
# `alloc_temp(args...)` : allocate temporary arrays for the evaluation of
# some calculator or potential; see developer docs for more information
# """
# alloc_temp(args...) = nothing

# """
# `alloc_temp_d(args...)` : allocate temporary arrays for the evaluation of
# some calculator or potential; see developer docs for more information
# """
# alloc_temp_d(args...) = nothing

# alloc_temp_dd(args...) = nothing


# alloc_temp_d(basis::ACEBasis, cfg::AbstractConfiguration) =
#             alloc_temp_d(basis, length(cfg))


# function alloc_B end

# alloc_B(basis::ACEBasis, x) = zeros(valtype(basis, x), length(basis))

# function alloc_dB end

# alloc_dB(B::ACEBasis, args...) = zeros(gradtype(B, args...), length(B))

# alloc_dB(B::ACEBasis, cfg::AbstractConfiguration) = 
#             zeros( gradtype(B, zero(eltype(cfg))), 
#                    (length(B), length(cfg) ) )


function combine end

function evaluate end
function evaluate_d end
function evaluate_dd end
function evaluate_ed end
function evaluate! end
function evaluate_d! end
function evaluate_dd! end
function evaluate_ed! end
function precon! end


# evaluate(basis::ACEBasis, args...) =
#       evaluate!(alloc_B(basis, args...), alloc_temp(basis, args...), basis, args...)

# evaluate_d(basis::ACEBasis, args...) =
#       evaluate_d!(alloc_dB(basis, args...),
#                   alloc_temp_d(basis, args...), basis, args...)

# function evaluate_ed(basis::ACEBasis, args...)
#    B = alloc_B(basis, args...)
#    dB = alloc_dB(basis, args...)
#    evaluate_ed!(B, dB, alloc_temp_d(basis, args...), basis, args...)
#    return B, dB
# end


# Some overloading to enable AD for some cases
#   TODO -> this needs to be extended..

# evaluate(basis::ScalarACEBasis, args...) =
#       evaluate!(alloc_B(basis, args...), alloc_temp(basis), basis, args...)

# evaluate_d(basis::ScalarACEBasis, args...) =
#       evaluate_d!(alloc_B(basis, args...), alloc_temp(basis), basis, args...)



"""
a simple utility function to check whether two objects are equal
"""
_allfieldsequal(x1, x2) =
      all( getfield(x1, sym) == getfield(x2, sym)
           for sym in union(fieldnames(typeof(x1)), fieldnames(typeof(x2))) )


include("fio.jl")
@reexport using ACEbase.FIO

include("testing.jl")


# ---- some utility code for discrete basis sets to avoid any attempt to 
#      differentiate it

# TODO: blech - nasty hack...
alloc_dB(basis::Discrete1pBasis{T}, N::Integer) where {T} = 
            zeros(SVector{3, T}, length(basis), N)





end
