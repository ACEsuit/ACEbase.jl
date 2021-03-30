module ACEbase

using Reexport



"""
`function fltype`

Return the output floating point type employed by some object, typically a
calculator or basis.
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


"""
`alloc_temp(args...)` : allocate temporary arrays for the evaluation of
some calculator or potential; see developer docs for more information
"""
alloc_temp(args...) = nothing

"""
`alloc_temp_d(args...)` : allocate temporary arrays for the evaluation of
some calculator or potential; see developer docs for more information
"""
alloc_temp_d(args...) = nothing

alloc_temp_dd(args...) = nothing



abstract type AbstractBasis end
abstract type ACEBasis <: AbstractBasis end
abstract type OneParticleBasis{T} <: ACEBasis end
abstract type ScalarACEBasis <: ACEBasis end


function evaluate end
function evaluate_d end
function evaluate_dd end
function evaluate_ed end
function evaluate! end
function evaluate_d! end
function evaluate_dd! end
function evaluate_ed! end
function precon! end


evaluate(basis::ACEBasis, args...) =
      evaluate!(alloc_B(basis), alloc_temp(basis), basis, args...)

evaluate_d(basis::ACEBasis, args...) =
      evaluate_d!(alloc_dB(basis), alloc_temp_d(basis), basis, args...)


# Some overloading to enable AD for some cases
#   TODO -> this needs to be extended..

evaluate(basis::ScalarACEBasis, args...) =
      evaluate!(alloc_B(basis, args...), alloc_temp(basis), basis, args...)

evaluate_d(basis::ScalarACEBasis, args...) =
      evaluate_d!(alloc_B(basis, args...), alloc_temp(basis), basis, args...)



function alloc_B end
function alloc_dB end


include("fio.jl")
@reexport using ACEbase.FIO


end
