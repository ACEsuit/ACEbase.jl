module ACEbase

using Reexport


function alloc_B end
function alloc_dB end
function fltype end
function alloc_temp end
function alloc_temp_d end
function evaluate! end
function evaluate_d! end
function write_dict end
function read_dict end

abstract type AbstractBasis end
abstract type ACEBasis <: AbstractBasis end
abstract type OneParticleBasis{T} <: ACEBasis end
abstract type ScalarACEBasis <: ACEBasis end

alloc_temp(::ACEBasis, args...) = nothing
alloc_temp_d(::ACEBasis, args...) = nothing

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



include("fio.jl")
@reexport using ACEbase.FIO


end
