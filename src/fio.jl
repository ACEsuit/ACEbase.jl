
"""
`ACEbase.FIO` : provides some basis file IO. All ACEsuit 
packages should use these interface functions so that the file formats
can be changed later. """
module FIO

using JSON, YAML
using ZipArchives: ZipReader, ZipWriter, zip_newfile, zip_nentries, zip_name, zip_readentry
using SparseArrays: SparseMatrixCSC
using StaticArrays


export read_dict, write_dict,
       zip_dict, unzip_dict,
       load_json, save_json,
       load_dict, save_dict,
       load_yaml, save_yaml


#######################################################################
#                     Conversions to and from Dict
#######################################################################

# eventually we want to be able to serialise all ACE types to Dict
# and back and those Dicts may only contain elementary data types
# this will then allow us to load them back via read_dict without
# having to know the type in the code

"""
`read_dict(D::Dict) -> ?`

Looks for a key `__id__` in `D` whose value is used to dynamically dispatch
the decoding to
```julia
read_dict(Val(Symbol(D["__id__"])), D)
```
That is, a user defined type must implement this `convert(::Val, ::Dict)`
utility function. The typical code would be
```julia
module MyModule
    struct MyStructA
        a::Float64
    end
    Dict(A::MyStructA) = Dict( "__id__" -> "MyModule_MyStructA",
                               "a" -> a )
    MyStructA(D::Dict) = A(D["a"])
    Base.convert(::Val{:MyModule_MyStructA})
end
```
The user is responsible for choosing an id that is sufficiently unique.

The purpose of this function is to enable the loading of more complex JSON files
and automatically converting the parsed Dict into the appropriate composite
types. It is intentional that a significant burden is put on the user code
here. This will maximise flexibiliy, e.g., by introducing version numbers,
and being able to read old version in the future.
"""
function read_dict(D::Dict)
    if !haskey(D, "__id__")
        error("ACEbase.FIO.read_dict: `D` has no key `__id__`")
    end
    if haskey(D, "__v__")
      return read_dict(Val(Symbol(D["__id__"])),
                       Val(Symbol(D["__v__"])),
                       D)
    end
    return read_dict(Val(Symbol(D["__id__"])), D)
end

read_dict(::Val{sym}, D::Dict) where {sym} =
    error("ACEbase.FIO.read_dict no implemented for symbol $(sym)")


#######################################################################
#                     JSON & YAML 
#######################################################################


function load_json(fname::AbstractString)
    return JSON.parsefile(fname)
end

function save_json(fname::AbstractString, D::Dict; indent=0)
    f = open(fname, "w")
    JSON.print(f, D, indent)
    close(f)
    return nothing
end

function load_yaml(fname::AbstractString)
   return YAML.load_file(fname)
end

function save_yaml(fname::AbstractString, D::Dict)
   YAML.write_file(fname, D) 
   return nothing
end

function load_dict(fname::AbstractString)
   if endswith(fname, ".json")
      return load_json(fname)
   elseif endswith(fname, ".yaml") || endswith(fname, ".yml")
      return load_yaml(fname)
   else
      @warn("Unrecognised file format. Expected: \"*.json\" or \"*.yaml\", got filename: $(fname); default to json format")
      return load_json(fname)# throw(error("Unrecognised file format. Expected: \"*.json\" or \"*.yaml\", got filename: $(fname)"))
      # throw(error("Unrecognised file format. Expected: \"*.json\" or \"*.yaml\", got filename: $(fname)"))
   end
end

function save_dict(fname::AbstractString, D::Dict; indent=0)
   if endswith(fname, ".json")
      return save_json(fname, D; indent=indent)
   elseif endswith(fname, ".yaml") || endswith(fname, ".yml")
      return save_yaml(fname, D)
   else
      @warn("Unrecognised file format. Expected: \"*.json\" or \"*.yaml\", got filename: $(fname); default to json format")
      # throw(error("Unrecognised file format. Expected: \"*.json\" or \"*.yaml\", got filename: $(fname)"))
      return save_json(fname, D; indent=indent)
   end
end


function zip_dict(fname::AbstractString, D::Dict; indent=0)
   ZipWriter(fname) do zipdir
      zip_newfile(zipdir, "dict.json"; compress=true)
      write(zipdir, JSON.json(D, indent))
   end
   return nothing
end

function unzip_dict(fname::AbstractString; verbose=false)
   zipdir = ZipReader(read(fname))
   if zip_nentries(zipdir) != 1
      error("Expected exactly one file in `$(fname)`")
   end
   verbose && @show zip_name(zipdir, 1)
   return JSON.parse(zip_readentry(zipdir, 1, String))
end

#######################################################################
#                     FIO for several standard objects  
#######################################################################

# Datatype
write_dict(T::Type) = Dict("__id__" => "Type", "T" => string(T))
read_dict(::Val{:Type}, D::Dict) = Main.eval(Meta.parse(D["T"]))

# Complex Vector

function write_dict(A::AbstractArray{T}) where {T <: Number}
   D = Dict("__id__" => "ACE_ArrayOfNumber",
                 "T" => write_dict(T),
              "size" => size(A), 
              "real" => real.(A[:]))
   if T <: Complex
      D["imag"] = imag.(A[:])
   end
   return D
end

function read_dict(::Val{:ACE_ArrayOfNumber}, D::Dict)
   T = read_dict(D["T"])
   sz = tuple(D["size"]...)
   data = T.(D["real"])
   if T <: Complex
      data[:] .+= im .* D["imag"]
   end
   return collect(reshape(data, sz))
end

# General Array 

write_dict(A::AbstractArray{T}) where {T} = 
         Dict("__id__" => "ACE_Array",
                "size" => size(A),
                "vals" => write_dict.(A[:]) )

read_dict(::Val{:ACE_Array}, D::Dict) = 
         collect(reshape(read_dict.(D["vals"]), tuple(D["size"]...)))


# SparseMatrixCSC

write_dict(A::SparseMatrixCSC{TF, TI}) where {TF, TI} =
   Dict("__id__" => "SparseMatrixCSC",
        "TI" => write_dict(TI),
        "colptr" => A.colptr,
        "rowval" => A.rowval,
        "nzval" => write_dict(A.nzval),
        "m" => A.m,
        "n" => A.n )

function read_dict(::Val{:SparseMatrixCSC}, D::Dict)
   TI = read_dict(D["TI"])
   return SparseMatrixCSC(Int(D["m"]), Int(D["n"]),
                           TI.(D["colptr"]), TI.(D["rowval"]), 
                           read_dict(D["nzval"]))
end


# ------------ Static Arrays 
# start with the standard simples ones and think later about how to generalize 
# this sensibly ... 

function write_dict(A::StaticArray{S, T, N}) where {S, T <: Real, N}
   D = Dict("__id__" => "ACE_StaticArray",
                 "T" => write_dict(typeof(A)),
              "data" => A.data, )
   D["T"]["T"] = "StaticArrays." * D["T"]["T"]
   return D
end

function read_dict(::Val{:ACE_StaticArray}, D::Dict)
   T = read_dict(D["T"])
   return T(tuple(D["data"]...))
end

end
