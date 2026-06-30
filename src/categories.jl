
# Categories : a small, sorted, isbits "list of categories" datastructure
# together with the bijection between a category and its 1-based index.
#
# This consolidates the `zlist` / `_i2z` / `_z2i` machinery that several ACEsuit
# packages (ACEpotentials, ACEradials, AtomicOrbitalKernels, ...) reimplement.
# It is fully generic over the category type `T` (chemical species are the case
# `T = Int` holding atomic numbers) and carries no chemistry-specific code: a
# consumer that needs e.g. `:Si -> 14` conversion performs that conversion
# *before* constructing a `Categories`.

using StaticArrays: SVector, SMatrix

# intentionally not exported - use explicit imports:
#   ACEbase.Categories, ACEbase.idx2cat, ACEbase.cat2idx, ACEbase.make_smatrix,
#   ACEbase.catcat2idx, ACEbase.catcat2idx_sym, ACEbase.symidx


# Below this number of categories `cat2idx` uses a linear scan; above it a
# binary search. For small `N` (the dominant MLIP case, 2-5 species) the linear
# scan is unrolled by the compiler and beats binary search.
const _NLINEAR_CAT = 8

# Above this number of categories `make_smatrix` returns a dense `Matrix`
# instead of an `SMatrix` to avoid pathological static-array sizes.
const _MAX_SMATRIX_NZ = 16


"""
`Categories{N, T}` : a **sorted** list of `N` distinct categories of type `T`,
together with the bijection between a category and its 1-based index. The sorted
position is the index used to address per-category data, so `cat2idx` is a binary
search (resp. a linear scan for small `N`) over the stored list. For chemical
species, `T = Int` holds atomic numbers. The struct is `isbits` whenever `T` is.

Construct with any iterable of distinct categories (order is irrelevant; the
entries are sorted on construction):
```julia
cats = Categories((14, 8))      # or Categories(14, 8), or Categories([8, 14])
@assert idx2cat(cats, 1) == 8   # sorted: O (Z=8) comes before Si (Z=14)
@assert cat2idx(cats, 14) == 2
```
"""
struct Categories{N, T}
   list::SVector{N, T}    # sorted, distinct categories;  index -> category (idx2cat)
end

# --- constructors ---------------------------------------------------------

function Categories(cats)
   v = collect(cats)
   sort!(v)
   for i = 2:length(v)
      v[i-1] == v[i] && error("Categories: categories must be distinct, got a \
                               repeated entry $(v[i])")
   end
   N = length(v)
   T = eltype(v)
   return Categories{N, T}(SVector{N, T}(v))
end

# varargs convenience: `Categories(8, 14)`, `Categories(:C, :H, :O)`
Categories(a, b, cs...) = Categories((a, b, cs...))

# --- forward / inverse maps -----------------------------------------------

"""
`idx2cat(cats, i)` : return the `i`-th category. Also defined for a raw list of
categories (tuple / vector), where it is just `getindex`.
"""
idx2cat(cats::Categories, i::Integer) = cats.list[i]
idx2cat(cats, i::Integer) = cats[i]

"""
`cat2idx(cats, z) -> Int` : return the 1-based index of category `z`, or `0` if
`z` is not a category. GPU-safe: no allocation and no error/boxing on the
not-found path. Uses a linear scan for small `N` and a binary search for large
`N`. Also defined for a raw list of categories (tuple / vector) via a linear
search, and for an `AbstractVector` of categories (returns a vector of indices,
mirroring a batched species-index lookup).
"""
@inline function cat2idx(cats::Categories{N}, z) where {N}
   list = cats.list
   if N <= _NLINEAR_CAT
      @inbounds for i = 1:N
         list[i] == z && return i
      end
      return 0
   else
      k = searchsortedfirst(list, z)
      (k <= N && @inbounds(list[k] == z)) && return k
      return 0
   end
end

# batched lookup (e.g. to feed a vector of indices into a kernel)
cat2idx(cats::Categories, zs::AbstractVector) = [ cat2idx(cats, z) for z in zs ]

# generic linear search over a raw list of categories
function cat2idx(cats, z)
   @inbounds for i in eachindex(cats)
      cats[i] == z && return i
   end
   return 0
end

# --- collection interface --------------------------------------------------

Base.length(::Categories{N}) where {N} = N
Base.eltype(::Categories{N, T}) where {N, T} = T
Base.getindex(cats::Categories, i::Integer) = cats.list[i]
Base.iterate(cats::Categories, state=1) =
      state > length(cats) ? nothing : (cats.list[state], state+1)
Base.in(z, cats::Categories) = (cat2idx(cats, z) != 0)
Base.:(==)(c1::Categories, c2::Categories) = (c1.list == c2.list)

Base.show(io::IO, cats::Categories) = print(io, "Categories(", Tuple(cats.list), ")")

# --- per-pair parameter matrices ------------------------------------------

"""
`make_smatrix(obj, NZ)` : build a square per-pair parameter container of size
`NZ x NZ` from `obj`:
- if `obj` is already an `NZ x NZ` matrix, it is converted (entries preserved);
- otherwise `obj` is treated as a scalar and broadcast to every entry.
For `NZ <= $(_MAX_SMATRIX_NZ)` an `SMatrix{NZ,NZ}` is returned; for larger `NZ`
a dense `Matrix` is returned instead, to avoid pathological static-array sizes.
"""
function make_smatrix(obj, NZ)
   if obj isa AbstractMatrix
      size(obj) == (NZ, NZ) ||
         error("make_smatrix: matrix input must have size ($NZ, $NZ)")
      return NZ <= _MAX_SMATRIX_NZ ? SMatrix{NZ, NZ}(obj) : Matrix(obj)
   end
   if obj isa AbstractArray
      error("make_smatrix: an array input must be an ($NZ, $NZ) matrix")
   end
   return NZ <= _MAX_SMATRIX_NZ ? SMatrix{NZ, NZ}(fill(obj, (NZ, NZ))) :
                                  fill(obj, (NZ, NZ))
end

# --- pair-index utilities --------------------------------------------------

"""
`symidx(j1, j2, n)` : index into a flattened symmetric `n x n` matrix for the
(unordered) pair `(j1, j2)`. E.g. for `n == 3` the layout is
```
   [ 1 2 3
     2 4 5
     3 5 6 ]
```
"""
function symidx(j1, j2, n)
   i1, i2 = minmax(j1, j2)
   return i2 + n * (i1 - 1) - (i1 * (i1 - 1)) ÷ 2
end

"""
`catcat2idx(cats1, cats2, a1, a2)` / `catcat2idx(cats, a1, a2)` : linear index of
the ordered category pair `(a1, a2)`, i.e. `(i1-1) * length(cats2) + i2` where
`i1 = cat2idx(cats1, a1)`, `i2 = cat2idx(cats2, a2)`.
"""
catcat2idx(cats1, cats2, a1, a2) =
      (cat2idx(cats1, a1) - 1) * length(cats2) + cat2idx(cats2, a2)

catcat2idx(cats, a1, a2) = catcat2idx(cats, cats, a1, a2)

"""
`catcat2idx_sym(cats, a1, a2)` : like `catcat2idx` but for symmetric pairs, i.e.
`(a1, a2)` and `(a2, a1)` map to the same index (see `symidx`).
"""
catcat2idx_sym(cats, a1, a2) =
      symidx(cat2idx(cats, a1), cat2idx(cats, a2), length(cats))
