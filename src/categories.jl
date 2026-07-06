
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

Construct from a single iterable of distinct categories (order is irrelevant; the
entries are sorted on construction):
```julia
cats = Categories((14, 8))      # or Categories([8, 14])
@assert idx2cat(cats, 1) == 8   # sorted: O (Z=8) comes before Si (Z=14)
@assert cat2idx(cats, 14) == 2
```
"""
struct Categories{N, T}
   list::SVector{N, T}    # sorted, distinct categories;  index -> category (idx2cat)
end

# --- constructors ---------------------------------------------------------

function Categories(cats; allow_nonbits = false)
   v = collect(cats)
   T = eltype(v)
   (allow_nonbits || isbitstype(T)) ||
      error("Categories: the category type `$T` is not an `isbits` type, so the \
             resulting `Categories` could not be moved to a GPU or passed into a \
             kernel. Convert the categories to an `isbits` code first (e.g. integer \
             atomic numbers), or pass `allow_nonbits = true` to override.")
   sort!(v)
   for i = 2:length(v)
      v[i-1] == v[i] && error("Categories: categories must be distinct, got a \
                               repeated entry $(v[i])")
   end
   N = length(v)
   return Categories{N, T}(SVector{N, T}(v))
end


# --- forward / inverse maps -----------------------------------------------

"""
`idx2cat(cats, i)` : return the `i`-th category.
"""
idx2cat(cats::Categories, i::Integer) = cats.list[i]

"""
`cat2idx(cats, z) -> Int` : return the 1-based index of category `z`, or `0` if
`z` is not a category. GPU-safe: no allocation and no error/boxing on the
not-found path. Uses a linear scan for small `N` and a binary search for large
`N`. To look up a batch of categories, broadcast: `cat2idx.(Ref(cats), zs)`.
"""
@inline function cat2idx(cats::Categories{N}, z) where {N}
   list = cats.list
   if N <= _NLINEAR_CAT
      @inbounds for i = 1:N
         list[i] == z && return i
      end
   else
      k = searchsortedfirst(list, z)
      (k <= N && @inbounds(list[k] == z)) && return k
   end
   # if the search failed, return 0 (GPU-safe)
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

# treat a `Categories` as a scalar in broadcasting, so batched lookups can be
# written `cat2idx.(cats, zs)` (rather than `cat2idx.(Ref(cats), zs)`)
Base.broadcastable(cats::Categories) = Ref(cats)

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


# --- AtomsBase extension entry points --------------------------------------
# The methods for these live in `ext/ACEbaseAtomsBaseExt.jl` and are only
# available once `AtomsBase` is loaded (`using AtomsBase`).

"""
`chemical_species(x) -> ChemicalSpecies` : convert `x` to an AtomsBase
`ChemicalSpecies`, where `x` may be an `Integer` (atomic number), a `Symbol`, an
`AbstractString`, or a `ChemicalSpecies`. Throws for any other type.

Requires `AtomsBase` to be loaded (provided by the `ACEbaseAtomsBaseExt`
extension).
"""
function chemical_species end

"""
`chemical_categories(list) -> Categories{N, ChemicalSpecies}` : build a
[`Categories`](@ref) from a list of chemical species. Each entry of `list` is
passed through [`chemical_species`](@ref), so entries may be atomic numbers,
symbols, strings, or `ChemicalSpecies` (mixed). Duplicate species are rejected.

Requires `AtomsBase` to be loaded (provided by the `ACEbaseAtomsBaseExt`
extension).
"""
function chemical_categories end
