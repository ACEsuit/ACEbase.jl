
using ACEbase, Test, StaticArrays
using ACEbase: Categories, idx2cat, cat2idx, make_smatrix,
               catcat2idx, catcat2idx_sym, symidx

##

@info("Testing ACEbase.Categories")

# brute-force reference for the inverse map
_ref_idx(list, z) = findfirst(==(z), list)

# round-trip + sorted invariant for a given (unsorted) input list
function test_roundtrip(input; allow_nonbits = false)
   cats = Categories(input; allow_nonbits)
   @test issorted(cats.list)
   @test length(cats) == length(input)
   # construction is order-independent
   @test cats == Categories(reverse(collect(input)); allow_nonbits)
   sorted = sort(collect(input))
   for i = 1:length(cats)
      @test idx2cat(cats, i) == sorted[i]
      @test cats[i] == sorted[i]
   end
   for z in input
      i = cat2idx(cats, z)
      @test idx2cat(cats, i) == z
      @test i == _ref_idx(cats.list, z)
      @test z in cats
   end
   return cats
end

##

@testset "construction & bijection" begin
   # few species (linear path)  and  integer categories
   c3 = test_roundtrip([14, 8, 1])
   @test idx2cat(c3, 1) == 1 && idx2cat(c3, 3) == 14   # sorted order
   @test cat2idx(c3, 8) == 2

   # non-integer but isbits categories (Char) work without the opt-out
   cc = test_roundtrip(['c', 'a', 'b'])
   @test idx2cat(cc, 1) == 'a'
   @test isbits(cc)

   # non-integer, non-isbits categories (Symbols) need `allow_nonbits`
   cs = test_roundtrip([:Si, :O, :H, :C]; allow_nonbits = true)
   @test idx2cat(cs, 1) == :C    # sorted alphabetically

   # construction is from a single iterable (tuple or vector)
   @test Categories((8, 14, 1)) == Categories([1, 8, 14])

   # many species (binary-search path)
   cN = test_roundtrip(collect(reverse(1:100)))
   @test cat2idx(cN, 1) == 1 && cat2idx(cN, 100) == 100
   @test cat2idx(cN, 57) == 57

   # duplicates rejected
   @test_throws ErrorException Categories([1, 2, 2, 3])
end

@testset "isbits requirement" begin
   # non-isbits category types are rejected by default ...
   @test !isbitstype(Symbol)
   @test_throws ErrorException Categories([:Si, :O])
   @test_throws ErrorException Categories(["a", "b"])         # String also non-bits
   # ... but can be opted into for CPU-only use
   cs = Categories([:Si, :O]; allow_nonbits = true)
   @test cs isa Categories{2, Symbol}
   @test !isbits(cs)
   @test cat2idx(cs, :Si) == 2
   # isbits category types are accepted (and produce an isbits Categories)
   @test isbits(Categories([14, 8]))            # Int
   @test isbits(Categories(['a', 'c', 'b']))    # Char
   @test isbits(Categories([(1,2), (0,9)]))     # Tuple{Int,Int}
end

@testset "absent categories" begin
   # cat2idx is non-throwing (GPU-safe): returns 0 when absent
   c = Categories([8, 14])
   @test cat2idx(c, 6) == 0
   @test !(6 in c)
   # large-N (binary-search) path
   cN = Categories(collect(1:100))
   @test cat2idx(cN, 0) == 0
   @test cat2idx(cN, 101) == 0
   @test cat2idx(cN, 200) == 0
end

@testset "isbits & allocation-free lookup" begin
   c = Categories([14, 8, 1])
   @test isbits(c)
   # warm up then check the small-N lookup does not allocate (and never throws)
   cat2idx(c, 8); cat2idx(c, 999)
   @test (@allocated cat2idx(c, 8)) == 0
   @test (@allocated cat2idx(c, 999)) == 0
end

@testset "pair-index utilities" begin
   c = Categories([1, 8, 14])
   # batched lookup via broadcast (Categories is treated as a scalar)
   @test cat2idx.(c, [14, 1, 8, 99]) == [3, 1, 2, 0]

   # catcat2idx against a brute-force reference (row-major linear index)
   n = length(c)
   for (k1, z1) in enumerate((1, 8, 14)), (k2, z2) in enumerate((1, 8, 14))
      @test catcat2idx(c, z1, z2) == (k1 - 1) * n + k2
   end

   # symmetric pair index: (a1,a2) and (a2,a1) collide, and the n(n+1)/2
   # distinct unordered pairs map onto 1:n(n+1)/2
   seen = Int[]
   for z1 in (1, 8, 14), z2 in (1, 8, 14)
      @test catcat2idx_sym(c, z1, z2) == catcat2idx_sym(c, z2, z1)
      push!(seen, catcat2idx_sym(c, z1, z2))
   end
   @test sort(unique(seen)) == collect(1:(n*(n+1))÷2)
   @test symidx(2, 3, 3) == symidx(3, 2, 3) == 5
end

@testset "make_smatrix" begin
   # scalar broadcast -> SMatrix for small NZ
   M = make_smatrix(2.0, 3)
   @test M isa SMatrix{3, 3}
   @test all(M .== 2.0)
   # matrix passthrough (entries preserved)
   A = reshape(collect(1:9), 3, 3)
   @test make_smatrix(A, 3) == A
   @test make_smatrix(A, 3) isa SMatrix{3, 3}
   # wrong-size matrix errors
   @test_throws ErrorException make_smatrix(reshape(collect(1:6), 2, 3), 3)
   # large NZ -> dense Matrix fallback (no SMatrix)
   Mbig = make_smatrix(1.0, 32)
   @test Mbig isa Matrix && size(Mbig) == (32, 32) && all(Mbig .== 1.0)
   Abig = ones(32, 32)
   @test make_smatrix(Abig, 32) isa Matrix
end
