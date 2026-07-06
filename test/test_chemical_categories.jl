
using ACEbase, Test
using AtomsBase: ChemicalSpecies, atomic_number   # activates ACEbaseAtomsBaseExt
using ACEbase: Categories, idx2cat, cat2idx, chemical_species, chemical_categories

##

@info("Testing ACEbase AtomsBase extension (chemical species)")

@testset "chemical_species converter" begin
   Si = ChemicalSpecies(:Si)
   @test chemical_species(14) == Si          # atomic number
   @test chemical_species(:Si) == Si         # symbol
   @test chemical_species("Si") == Si        # string
   @test chemical_species(Si) === Si         # identity / passthrough
   # unsupported types throw a clear error
   @test_throws ErrorException chemical_species(1.5)
   @test_throws ErrorException chemical_species([1, 2])
   @test_throws ErrorException chemical_species(nothing)
end

@testset "chemical_categories" begin
   # mixed input types (symbol / symbol / string), unsorted -> converted + sorted by Z
   c = chemical_categories([:Si, :O, "H"])
   @test c isa Categories{3, ChemicalSpecies}
   @test isbits(c)
   @test length(c) == 3
   # sorted by Z: H(1), O(8), Si(14)
   @test atomic_number.(idx2cat.(Ref(c), 1:3)) == [1, 8, 14]
   @test idx2cat(c, 1) == ChemicalSpecies(:H)
   # round-trip lookup (mixed query forms go through the stored ChemicalSpecies)
   @test cat2idx(c, ChemicalSpecies(:Si)) == 3
   @test cat2idx(c, ChemicalSpecies(:O))  == 2
   @test cat2idx(c, ChemicalSpecies(:C))  == 0     # absent -> 0
   @test ChemicalSpecies(:H) in c
   @test !(ChemicalSpecies(:C) in c)

   # duplicate species (":Si" and atomic number 14 are the same element) are rejected
   @test_throws ErrorException chemical_categories([:Si, 14])
end

@testset "many species (binary-search path)" begin
   # Z = 1:30 exercises the binary-search branch of cat2idx, which relies on
   # the `isless` ordering of ChemicalSpecies provided by AtomsBase
   c = chemical_categories(collect(1:30))
   @test c isa Categories{30, ChemicalSpecies}
   @test issorted(atomic_number.(c.list))
   for z = 1:30
      @test cat2idx(c, ChemicalSpecies(z)) == z
      @test idx2cat(c, z) == ChemicalSpecies(z)
   end
   @test cat2idx(c, ChemicalSpecies(31)) == 0
end
