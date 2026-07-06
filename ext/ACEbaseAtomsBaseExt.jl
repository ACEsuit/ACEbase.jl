module ACEbaseAtomsBaseExt

# Chemical-species conveniences for `ACEbase.Categories`, available when
# `AtomsBase` is loaded. Provides:
#   - `chemical_species` : ? -> ChemicalSpecies  (autoconverter)
#   - `chemical_categories` : list of species -> Categories{N, ChemicalSpecies}

using ACEbase: Categories
import ACEbase: chemical_species, chemical_categories
import AtomsBase: ChemicalSpecies


# -------------------- ? -> ChemicalSpecies autoconverter

chemical_species(s::ChemicalSpecies) = s
chemical_species(z::Integer) = ChemicalSpecies(z)          # z = atomic number
chemical_species(s::Symbol) = ChemicalSpecies(s)
# AtomsBase has no String constructor, so route strings through a Symbol
chemical_species(s::AbstractString) = ChemicalSpecies(Symbol(s))

chemical_species(x) = error(
      "chemical_species: cannot convert `$(repr(x))::$(typeof(x))` to a \
       `ChemicalSpecies`; expected an `Integer` (atomic number), `Symbol`, \
       `AbstractString`, or `ChemicalSpecies`.")


# -------------------- list of species -> Categories{N, ChemicalSpecies}

chemical_categories(list) = Categories([ chemical_species(x) for x in list ])


end
