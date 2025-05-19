# Superselection sectors (quantum numbers):
# for defining graded vector spaces and invariant subspaces of tensor products
#==========================================================================================#

module TensorKitSectors

# exports
# -------
export FusionSector, BraidedSector, SymmetricSector, ModularSector
export AbelianGroup
export elements, order
export VecGωIrr, VecGIrr, TYIrr, QDAb, ModularDouble
export Sector, Group, AbstractIrrep
export Irrep
export take_center

export Nsymbol, Fsymbol, Rsymbol, Asymbol, Bsymbol
export sectorscalartype
export dim, sqrtdim, invsqrtdim, frobeniusschur, twist, fusiontensor, dual, sector_rank, S_matrix, T_vector
export otimes, deligneproduct, times, TimeReversed
export FusionStyle, UniqueFusion, MultipleFusion, SimpleFusion, GenericFusion,
       MultiplicityFreeFusion
export BraidingStyle, NoBraiding, SymmetricBraiding, Bosonic, Fermionic, Anyonic
export SectorSet, SectorValues, findindex
export rightone, leftone

export pentagon_equation, hexagon_equation

export Trivial, Z2Irrep, Z3Irrep, Z4Irrep, ZNIrrep, U1Irrep, SU2Irrep, CU1Irrep
export ProductSector
export FermionParity, FermionNumber, FermionSpin
export PlanarTrivial, FibonacciAnyon, IsingAnyon, Semion

export is_abelian, is_modular

# unicode exports
# ---------------
export ⊠, ⊗, ×
export ℤ, D, ℤ₂, ℤ₃, ℤ₄, D₃, D₄, D₅, U₁, SU, SU₂, CU₁, χ
export fℤ₂, fU₁, fSU₂
export ℨ
export QDℤ
export 𝔐at


# imports
# -------
using Base: SizeUnknown, HasLength, IsInfinite
using Base: HasEltype, EltypeUnknown
using Base.Iterators: product, filter
using Base: tuple_type_head, tuple_type_tail

using LinearAlgebra: tr
using TensorOperations
using HalfIntegers
using WignerSymbols

# includes
# --------
include("auxiliary.jl")
include("sectors.jl")
include("trivial.jl")
include("groups.jl")
include("irreps.jl")    # irreps of symmetry groups, with bosonic braiding
include("product.jl")   # direct product of different sectors
include("fermions.jl")  # irreps with defined fermionparity and fermionic braiding
include("anyons.jl")    # non-group sectors
include("group_related_cats/TambaraYamagami.jl")
include("group_related_cats/pointedfusion.jl")
include("group_related_cats/display.jl")
include("braidedsectors.jl")
include("drinfeldcenter.jl")
include("2Mat.jl")

# precompile
# ----------
include("precompile.jl")

function __precompile__()
    for I in (Trivial, Z2Irrep, Z3Irrep, Z4Irrep, ZNIrrep, U1Irrep, SU2Irrep, CU1Irrep,
              FermionParity, FermionNumber, FermionSpin, PlanarTrivial, FibonacciAnyon,
              IsingAnyon)
        precompile_sector(I)
    end
end

end # module TensorKitSectors
