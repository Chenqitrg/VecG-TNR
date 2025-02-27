module VecG_TNR
using LinearAlgebra
# using Revise
export Group, CyclicGroup, DihedralGroup, GroupElement, elements, identity_element, inverse, multiply, group_tree, group_iter
export block_matrix_svd, block_matrix_qr, partial_trace, outer_product
export Obj, Sector, Mor, zero_obj, dual_obj, random_mor, zero_mor, identity_mor, get_group, get_sector_size, VecG_permutedims, VecG_dag, max_abs, is_accend, is_descend
export VecG_svd, VecG_qr
export VecG_tensordot, VecG_partial_trace, VecG_factorize, VecG_square, VecG_outprod
export TubeSector, TubeMor

include("groups.jl")
include("block_matrix_calculation.jl")
include("VecGtensor.jl")
include("TubeG.jl")
include("display.jl")
include("VecGdecomposition.jl")
include("VecGcontraction.jl")


end