module VecG_TNR
using LinearAlgebra
export Group, CyclicGroup, DihedralGroup, GroupElement, elements, identity_element, inverse, multiply, group_tree, group_iter
export block_matrix_svd, block_matrix_qr, partial_trace
export Obj, Sector, Mor, zero_obj, dual_obj, random_mor, zero_mor, identity_mor, get_group, get_sector_size, is_accend, is_cyclic, is_decend, to_perm, VecG_permutedims, VecG_dag, max_abs
export VecG_svd, VecG_qr
export VecG_tensordot, VecG_partial_trace, VecG_factorize

include("groups.jl")
include("block_matrix_calculation.jl")
include("VecGtensor.jl")
include("display.jl")
include("VecGdecomposition.jl")
include("VecGcontraction.jl")

end