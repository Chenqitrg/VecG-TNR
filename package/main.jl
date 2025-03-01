module VecG_TNR
using LinearAlgebra
using Revise
export Group, CyclicGroup, DihedralGroup, ProductGroup, IntegerGroup, GroupElement, elements, identity_element, inverse, multiply, group_tree, group_iter
export block_matrix_svd, block_matrix_qr, partial_trace, outer_product, concatenate_matrices_with_metadata
export Obj, Sector, Mor, zero_obj, dual_obj, random_mor, zero_mor, identity_mor, get_group, get_sector_size, VecG_permutedims, VecG_dag, max_abs, is_accend, is_descend, is_cyclic, to_perm, VecG_permutesectors
export VecG_svd, VecG_qr, to_sector_outin
export VecG_tensordot, VecG_partial_trace, VecG_factorize, VecG_square, VecG_outprod
export TubeSector, TubeMor

includet("groups.jl")
includet("block_matrix_calculation.jl")
includet("VecGtensor.jl")
includet("TubeG.jl")
includet("display.jl")
includet("VecGdecomposition.jl")
includet("VecGcontraction.jl")


end