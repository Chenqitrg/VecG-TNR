using Test

# Mock data
D4 = DihedralGroup(4)  # 假设这是一个群对象
e = GroupElement((0, 0), D4)
s = GroupElement((1, 0), D4)
r = GroupElement((0, 1), D4)
r2 = GroupElement((0, 2), D4)
sr2 = s * r2

# Create Objects
A = Object(Dict(e => 2, s*r => 3, s => 1, r2 => 2, sr2 => 1, r => 5), D4)

# Create BlockTensor
T = BlockTensor(D4, [A, A, A, A])

# Test set_block!
@testset "set_block! tests" begin
    valid_sector = (e, s*r, s, r)
    valid_value = rand(2, 3, 1, 5)
    set_block!(T, valid_sector, valid_value)
    @test get_block(T, valid_sector) == valid_value

    invalid_sector = (e, s, s, r)  # Does not multiply to identity
    invalid_value = rand(2, 1, 1, 5)
    @test_throws ArgumentError set_block!(T, invalid_sector, valid_value)
    @test_throws ArgumentError set_block!(T, valid_sector, invalid_value)
end

# Test addition
@testset "BlockTensor addition" begin
    T2 = BlockTensor(D4, [A, A, A, A])
    set_block!(T2, valid_sector, valid_value)
    T_sum = T + T2
    @test get_block(T_sum, valid_sector) == 2 * valid_value
end

# Test display
println("\nDisplaying BlockTensor T:")
show(stdout, MIME"text/plain", T)