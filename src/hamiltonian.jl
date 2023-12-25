# const AllowedCoefficient = Union{Number, Function, Coefficient}

"""
	hamiltonian(lattice::FermionicLattice, t::Array{<:Real, 2}, v::Array{<:Real, 4}; atol::Real=1.0e-14)

Return the quantum chemistry Hamiltonian with the correct symmetry set in lattice, 
given the twobody coefficients t and fourbody coefficients v.

Definiton of the Hamiltonian can be found in docs/Hamiltonian.jpeg, for a molecule 
with L spatial orbital, t should have size (L, L) and v should have size (L,L,L,L).

Three different symmetries are supported for the Fermionic lattice:
* NoSymmetry, no symmetry
* ChargeCharge, U₁ ⊗ U₁ symmetry for up and down charge conservation
* SpinCharge, U₁ ⊗ SU₂ spin SU₂ symmetry and total U₁ chagre conservation
"""
function hamiltonian(t::Array{<:Real, 2}, v::Array{<:Real, 4}, symmetry::FermionicSymmetry; atol::Real=1.0e-14)
	@assert size(t, 1) == size(t, 2) == size(v, 1) == size(v, 2) == size(v, 3) == size(v, 4)
	L = size(t, 1)
	terms = []
	for pos1 in 1:L
		for pos2 in 1:L
			tmp = t[pos1, pos2]
			if abs(tmp) > atol
				m = twobody(pos1, pos2, tmp, symmetry)
				push!(terms, m)
			end
		end
	end
	for p in 1:L
		for q in 1:L
			for r in 1:L
				for s in 1:L
					tmp = v[p, q, r, s]
					if abs(tmp) > atol
						m = fourbody(p, q, r, s, tmp, symmetry)
						push!(terms, m)
					end
				end
			end
		end
	end
	physpace = physical_space(symmetry)
	return QuantumOperator([physpace for i in 1:L], terms)
end 



# """
# 	twobody(lattice::FermionicLattice, i::Int, j::Int, t::Number)

# Two-body term tĉᵢ†ĉⱼ	
# """
# function twobody(pos1::Int, pos2::Int, c::AllowedCoefficient, symmetry::FermionicSymmetry)
# 	min_pos = min(pos1, pos2) - 1
# 	t = twobody_util(pos1 - min_pos, pos2 - min_pos, c, symmetry)
# 	return Hamiltonians.shift(t, min_pos)
# end

# function twobody_util(pos1::Int, pos2::Int, c::AllowedCoefficient, symmetry::FermionicSymmetry)
# 	p = _select_sp(symmetry)
# 	adag = a_dagger(pos1, p)
# 	a = a_dagger(pos2, p)'
# 	r = adag * a
# 	return simplify(r) * Coefficient(c)
# end

# """
# 	fourbody(lattice::FermionicLattice, i::Int, j::Int, k::Int, l::Int, t::Number)

# Four-body term tĉᵢ†ĉⱼ†ĉₖĉₗ
# """
# function fourbody(pos1::Int, pos2::Int, pos3::Int, pos4::Int, v::AllowedCoefficient, symmetry::FermionicSymmetry)
# 	min_pos = min(pos1, pos2, pos3, pos4) - 1
# 	t = fourbody_util(pos1 - min_pos, pos2 - min_pos, pos3 - min_pos, pos4 - min_pos, v, symmetry)
# 	return Hamiltonians.shift(t, min_pos)
# end

# function fourbody_util(pos1::Int, pos2::Int, pos3::Int, pos4::Int, v::AllowedCoefficient, symmetry::FermionicSymmetry)
# 	p = _select_sp(symmetry)
# 	adag1 = a_dagger(pos1, p)
# 	adag2 = a_dagger(pos2, p)
# 	a3 = a_dagger(pos3, p)'
# 	a4 = a_dagger(pos4, p)'

# 	r = (adag1 * adag2) * (a3 * a4)
# 	return simplify(r) * Coefficient(v)
# end

# """
# 	creation(pos::Int, symmetry::FermionicSymmetry)

# Fermionic spinal creationa operator with SpinCharge symmetry.

# pos: the position of the creational operator
# """
# creation(pos::Int, symmetry::FermionicSymmetry) = a_dagger(pos, _select_sp(symmetry))

# a_dagger(pos::Int, p) = a_dagger_util(pos, p["+"], p["JW"])
# function a_dagger_util(pos::Int, adag, JW)
# 	_pos = collect(1:pos)
# 	_ops = Vector{Any}(undef, length(_pos)) 	
# 	for i in 1:length(_pos)-1
# 		_ops[i] = JW
# 	end
# 	_ops[end] = adag
# 	return QTerm(_pos, [_ops...])

# end

# """
# 	The convention is that the creation operator on the left of the annihilation operator

# By convention space_l of all the operators are vacuum
# """
# function spinal_fermion_site_ops_u1_su2()
# 	ph = Rep[U₁×SU₂]((0, 0)=>1, (2, 0)=>1, (1, 0.5)=>1)
# 	bh = Rep[U₁×SU₂]((1, 0.5)=>1)
# 	vh = oneunit(ph)
# 	# JW operator is contained in adag
# 	adag = TensorMap(zeros, Float64, vh ⊗ ph ← bh ⊗ ph)
# 	blocks(adag)[Irrep[U₁](1) ⊠ Irrep[SU₂](0.5)] = ones(1,1)
# 	blocks(adag)[Irrep[U₁](2) ⊠ Irrep[SU₂](0)] = sqrt(2) * ones(1,1) 

# 	bh = Rep[U₁×SU₂]((-1, 0.5)=>1)
# 	a = TensorMap(zeros, Float64, vh ⊗ ph ← bh ⊗ ph)
# 	blocks(a)[Irrep[U₁](1) ⊠ Irrep[SU₂](0.5)] = ones(1,1)
# 	blocks(a)[Irrep[U₁](0) ⊠ Irrep[SU₂](0)] = -sqrt(2) * ones(1,1) 


# 	onsite_interact = TensorMap(zeros, Float64, ph ← ph)
# 	blocks(onsite_interact)[Irrep[U₁](2) ⊠ Irrep[SU₂](0)] = ones(1, 1)

# 	JW = TensorMap(ones, Float64, ph ← ph)
# 	blocks(JW)[Irrep[U₁](1) ⊠ Irrep[SU₂](0.5)] = -ones(1, 1)

# 	# adagJW = TensorMap(zeros, Float64, vh ⊗ ph ← bh ⊗ ph)
# 	# blocks(adagJW)[Irrep[U₁](0) ⊠ Irrep[SU₂](0.5)] = ones(1,1)
# 	# blocks(adagJW)[Irrep[U₁](0.5) ⊠ Irrep[SU₂](0)] = -sqrt(2) * ones(1,1) 

# 	# hund operators
# 	# c↑† ⊗ c↓†
# 	bhr = Rep[U₁×SU₂]((2, 0)=>1)
# 	adagadag = TensorMap(ones, Float64, vh ⊗ ph ← bhr ⊗ ph)

# 	# c↑† ⊗ c↓, this is a spin 1 sector operator!!!
# 	bhr = Rep[U₁×SU₂]((1, 1)=>1)
# 	adaga = TensorMap(zeros, Float64, vh ⊗ ph ← bhr ⊗ ph)
# 	blocks(adaga)[Irrep[U₁](1) ⊠ Irrep[SU₂](0.5)] = ones(1, 1) * (-sqrt(3) / 2)

# 	n = TensorMap(ones, Float64, ph ← ph)
# 	blocks(n)[Irrep[U₁](0) ⊠ Irrep[SU₂](0)] = zeros(1, 1)
# 	blocks(n)[Irrep[U₁](2) ⊠ Irrep[SU₂](0)] = 2 * ones(1, 1)

# 	return Dict("+"=>adag, "-"=>a, "++"=>adagadag, "+-"=>adaga, "n↑n↓"=>onsite_interact, "JW"=>JW, "n"=>n)
# end

# function spinal_fermion_site_ops_u1_u1()
# 	ph = Rep[U₁×U₁]((1, 1)=>1, (1,0)=>1, (0,1)=>1, (0,0)=>1)
# 	vacuum = oneunit(ph)

# 	# adag
# 	adagup = TensorMap(zeros, Float64, vacuum ⊗ ph ← Rep[U₁×U₁]((1,0)=>1) ⊗ ph )
# 	blocks(adagup)[Irrep[U₁](1) ⊠ Irrep[U₁](0)] = ones(1,1)
# 	blocks(adagup)[Irrep[U₁](1) ⊠ Irrep[U₁](1)] = ones(1,1)

# 	# JW operator is not contained in adagdown
# 	adagdown = TensorMap(zeros, Float64, vacuum ⊗ ph ← Rep[U₁×U₁]((0,1)=>1) ⊗ ph)
# 	blocks(adagdown)[Irrep[U₁](0) ⊠ Irrep[U₁](1)] = ones(1,1)
# 	blocks(adagdown)[Irrep[U₁](1) ⊠ Irrep[U₁](1)] = -ones(1,1)

# 	adag = cat(adagup, adagdown, dims=3)

# 	# a
# 	aup = TensorMap(zeros, Float64, vacuum ⊗ ph ← Rep[U₁×U₁]((-1,0)=>1) ⊗ ph)
# 	blocks(aup)[Irrep[U₁](0) ⊠ Irrep[U₁](0)] = ones(1,1)
# 	blocks(aup)[Irrep[U₁](0) ⊠ Irrep[U₁](1)] = ones(1,1)

# 	adown = TensorMap(zeros, Float64, vacuum ⊗ ph ← Rep[U₁×U₁]((0,-1)=>1) ⊗ ph)
# 	blocks(adown)[Irrep[U₁](0) ⊠ Irrep[U₁](0)] = ones(1,1)
# 	blocks(adown)[Irrep[U₁](1) ⊠ Irrep[U₁](0)] = -ones(1,1)

# 	a = cat(aup, - adown, dims=3)

# 	# hund operators
# 	adagadag = TensorMap(zeros, Float64, vacuum ⊗ ph ← Rep[U₁×U₁]((1,1)=>1) ⊗ ph)
# 	blocks(adagadag)[Irrep[U₁](1) ⊠ Irrep[U₁](1)] = ones(1, 1)

# 	# c↑† ⊗ c↓, this is a spin 1 sector operator!!!
# 	up = TensorMap(zeros, Float64, vacuum ⊗ ph ← Rep[U₁×U₁]((1,-1)=>1) ⊗ ph)
# 	blocks(up)[Irrep[U₁](1) ⊠ Irrep[U₁](0)] = ones(1,1) / (-sqrt(2))
# 	middle = TensorMap(zeros, Float64, vacuum ⊗ ph ← vacuum ⊗ ph )
# 	blocks(middle)[Irrep[U₁](1) ⊠ Irrep[U₁](0)] = 0.5 * ones(1,1)
# 	blocks(middle)[Irrep[U₁](0) ⊠ Irrep[U₁](1)] = -0.5 * ones(1,1)
# 	down = TensorMap(zeros, Float64, vacuum ⊗ ph ← Rep[U₁×U₁]((-1,1)=>1) ⊗ ph)
# 	blocks(down)[Irrep[U₁](0) ⊠ Irrep[U₁](1)] = ones(1,1) / sqrt(2)
# 	adaga = cat(cat(up, middle, dims=3), down, dims=3)

# 	onsite_interact = TensorMap(zeros, Float64, ph ← ph)
# 	blocks(onsite_interact)[Irrep[U₁](1) ⊠ Irrep[U₁](1)]= ones(1,1)

# 	JW = TensorMap(ones, Float64, ph ← ph)
# 	blocks(JW)[Irrep[U₁](1) ⊠ Irrep[U₁](0)] = -ones(1, 1)
# 	blocks(JW)[Irrep[U₁](0) ⊠ Irrep[U₁](1)] = -ones(1, 1)

# 	occupy = TensorMap(ones, Float64, ph ← ph)
# 	blocks(occupy)[Irrep[U₁](0) ⊠ Irrep[U₁](0)] = zeros(1, 1)
# 	blocks(occupy)[Irrep[U₁](1) ⊠ Irrep[U₁](1)] = 2 * ones(1, 1)
# 	return Dict("+"=>adag, "-"=>a, "++"=>adagadag, "+-"=>adaga, "n↑n↓"=>onsite_interact, 
# 		"JW"=>JW, "n"=>occupy)
# end

# const sp_u1u1 = spinal_fermion_site_ops_u1_u1()
# const sp_u1su2 = spinal_fermion_site_ops_u1_su2()

# _select_sp(symmetry::ChargeCharge) = sp_u1u1
# _select_sp(symmetry::SpinCharge) = sp_u1su2


# function _physical_space(symmetry::FermionicSymmetry)
# 	p = _select_sp(symmetry)
# 	return space(p["+"], 2)
# end
