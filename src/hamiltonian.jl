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

