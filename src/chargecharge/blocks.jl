# function spin_half_matrices()
# 	s_SP = Array{Float64, 2}([0 0; 1 0])
# 	s_SM = Array{Float64, 2}([0 1; 0 0])
# 	s_Z = Array{Float64, 2}([-1 0; 0 1])
# 	s_x = s_SP+s_SM
# 	s_y = -im*(s_SP-s_SM)
# 	return Dict("x"=>s_x, "y"=>s_y, "z"=>s_Z, "+"=>s_SP, "-"=>s_SM)
# end
function _spinless_operators()
	p = GeneralHamiltonians.spin_half_matrices()
	p = Dict(k=>convert(AbelianMatrix, v) for (k, v) in p)
	sp, sm, z = p["+"], p["-"], p["z"]
	I2 = one(z)
	JW = -z
	return sp, sm, I2, JW
end

const u1u1_space = Rep[U₁×U₁]((1, 1)=>1, (1,0)=>1, (0,1)=>1, (0,0)=>1)
const adag, a, I2, JW2 = _spinless_operators()
const adagup = kron(adag, I2)
const aup = adagup'
const adagdown = kron(JW2, adag)
const adown = adagdown'
const I4 = kron(I2, I2)
const JW4 = kron(JW2, JW2)

function empty_u1u1_matrix(m::Int, n::Int)
	tmp = Matrix{AbelianMatrix{Rep[U₁×U₁], Float64, Irrep[U₁×U₁]}}(undef, m, n)
	for i in 1:length(tmp)
		tmp[i] = AbelianMatrix(Float64, u1u1_space, Dict())
	end
	return tmp
end

idnC(sc) = I4
sgnC(sc) = JW4
function sqC(sc, ic::Int, iop::Bool)
	@assert (length(sc) == 2) && (sc[2] == sc[1]+1)
	@assert (ic == 1) || (ic == 2)
	if ic == 1
		return iop ? adagup : aup
	else
		return iop ? adagdown : adown
	end
end

# row major traversal
function uppertriangular_to_linear(m::Int, idx1::Int, idx2::Int)
	@assert idx1 < idx2 <= m
	# column major traversal
	# return div((idx2-1) * (idx2-2), 2) + idx1
	return div((idx1-1) * (2*m - idx1), 2) + idx2 - idx1
end
uppertriangular_dim(c::Int) = div(c * (c-1), 2)
function square_to_linear(m::Int, n::Int, idx1::Int, idx2::Int)
	@assert (idx1 <= m) && (idx2 <= n)
	return (idx1-1) * n + idx2
end

# sc must be a pair of nearest-neighbour sites
function l1r1(h1e, h2e, sl, sc, sr)
	tmp = empty_u1u1_matrix(1, 1)
	tmp[1, 1] = idnC(sc)
	return tmp
end 
function l1r2(h1e, h2e, sl, sc, sr)
	tmp = empty_u1u1_matrix(1, length(sr))
	for (idxp, orbp) in enumerate(sc)
		op_p = sqC(sc, idxp, true)
		for (idxr, orbr) in enumerate(sc)
			op_r = sqC(sc, idxr, false)
			for (idxs, orbs) in enumerate(sc)
				op_s = sqC(sc, idxs, false)
				if orbr < orbs
					op = op_p * op_r * op_s
					for (idxq, orbq) in enumerate(sr)
						tmp[1, idxq] += h2e[orbp,orbq,orbr,orbs]*op
					end
				end
			end
		end
	end
	sgn = sgnC(sc)
	for (idxq, orbq) in enumerate(sr)
		tmp[1, idxq] = tmp[1, idxq] * sgn
	end
	return tmp
end
function l1r3(h1e, h2e, sl, sc, sr)
	phy = 2^(length(sc))
	tmp = empty_u1u1_matrix(1, length(sr))
	sgn = sgnC(sc)
	for (idxp, orbp) in enumerate(sc)
		op_p = sqC(sc, idxp, true)
		for (idxq, orbq) in enumerate(sr)
			tmp[1, idxq] += h1e[orbp, orbq] * op_p
		end
		for (idxq, orbq) in enumerate(sc)
			op_q = sqC(sc, idxq, true)
			if orbp < orbq
				for (idxr, orbr) in enumerate(sc)
					op_r = sqC(sc, idxr, false)
					op = op_p * op_q * op_r
					for (idxs, orbs) in enumerate(sr)
						tmp[1, idxs] += h2e[orbp,orbq,orbr,orbs] * op
					end
				end
			end
		end
	end
	for (idxq, orbq) in enumerate(sr)
		tmp[1, idxq] = tmp[1, idxq] * sgn
	end
	return tmp
end
function l1r4(h1e, h2e, sl, sc, sr)
	rg = length(sr)
	tmp = empty_u1u1_matrix(1, uppertriangular_dim(rg))
	for (idxr, orbr) in enumerate(sc)
		op_r = sqC(sc, idxr, false)
		for (idxs, orbs) in enumerate(sc)
			op_s = sqC(sc, idxs, false)
			if orbr < orbs
				op = op_r * op_s
				for (idxp, orbp) in enumerate(sr)
					for (idxq, orbq) in enumerate(sr)
						if orbp < orbq
							ipq = uppertriangular_to_linear(rg, idxp, idxq)
							tmp[1, ipq] += h2e[orbp,orbq,orbr,orbs] * op
						end
					end
				end
			end
		end
	end	
	return tmp
end
function l1r5(h1e, h2e, sl, sc, sr)
	rg = length(sr)
	tmp = empty_u1u1_matrix(1, uppertriangular_dim(rg))
	for (idxr, orbr) in enumerate(sc)
		op_r = sqC(sc, idxr, true)
		for (idxs, orbs) in enumerate(sc)
			op_s = sqC(sc, idxs, true)
			if orbr < orbs
				op = op_r * op_s
				for (idxp, orbp) in enumerate(sr)
					for (idxq, orbq) in enumerate(sr)
						if orbp < orbq
							ipq = uppertriangular_to_linear(rg, idxp, idxq)
							tmp[1, ipq] += h2e[orbr,orbs,orbp,orbq] * op
						end
					end
				end
			end
		end
	end	
	return tmp
end
function l1r7(h1e, h2e, sl, sc, sr)
	nc = length(sc)
	tmp = empty_u1u1_matrix(1, nc)
	sgn = sgnC(sc)
	for (idxr, orbr) in enumerate(sc)
		tmp[1, idxr] = -sqC(sc, idxr, false) * sgn
	end
	return tmp
end
function l1r11(h1e, h2e, sl, sc, sr)
	nc = length(sc)
	tmp = empty_u1u1_matrix(1, nc * nc)
	for (idxr, orbr) in enumerate(sc)
		for (idxs, orbs) in enumerate(sc)
			idx = square_to_linear(nc, nc, idxr, idxs)
			tmp[1, idx] = -sqC(sc, idxr, true) * sqC(sc, idxs, false)
		end
	end
	return tmp
end
function l1r13(h1e, h2e, sl, sc, sr)
	nc = length(sc)
	tmp = empty_u1u1_matrix(1, nc)
	sgn = sgnC(sc)
	for (idxr, orbr) in enumerate(sc)
		tmp[1, idxr] = sqC(sc, idxr, true) * sgn
	end
	return tmp
end
function l1r15(h1e, h2e, sl, sc, sr)
	nc = length(sc)
	tmp = empty_u1u1_matrix(1, nc)
	sgn = sgnC(sc)
	for (idxr, orbr) in enumerate(sc)
		tmp[1, idxr] = sqC(sc, idxr, false) * sgn
	end
	return tmp
end
function l1r16(h1e, h2e, sl, sc, sr)
	tmp = empty_u1u1_matrix(1, 1)
	for (idxp, orbp) in enumerate(sc)
		op_p = sqC(sc, idxp, true)
		for (idxq, orbq) in enumerate(sc)
			op_q = sqC(sc, idxq, false)
			op = op_p * op_q
			tmp[1,1] += h1e[orbp,orbq]*op
		end
	end
	for (idxp, orbp) in enumerate(sc)
		op_p = sqC(sc, idxp, true)
		for (idxq, orbq) in enumerate(sc)
			op_q = sqC(sc, idxq, true)
			if orbp < orbq
				op_pq = op_p * op_q
				for (idxr, orbr) in enumerate(sc)
					op_r = sqC(sc, idxr, false)
					for (idxs, orbs) in enumerate(sc)
						op_s = sqC(sc, idxs, false)
						if orbr < orbs
							op_rs = op_r * op_s
							op = op_pq * op_rs
							tmp[1,1] += h2e[orbp,orbq,orbr,orbs]*op
						end
					end
				end
			end
		end
	end
	return tmp
end
function l2r16(h1e, h2e, sl, sc, sr)
	nc = length(sc)
	tmp = empty_u1u1_matrix(nc, 1)
	for (idxr, orbr) in enumerate(sc)
		tmp[idxr, 1] = sqC(sc, idxr, true)
	end
	return tmp
end
function l3r2(h1e, h2e, sl, sc, sr)
	nr = length(sr)
	tmp = empty_u1u1_matrix(nr, nr)
	sgn = sgnC(sc)
	for (idxr, orbr) in enumerate(sr)
		tmp[idxr, idxr] = sgn
	end
	return tmp
end
function l4r16(h1e, h2e, sl, sc, sr)
	nc = length(sc)
	tmp = empty_u1u1_matrix(nc, 1)
	for (idxr, orbr) in enumerate(sc)
		tmp[idxr, 1] = sqC(sc, idxr, false)
	end
	return tmp
end
function l5r3(h1e, h2e, sl, sc, sr)
	nr = length(sr)
	tmp = empty_u1u1_matrix(nr, nr)
	sgn = sgnC(sc)
	for (idxr, orbr) in enumerate(sr)
		tmp[idxr, idxr] = sgn
	end
	return tmp
end
function l6r16(h1e, h2e, sl, sc, sr)
	nc = length(sc)
	tmp = empty_u1u1_matrix(uppertriangular_dim(nc), 1)
	for (idxr, orbr) in enumerate(sc)
		for (idxs, orbs) in enumerate(sc)
			if orbr < orbs
				idx = uppertriangular_to_linear(nc, idxr, idxs)
				tmp[idx, 1] = sqC(sc, idxr, true) * sqC(sc, idxs, true)
			end
		end
	end
	return tmp
end
function l7r2(h1e, h2e, sl, sc, sr)
	nc = length(sc)
	nr = length(sr)
	sgn = sgnC(sc)
	tmp = empty_u1u1_matrix(nc * nr, nr)
	for (idxr, orbr) in enumerate(sc)
		for (idxs, orbs) in enumerate(sr)
			idx = square_to_linear(nc, nr, idxr, idxs)
			tmp[idx, idxs] = sqC(sc, idxr, true) * sgn
		end
	end
	return tmp
end
function l8r4(h1e, h2e, sl, sc, sr)
	nr = length(sr)
	nr2 = uppertriangular_dim(nr)	
	tmp = empty_u1u1_matrix(nr2, nr2)
	idn = idnC(sc)
	for irs in 1:nr2
		tmp[irs, irs] = idn
	end
	return tmp
end
function l9r16(h1e, h2e, sl, sc, sr)
	nc = length(sc)
	nr = length(sr)
	tmp = empty_u1u1_matrix(uppertriangular_dim(nc), 1)
	for (idxr, orbr) in enumerate(sc)
		for (idxs, orbs) in enumerate(sc)
			if orbr < orbs
				idx = uppertriangular_to_linear(nc, idxr, idxs)
				tmp[idx, 1] = sqC(sc, idxr, false) * sqC(sc, idxs, false)
			end
		end
	end
	return tmp
end
function l10r3(h1e, h2e, sl, sc, sr)
	nc = length(sc)
	nr = length(sr)
	tmp = empty_u1u1_matrix(nc *nr, nr)
	sgn = sgnC(sc)
	for (idxr, orbr) in enumerate(sc)
		for (idxs, orbs) in enumerate(sr)
			idx = square_to_linear(nc, nr, idxr, idxs)
			tmp[idx, idxs] = sqC(sc, idxr, false) * sgn
		end
	end
	return tmp
end
function l11r5(h1e, h2e, sl, sc, sr)
	nc = length(sc)
	nr = length(sr)
	nr2 = uppertriangular_dim(nr)	
	tmp = empty_u1u1_matrix(nr2, nr2)
	idn = idnC(sc)
	for irs in 1:nr2
		tmp[irs, irs] = idn
	end
	return tmp
end
function l12r6(h1e, h2e, sl, sc, sr)
	nl = length(sl)
	tmp = empty_u1u1_matrix(nl, nl)
	sgn = sgnC(sc)
	for (idxr, orbr) in enumerate(sl)
		tmp[idxr, idxr] = sgn
	end
	return tmp
end
function l12r16(h1e, h2e, sl, sc, sr)
	nl = length(sl)
	tmp = empty_u1u1_matrix(nl, 1)
	sgn = sgnC(sc)
	for (idxp, orbp) in enumerate(sc)
		op_p = sqC(sc, idxp, true)
		for (idxq, orbq) in enumerate(sl)
			tmp[idxq, 1] += h1e[orbp,orbq]*op_p
		end
	end
	return tmp
end
function l13r2(h1e, h2e, sl, sc, sr)
	nl = length(sl)
	nr = length(sr)
	tmp = empty_u1u1_matrix(nl * nl, nr)
	sgn = sgnC(sc)
	for (idxs, orbs) in enumerate(sc)
		op_s = sqC(sc, idxs, false)
		op = op_s * sgn
		for (idxp, orbp) in enumerate(sl)
			for (idxr, orbr) in enumerate(sl)
				for (idxq, orbq) in enumerate(sr)
					# idx = dim[idxp, idxr]
					idx = square_to_linear(nl, nl, idxp, idxr)
					tmp[idx, idxq] += -h2e[orbp,orbq,orbr,orbs]*op
				end
			end
		end
	end
	return tmp
end
function l13r3(h1e, h2e, sl, sc, sr)
	nl = length(sl)
	nr = length(sr)
	tmp = empty_u1u1_matrix(nl*nl, nr)
	sgn = sgnC(sc)
	for (idxq, orbq) in enumerate(sc)
		op_q = sqC(sc, idxq, true)
		op = op_q * sgn
		for (idxp, orbp) in enumerate(sl)
			for (idxr, orbr) in enumerate(sl)
				for (idxs, orbs) in enumerate(sr)
					# idx = dim[idxp, idxr]
					idx = square_to_linear(nl, nl, idxp, idxr)
					tmp[idx, idxs] += h2e[orbp,orbq,orbr,orbs]*op
				end
			end
		end
	end
	return tmp
end
function l13r8(h1e, h2e, sl, sc, sr)
	nl = length(sl)
	tmp = empty_u1u1_matrix(nl*nl, nl*nl)
	idn = idnC(sc)
	for i in 1:nl*nl
		tmp[i, i] = idn
   	end
   	return tmp
end
function l13r16(h1e, h2e, sl, sc, sr)
	nl = length(sl)
   	tmp = empty_u1u1_matrix(nl*nl, 1)
   	for (idxq, orbq) in enumerate(sc)
   		op_q = sqC(sc, idxq, true)
   		for (idxs, orbs) in enumerate(sc)
   			op_s = sqC(sc, idxs, false)
   			op = op_q * op_s
   			for (idxp, orbp) in enumerate(sl)
   				for (idxr, orbr) in enumerate(sl)
   					idx = square_to_linear(nl, nl, idxp, idxr)
   					tmp[idx, 1] += h2e[orbp,orbq,orbr,orbs]*op
   				end
   			end
   		end
   	end
   	return tmp
end
function l14r2(h1e, h2e, sl, sc, sr)
	nl = length(sl)
	nr = length(sr)	
	sgn = sgnC(sc)
	tmp = empty_u1u1_matrix(nl, nr)
	for (idxr, orbr) in enumerate(sc)
		op_r = sqC(sc, idxr, false)
		for (idxs, orbs) in enumerate(sc)
			op_s = sqC(sc, idxs, false)
			if orbr < orbs
				op = op_r * op_s * sgn
				for (idxp, orbp) in enumerate(sl)
					for (idxq, orbq) in enumerate(sr)
						tmp[idxp, idxq] += h2e[orbp,orbq,orbr,orbs]*op
					end
				end
			end
		end
	end
	return tmp
end
function l14r3(h1e, h2e, sl, sc, sr)
	nl = length(sl)
	nr = length(sr)	
	sgn = sgnC(sc)
	tmp = empty_u1u1_matrix(nl, nr)
	for (idxq, orbq) in enumerate(sc)
		op_q = sqC(sc, idxq, true)		
		for (idxr, orbr) in enumerate(sc)
			op_r = sqC(sc, idxr, false)
			op = op_q * op_r * sgn
			for (idxp, orbp) in enumerate(sl)
				for (idxs, orbs) in enumerate(sr)
					tmp[idxp, idxs] += h2e[orbp,orbq,orbr,orbs]*op
				end
			end
		end
	end	
	return tmp
end
function l14r5(h1e, h2e, sl, sc, sr)
	nl = length(sl)
	nr = length(sr)	
	tmp = empty_u1u1_matrix(nl, uppertriangular_dim(nr))
	for (idxq, orbq) in enumerate(sc)
		op_q = sqC(sc, idxq, true)
		for (idxp, orbp) in enumerate(sl)
			for (idxr, orbr) in enumerate(sr)
				for (idxs, orbs) in enumerate(sr)
					if orbr < orbs
						idx = uppertriangular_to_linear(nr, idxr, idxs)
						tmp[idxp, idx] += h2e[orbp,orbq,orbr,orbs]*op_q
					end
				end
			end
		end
	end
	return tmp
end
function l14r9(h1e, h2e, sl, sc, sr)
	nl = length(sl)
	nc = length(sc)
	tmp = empty_u1u1_matrix(nl, nl*nc)
	for (idxr, orbr) in enumerate(sc)
		op_r = sqC(sc, idxr, false)
		for (idxp, orbp) in enumerate(sl)
			# idx = dim[idxp, idxr]
			idx = square_to_linear(nl, nc, idxp, idxr)
			tmp[idxp, idx] = -op_r
		end
	end
	return tmp
end
function l14r12(h1e, h2e, sl, sc, sr)
	nl = length(sl)
	tmp = empty_u1u1_matrix(nl, nl)
	sgn = sgnC(sc)
	for (idxq, orbq) in enumerate(sl)
		tmp[idxq, idxq] = sgn
	end
	return tmp
end
function l14r16(h1e, h2e, sl, sc, sr)
	nl = length(sl)
   	tmp = empty_u1u1_matrix(nl, 1)
   	for (idxq, orbq) in enumerate(sc)
   		op_q = sqC(sc, idxq, true)
   		for (idxr, orbr) in enumerate(sc)
   			op_r = sqC(sc, idxr, false)
   			for (idxs, orbs) in enumerate(sc)
   				op_s = sqC(sc, idxs, false)
   				if orbr < orbs
   					op = op_q * op_r * op_s
   					for (idxp, orbp) in enumerate(sl)
   						tmp[idxp, 1] += h2e[orbp,orbq,orbr,orbs]*op
   					end
   				end
   			end
   		end
   	end
   	return tmp
end
function l15r2(h1e, h2e, sl, sc, sr)
	nl = length(sl)
   	nr = length(sr)	
   	tmp = empty_u1u1_matrix(nl, nr)
	sgn = sgnC(sc)
	for (idxp, orbp) in enumerate(sc)
		op_p = sqC(sc, idxp, true)
		for (idxs, orbs) in enumerate(sc)
			op_s = sqC(sc, idxs, false)
			op = op_p * op_s * sgn
			for (idxr, orbr) in enumerate(sl)
				for (idxq, orbq) in enumerate(sr)
					tmp[idxr, idxq] += -h2e[orbp,orbq,orbr,orbs]*op
				end
			end
		end
	end
	return tmp
end
function l15r3(h1e, h2e, sl, sc, sr)
	nl = length(sl)
   	nr = length(sr)	
   	tmp = empty_u1u1_matrix(nl, nr)
   	sgn = sgnC(sc)
   	for (idxq, orbq) in enumerate(sc)
   		op_q = sqC(sc, idxq, true)
   		for (idxp, orbp) in enumerate(sc)
   			op_p = sqC(sc, idxp, true)
   			if orbp < orbq
   				op = op_p * op_q * sgn
   				for (idxr, orbr) in enumerate(sl)
   					for (idxs, orbs) in enumerate(sr)
   						tmp[idxr, idxs] += h2e[orbp,orbq,orbr,orbs]*op
   					end
   				end
   			end
   		end
   	end
   	return tmp
end
function l15r4(h1e, h2e, sl, sc, sr)
	nl = length(sl)
   	nr = length(sr)	
	tmp = empty_u1u1_matrix(nl, uppertriangular_dim(nr))
	for (idxs, orbs) in enumerate(sc)
		op_s = sqC(sc, idxs, false)
		for (idxr, orbr) in enumerate(sl)
			for (idxp, orbp) in enumerate(sr)
				for (idxq, orbq) in enumerate(sr)
					if orbp < orbq
						idx = uppertriangular_to_linear(nr, idxp, idxq)
						tmp[idxr, idx] += h2e[orbp,orbq,orbr,orbs]*op_s
					end
				end
			end
		end
	end
	return tmp
end
function l15r10(h1e, h2e, sl, sc, sr)
	nl = length(sl)
	nc = length(sc)
	tmp = empty_u1u1_matrix(nl, nc*nl)
	for (idxp, orbp) in enumerate(sc)
		op_p = sqC(sc, idxp, true)
		for (idxr, orbr) in enumerate(sl)
			idx = square_to_linear(nc, nl, idxp, idxr)
			tmp[idxr, idx] = op_p
		end
	end
	return tmp
end
function l15r14(h1e, h2e, sl, sc, sr)
	nl = length(sl)
   	tmp = empty_u1u1_matrix(nl, nl)
   	sgn = sgnC(sc)
   	for (idxq, orbq) in enumerate(sl)
   		tmp[idxq, idxq] = sgn
   	end
   	return tmp
end
function l15r16(h1e, h2e, sl, sc, sr)
	nl = length(sl)
   	tmp = empty_u1u1_matrix(nl, 1)
   	for (idxp, orbp) in enumerate(sc)
   		op_p = sqC(sc, idxp, true)
   		for (idxq, orbq) in enumerate(sc)
   			op_q = sqC(sc, idxq, true)
   			if orbp < orbq
   				for (idxs, orbs) in enumerate(sc)
   					op_s = sqC(sc, idxs, false)
   					op = op_p * op_q * op_s
   					for (idxr, orbr) in enumerate(sl)
   						tmp[idxr, 1] += h2e[orbp,orbq,orbr,orbs]*op
   					end
   				end
   			end
   		end
   	end
   	return tmp
end
l16r16(h1e, h2e, sl, sc, sr) = l1r1(h1e, h2e, sl, sc, sr)



