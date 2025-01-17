

function convert_to_mpotensor(data, leftspaces)
	@assert length(leftspaces) == size(data, 1)
	m, n = size(data)
	newdata = Matrix{Any}(undef, m, n)
	rightspaces = Vector{Union{Missing, eltype(leftspaces)}}(undef, n)
	rightspaces .= missing

	for i in 1:m
		for j in 1:n
			dj = data[i, j]
			if iszero(dj)
				newdata[i, j] = 0
			else
				newdata[i, j] = tompotensor(dj, left=leftspaces[i])
				right = space(newdata[i, j], 3)'
				if !ismissing(rightspaces[j])
					(rightspaces[j] == right) || throw(SpaceMismatch())
				else
					rightspaces[j] = right
				end
			end
			# newdata[i, j], right = _convert_to_tensor_map(dj, leftspaces[i])
			# if !ismissing(rightspaces[j])
			# 	(rightspaces[j] == right) || throw(SpaceMismatch())
			# else
			# 	rightspaces[j] = right
			# end
		end
	end
	# for i in 1:n
	# 	ismissing(rightspaces[i]) && error("rightspaces $i is missing")
	# end
	any(ismissing, rightspaces) && error("Missing spaces unresolved")
	return newdata, convert(typeof(leftspaces), rightspaces)
end
function convert_to_mpotensors(data)
	data = _prune!(data)
	# @assert !isempty(data)
	# K = length(data)
	# @assert K > 1
	leftspaces = [oneunit(spacetype(data[1][1,1]))]
	newdata = Vector{Any}(undef, length(data))
	for i in 1:length(data)
		# println("we are on orbital $i......")
		newdata[i], leftspaces = convert_to_mpotensor(data[i], leftspaces)
	end
	(length(leftspaces) == 1) || error("something wrong")
	(leftspaces[1] == oneunit(leftspaces[1])) || throw(SpaceMismatch("space_r must be vacuum"))
	return newdata
end
function _prune!(data)
	@assert !isempty(data)
	K = length(data)
	@assert K > 1
	for i in 1:K-1
		m, n = size(data[i])
		validcols = [true for l in 1:n]
		for s2 in 1:n
			if all(iszero, data[i][:, s2])
				# println("col $s2 on orbital $i will be pruned")
				validcols[s2] = false
			end
		end
		any(validcols) || error("something wrong")
		if !all(validcols)
			newdata1 = data[i][:, validcols]
			newdata2 = data[i+1][validcols, :]
			data[i], data[i+1] = newdata1, newdata2
		end
	end
	for i in K:-1:2
		m, n = size(data[i])
		validrows = [true for l in 1:m]
		for s1 in 1:m
			if all(iszero, data[i][s1, :])
				# println("row $s1 on orbital $i will be pruned")
				validrows[s1] = false
			end
		end
		any(validrows) || error("something wrong")
		if !all(validrows)
			newdata1 = data[i-1][:, validrows]
			newdata2 = data[i][validrows, :]
			data[i-1], data[i] = newdata1, newdata2
		end
	end
	return data
end

# dims1l = [1,cg+rg,cg+rg,cg2+cg*rg+rg2,cg2+cg*rg+rg2,lg,lg**2,lg,lg,1]
# dims2l = [1,rg,rg,rg2,rg2,lg+cg,lg**2+lg*cg+lg*cg+cg**2,lg+cg,lg+cg,1]

# function qcmpo(h1e::AbstractMatrix{<:Real}, h2e::AbstractArray{<:Real, 4}) 
# 	k = size(h1e, 1)
# 	@assert k == size(h1e, 2) == size(h2e, 1) == size(h2e, 2) == size(h2e, 3) == size(h2e, 4)
# 	h1e_s = _to_spin_orbitals(h1e)
# 	h2e_s = _to_spin_orbitals(h2e)
# 	partition = [[2*i-1, 2*i] for i in 1:k]
# 	return qcmpo_spin_orbitals(h1e_s, h2e_s, partition)
# end

function qcmpo(h1e::AbstractMatrix{<:Real}, h2e::AbstractArray{<:Real, 4}, symmetry::ChargeCharge)
	K = size(h1e, 1)
	partition = [[2*i-1, 2*i] for i in 1:K]
	@assert K > 1
	h1e, h2e = get_spin_orbitals(h1e, h2e) 
	h2e = remove_antisymmetric(antisymmetrize(h2e))
	mpo = Vector{Any}(undef, K)
	for i in 1:K
		sl, sc, sr = vcat(partition[1:i-1]...), partition[i], vcat(partition[i+1:end]...)
		lg, cg, rg = length(sl), length(sc), length(sr)
		# lg2, cg2, rg2 = uppertriangular_dim(lg), uppertriangular_dim(cg), uppertriangular_dim(rg)
		dimls = [1, cg+rg, cg+rg, uppertriangular_dim(cg+rg), uppertriangular_dim(cg+rg), lg, lg^2, lg, lg, 1]
		dimrs = [1, rg, rg, uppertriangular_dim(rg), uppertriangular_dim(rg), lg+cg, (lg+cg)^2, lg+cg, lg+cg, 1]
		a1 = accumulate(+, dimls)
		a2 = accumulate(+, dimrs)
		# println("site ========== $i")
		# println(sl, " ", sc, " ", sr)
		if i == 1
			# the first row
			Wj = empty_u1u1_matrix(1, a2[end]) 
			Wj[:, 1:a2[1]] = l1r1(h1e, h2e, sl, sc, sr)
			Wj[:, a2[1]+1:a2[2]] = l1r2(h1e, h2e, sl, sc, sr)
			Wj[:, a2[2]+1:a2[3]] = l1r3(h1e, h2e, sl, sc, sr)
			Wj[:, a2[3]+1:a2[4]] = l1r4(h1e, h2e, sl, sc, sr)
			Wj[:, a2[4]+1:a2[5]] = l1r5(h1e, h2e, sl, sc, sr)
			Wj[:, a2[5]+lg+1:a2[6]] = l1r7(h1e, h2e, sl, sc, sr)
			idx = [square_to_linear(lg+cg, lg+cg, idxr, idxs) for idxr in lg+1:lg+cg for idxs in lg+1:lg+cg]
			Wj[:, a2[6] .+ idx] = l1r11(h1e, h2e, sl, sc, sr)
			Wj[:, a2[7] + (lg+1):a2[8]] = l1r13(h1e, h2e, sl, sc, sr)
			Wj[:, a2[8] + (lg+1):a2[9]] = l1r15(h1e, h2e, sl, sc, sr)
			Wj[:, a2[9]+1:a2[10]] = l1r16(h1e, h2e, sl, sc, sr)
		elseif i == K
			# the last column
			Wj = empty_u1u1_matrix(a1[end], 1)
			Wj[1:a1[1], 1] = l1r16(h1e, h2e, sl, sc, sr)
			Wj[a1[1]+1:a1[1]+cg, 1] = l2r16(h1e, h2e, sl, sc, sr)
			Wj[a1[2]+1:a1[2]+cg, 1] = l4r16(h1e, h2e, sl, sc, sr)
			idx = [uppertriangular_to_linear(cg+rg, idxr, idxs) for idxr in 1:cg for idxs in idxr+1:cg]
			Wj[a1[3] .+ idx, 1] = l6r16(h1e, h2e, sl, sc, sr)
			Wj[a1[4] .+ idx, 1] = l9r16(h1e, h2e, sl, sc, sr)
			Wj[a1[5]+1:a1[6], 1]= l12r16(h1e, h2e, sl, sc, sr)
			Wj[a1[6]+1:a1[7], 1]= l13r16(h1e, h2e, sl, sc, sr)
			Wj[a1[7]+1:a1[8], 1]= l14r16(h1e, h2e, sl, sc, sr)
			Wj[a1[8]+1:a1[9], 1]= l15r16(h1e, h2e, sl, sc, sr)
			Wj[a1[9]+1:a1[10], 1]= l16r16(h1e, h2e, sl, sc, sr)
		else
			Wj = empty_u1u1_matrix(a1[end], a2[end])
			# row 1
			Wj[1:a1[1], 1:a2[1]] = l1r1(h1e, h2e, sl, sc, sr)
			Wj[1:a1[1], a2[1]+1:a2[2]] = l1r2(h1e, h2e, sl, sc, sr)
			Wj[1:a1[1], a2[2]+1:a2[3]] = l1r3(h1e, h2e, sl, sc, sr)
			Wj[1:a1[1], a2[3]+1:a2[4]] = l1r4(h1e, h2e, sl, sc, sr)
			Wj[1:a1[1], a2[4]+1:a2[5]] = l1r5(h1e, h2e, sl, sc, sr)
			Wj[1:a1[1], a2[5]+lg+1:a2[6]] = l1r7(h1e, h2e, sl, sc, sr)
			idx = [square_to_linear(lg+cg, lg+cg, idxr, idxs) for idxr in lg+1:lg+cg for idxs in lg+1:lg+cg]
			Wj[1:a1[1], a2[6] .+ idx] = l1r11(h1e, h2e, sl, sc, sr)
			Wj[1:a1[1], a2[7] + (lg+1):a2[8]] = l1r13(h1e, h2e, sl, sc, sr)
			Wj[1:a1[1], a2[8] + (lg+1): a2[9]] = l1r15(h1e, h2e, sl, sc, sr)
			Wj[1:a1[1], a2[9]+1:a2[10]] = l1r16(h1e, h2e, sl, sc, sr)
			# row 2
			Wj[a1[1]+1:a1[1]+cg, a2[9]+1:a2[10]] = l2r16(h1e, h2e, sl, sc, sr)
			# row 3
			Wj[a1[1]+cg+1:a1[2], a2[1]+1:a2[2]] = l3r2(h1e, h2e, sl, sc, sr)
			# row 4
			Wj[a1[2]+1:a1[2]+cg, a2[9]+1:a2[10]] = l4r16(h1e, h2e, sl, sc, sr)
			# row 5
			Wj[a1[2]+cg+1:a1[3], a2[2]+1:a2[3]] = l5r3(h1e, h2e, sl, sc, sr)
			# row 6
			idx1 = [uppertriangular_to_linear(cg+rg, idxr, idxs) for idxr in 1:cg for idxs in (idxr+1):cg]
			Wj[a1[3] .+ idx1, a2[9]+1:a2[10]] = l6r16(h1e, h2e, sl, sc, sr)
			# row 7
			idx2 = [uppertriangular_to_linear(cg+rg, idxr, idxs) for idxr in 1:cg for idxs in cg+1:cg+rg]
			Wj[a1[3] .+ idx2, a2[1]+1:a2[2]] = l7r2(h1e, h2e, sl, sc, sr)
			# row 8
			idx3 = [uppertriangular_to_linear(cg+rg, idxr, idxs) for idxr in cg+1:cg+rg for idxs in (idxr+1):cg+rg]
			Wj[a1[3] .+ idx3, a2[3]+1:a2[4]] = l8r4(h1e, h2e, sl, sc, sr)
			@assert [sort(union(idx1, idx2, idx3))...] == collect(1:uppertriangular_dim(cg+rg))
			# row 9
			Wj[a1[4] .+ idx1, a2[9]+1:a2[10]] = l9r16(h1e, h2e, sl, sc, sr)
		    # row 10
		    Wj[a1[4] .+ idx2, a2[2]+1:a2[3]] = l10r3(h1e, h2e, sl, sc, sr)
		    # row 11
		    Wj[a1[4] .+ idx3, a2[4]+1:a2[5]] = l11r5(h1e, h2e, sl, sc, sr)
		    # row 12
		    Wj[a1[5]+1:a1[6], a2[5]+1:a2[5]+lg] = l12r6(h1e, h2e, sl, sc, sr)
		    Wj[a1[5]+1:a1[6], a2[9]+1:a2[10]] = l12r16(h1e, h2e, sl, sc, sr)
		    # row 13
		    Wj[a1[6]+1:a1[7], a2[1]+1:a2[2]] = l13r2(h1e, h2e, sl, sc, sr)
		    Wj[a1[6]+1:a1[7], a2[2]+1:a2[3]] = l13r3(h1e, h2e, sl, sc, sr)
		    # dim = LinearIndices((lg+cg, lg+cg))
		    # idx = [dim[i, j] for i in 1:lg for j in 1:lg]
		    idx = [square_to_linear(lg+cg, lg+cg, idxr, idxs) for idxr in 1:lg for idxs in 1:lg]
		    Wj[a1[6]+1:a1[7], a2[6] .+ idx] = l13r8(h1e, h2e, sl, sc, sr)
		    Wj[a1[6]+1:a1[7], a2[9]+1:a2[10]] = l13r16(h1e, h2e, sl, sc, sr)
		    # row 14
		    Wj[a1[7]+1:a1[8], a2[1]+1:a2[2]] = l14r2(h1e, h2e, sl, sc, sr)
		    Wj[a1[7]+1:a1[8], a2[2]+1:a2[3]] = l14r3(h1e, h2e, sl, sc, sr)
		    Wj[a1[7]+1:a1[8], a2[4]+1:a2[5]] = l14r5(h1e, h2e, sl, sc, sr)
		    # idx = [dim[i, j] for i in 1:lg for j in lg+1:lg+cg]
		    idx = [square_to_linear(lg+cg, lg+cg, idxr, idxs) for idxr in 1:lg for idxs in lg+1:lg+cg]
		    Wj[a1[7]+1:a1[8], a2[6] .+ idx] = l14r9(h1e, h2e, sl, sc, sr)
		    Wj[a1[7]+1:a1[8], a2[7]+1:a2[7]+lg] = l14r12(h1e, h2e, sl, sc, sr)
		    Wj[a1[7]+1:a1[8], a2[9]+1:a2[10]] = l14r16(h1e, h2e, sl, sc, sr)
		    # row 15
		    Wj[a1[8]+1:a1[9], a2[1]+1:a2[2]] = l15r2(h1e, h2e, sl, sc, sr)
		    Wj[a1[8]+1:a1[9], a2[2]+1:a2[3]] = l15r3(h1e, h2e, sl, sc, sr)
		    Wj[a1[8]+1:a1[9], a2[3]+1:a2[4]] = l15r4(h1e, h2e, sl, sc, sr)
		    # idx = [dim[i, j] for i in lg+1:lg+cg for j in 1:lg]
		    idx = [square_to_linear(lg+cg, lg+cg, idxr, idxs) for idxr in lg+1:lg+cg for idxs in 1:lg]
		    Wj[a1[8]+1:a1[9], a2[6] .+ idx] = l15r10(h1e, h2e, sl, sc, sr)
		    Wj[a1[8]+1:a1[9], a2[8]+1:a2[8]+lg] = l15r14(h1e, h2e, sl, sc, sr)
		    Wj[a1[8]+1:a1[9], a2[9]+1:a2[10]] = l15r16(h1e, h2e, sl, sc, sr)
		    # row 16
		    Wj[a1[9]+1:a1[10], a2[9]+1:a2[10]] = l16r16(h1e, h2e, sl, sc, sr)
		end
		mpo[i] = Wj
	end
	mpo = convert_to_mpotensors(mpo)
	mpo = [SparseMPOTensor(m, Float64, u1u1_space) for m in mpo]
	return MPOHamiltonian(mpo)
end
