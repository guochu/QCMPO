
using JSON

function read_data(pathname)
	data = JSON.parsefile(pathname)
	E0 = data["E0"]
	L = data["L"]
	t = data["t"]
	t = [t...]
	v = data["v"]
	v = [v...]
	return E0, reshape(t, (L, L)), reshape(v, (L, L, L, L))
end

read_lih_data() = read_data("lih.json")

const LiH_FCI_ENERGY = -7.78446028003123

@testset "LiH ground state with ED" begin
	E0, t, v = read_lih_data()

	h = hamiltonian(t, 0.5*v, ChargeCharge())
	mpo = MPO(h)
	sector = ChargeChargeSector(up=2, down=2)
	eigvalue, eigvector = ground_state(mpo, ED(), right=space(sector))
	energy = eigvalue + E0
	@test abs(energy - LiH_FCI_ENERGY) < 1.0e-8

	h = hamiltonian(t, 0.5*v, SpinCharge())
	mpo = MPO(h)
	sector = SpinChargeSector(spin=0, charge=4)
	eigvalue, eigvector = ground_state(mpo, ED(), right=space(sector))
	energy = eigvalue + E0
	@test abs(energy - LiH_FCI_ENERGY) < 1.0e-8

	mpo = qcmpo(t, 0.5 * v, ChargeCharge())
	sector = ChargeChargeSector(up=2, down=2)
	eigvalue, eigvector = ground_state(mpo, ED(), right=space(sector))
	energy = eigvalue + E0
	@test abs(energy - LiH_FCI_ENERGY) < 1.0e-8
end

