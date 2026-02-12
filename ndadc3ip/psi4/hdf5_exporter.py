"""
hdf5_exporter.py

Export Psi4 SCF results into an HDF5 file in the same format
that Serenity's WriteDataTask produces for TheADCcode.

Assumptions:
- C1 symmetry (no symmetry): nSym = 1, irreps[:] = 1
- Closed-shell RHF
- Cartesian basis (or at least, a consistent mapping of AO polynomials)
"""

import numpy as np
import h5py
import psi4


def double_factorial(n: int) -> int:
    """Double factorial n!! with the convention (-1)!! = 1, 0!! = 1."""
    if n <= 0:
        return 1
    result = 1
    k = n
    while k > 1:
        result *= k
        k -= 2
    return result


def _enumerate_cartesian_exponents(L: int):
    """
    Generate (aX, aY, aZ) for Cartesian Gaussians with total angular momentum L
    in the same order as Serenity:

        for (int aX = l; aX >= 0; --aX)
          for (int aY = l - aX; aY >= 0; --aY)
            int aZ = l - aY - aX;

    """
    exps = []
    for aX in range(L, 0 - 1, -1):
        for aY in range(L - aX, 0 - 1, -1):
            aZ = L - aX - aY
            exps.append((aX, aY, aZ))
    return exps


def _build_ao_basis_serenity_style(basis: psi4.core.BasisSet, mol: psi4.core.Molecule):
    """
    Reproduce Serenity::prepAOBasis using Psi4's BasisSet & GaussianShell.

    Returns:
        nAO          (int)
        maxNmbCC     (int)
        ncc          (np.int32, shape (nAO,))
        cc           (np.float64, shape (nAO * maxNmbCC,))
        alpha        (np.float64, shape (nAO * maxNmbCC,))
        center       (np.int32, shape (nAO,))
        polynomial   (np.int32, shape (3 * nAO,))
    """
    nshell = basis.nshell()
    nbf = basis.nbf()

    # First pass: compute nAO and maxNmbCC
    nAO = 0
    maxNmbCC = 0
    for ish in range(nshell):
        sh = basis.shell(ish)
        L = sh.am
        nfunc = sh.nfunction  # number of Cartesian (or pure) functions in this shell
        nprim = sh.nprimitive
        nAO += nfunc
        if nprim > maxNmbCC:
            maxNmbCC = nprim

    if nAO != nbf:
        raise RuntimeError(f"Internal consistency error: nAO={nAO} but nbf={nbf}")

    # Allocate arrays (same layout as Serenity)
    ncc = np.zeros(nAO, dtype=np.int32)
    cc = np.zeros(nAO * maxNmbCC, dtype=np.float64)
    alpha = np.zeros(nAO * maxNmbCC, dtype=np.float64)
    center = np.zeros(nAO, dtype=np.int32)
    polynomial = np.zeros(3 * nAO, dtype=np.int32)

    # Build AO basis info
    ao_idx = 0
    for ish in range(nshell):
        sh = basis.shell(ish)
        L = sh.am
        nfunc = sh.nfunction
        nprim = sh.nprimitive

        # Primitive exponents and contraction coefficients
        # Psi4 docs: GaussianShell.exp(ip), GaussianShell.coef(ip)
        exps = [sh.exp(ip) for ip in range(nprim)]
        contr = [sh.coef(ip) for ip in range(nprim)]

        # Center index: Psi4 provides shell_to_center(ish) (0-based atom index)
        atom_idx = basis.shell_to_center(ish)  # 0-based
        # We will store 1-based index (Serenity uses atom_indx + 1)
        atom_idx_hdf5 = atom_idx + 1

        # Cartesian polynomial exponents in Serenity order
        exps_cart = _enumerate_cartesian_exponents(L)
        if len(exps_cart) < nfunc:
            raise RuntimeError(
                f"Not enough Cartesian functions for shell L={L}: "
                f"have {len(exps_cart)}, need {nfunc}"
            )

        # Loop over AO functions in this shell
        for ifn in range(nfunc):
            aX, aY, aZ = exps_cart[ifn]

            # polynomial[3*indx + 0..2] = (aX, aY, aZ)
            polynomial[3 * ao_idx + 0] = aX
            polynomial[3 * ao_idx + 1] = aY
            polynomial[3 * ao_idx + 2] = aZ

            # Assign number of contractions and center index
            ncc[ao_idx] = nprim
            center[ao_idx] = atom_idx_hdf5

            # Normalization factor (same as Serenity)
            # norm = [ (2aX-1)!! (2aY-1)!! (2aZ-1)!! ]^(-1/2) * [ (2L-1)!! ]^(1/2)
            norm = (
                (
                    double_factorial(2 * aX - 1)
                    * double_factorial(2 * aY - 1)
                    * double_factorial(2 * aZ - 1)
                )
                ** (-0.5)
                * double_factorial(2 * L - 1) ** 0.5
            )

            # Fill cc and alpha, padded up to maxNmbCC
            for icc in range(nprim):
                cc[ao_idx * maxNmbCC + icc] = norm * contr[icc]
                alpha[ao_idx * maxNmbCC + icc] = exps[icc]

            ao_idx += 1

    if ao_idx != nAO:
        raise RuntimeError(f"AO enumeration mismatch: ao_idx={ao_idx}, nAO={nAO}")

    return nAO, maxNmbCC, ncc, cc, alpha, center, polynomial


def export_adc_hdf5(
    wfn: psi4.core.Wavefunction,
    mol: psi4.core.Molecule = None,
    filename: str = "data.hdf5",
    iActOrb: int = None,
    fActOrb: int = None,
    eri_threshold: float = 1.0e-10,
):
    """
    Export Psi4 SCF data to an HDF5 file in Serenity/ADC format (C1 symmetry).

    Parameters
    ----------
    wfn : psi4.core.Wavefunction
        Wavefunction from psi4.energy("SCF", return_wfn=True).
    mol : psi4.core.Molecule, optional
        Molecule; if None, wfn.molecule() is used.
    filename : str
        Output HDF5 filename.
    iActOrb, fActOrb : int
        Active-space MO indices [iActOrb, fActOrb) in canonical order.
        Defaults: full space [0, nmo).
    eri_threshold : float
        Screening threshold for MO integrals (like settings.eriThreshold).
    """
    if mol is None:
        mol = wfn.molecule()

    basis = wfn.basisset()
    mints = psi4.core.MintsHelper(basis)

    nbf = basis.nbf()
    nmo = wfn.nmo()
    nocc = wfn.nalpha()  # closed-shell assumption

    # Total SCF energy (Serenity's "energy")
    energy = float(wfn.energy())

    # Active space (like Serenity's default: full basis)
    if iActOrb is None:
        iActOrb = 0
    else:
        iActOrb = iActOrb - 1
    if fActOrb is None:
        fActOrb = nmo

    if not (0 <= iActOrb < fActOrb <= nmo):
        raise ValueError(f"Invalid active space [{iActOrb}, {fActOrb}) for nmo={nmo}")

    nActOrbs = fActOrb - iActOrb
    nOccActOrbs = min(nocc - iActOrb, nActOrbs)
    if nOccActOrbs < 0:
        raise ValueError("Active space starts above HOMO, adjust iActOrb/fActOrb.")

    print(f"[hdf5_exporter] nbf        = {nbf}")
    print(f"[hdf5_exporter] nmo        = {nmo}")
    print(f"[hdf5_exporter] nocc       = {nocc}")
    print(f"[hdf5_exporter] iActOrb    = {iActOrb}")
    print(f"[hdf5_exporter] fActOrb    = {fActOrb}")
    print(f"[hdf5_exporter] nActOrbs   = {nActOrbs}")
    print(f"[hdf5_exporter] nOccActOrbs= {nOccActOrbs}")

    # === Orbital energies (like orbitalEnergies.segment(_iActOrb, _nActOrbs)) ===
    eps_all = wfn.epsilon_a().to_array(dense=True)
    energies = eps_all[iActOrb:fActOrb]

    # === Occupations (like occs.segment) ===
    occs = np.zeros(nActOrbs, dtype=float)
    occs[:nOccActOrbs] = 2.0

    # === Irreps (Serenity: all ones, C1 symmetry) ===
    irreps = np.ones(nActOrbs, dtype=np.int32)

    # === Coefficients (prepCoefficients) ===
    # _coeffs(j,i) in Serenity == C_AO→MO(j,i) here
    C_all = np.array(wfn.Ca().to_array(dense=True))  # (nbf, nmo)
    C_act = C_all[:, iActOrb:fActOrb]                # (nbf, nActOrbs)

    coeffs_flat = np.zeros(nActOrbs * nbf, dtype=float)
    for i in range(nActOrbs):
        for j in range(nbf):
            # i_ = i - iActOrb, but we've already shifted C_act, so i_ = i
            # coeffs(i_ * nBasisFunc + j) = (*_coeffs)(j, i) in Serenity
            coeffs_flat[i * nbf + j] = C_act[j, i]

    # === Two-electron integrals in MO (prepERIS analogue) ===
    print("[hdf5_exporter] Building AO ERIs...")
    eri_ao = np.array(mints.ao_eri())  # (nbf, nbf, nbf, nbf)

    print("[hdf5_exporter] Transforming AO ERIs to MO (active space)...")
    eri_mo = np.einsum(
        "ijkl,ip,jq,kr,ls->pqrs",
        eri_ao, C_act, C_act, C_act, C_act,
        optimize=True,
    )  # shape (nActOrbs, nActOrbs, nActOrbs, nActOrbs)

    print(f"[hdf5_exporter] Screening MO ERIs with threshold {eri_threshold:.1e}")
    ints = []
    idx = []
    for i in range(nActOrbs):
        for j in range(i, nActOrbs):
            for k in range(i, nActOrbs):
                for l in range(k, nActOrbs):
                    val = eri_mo[i, j, k, l]
                    if abs(val) >= eri_threshold:
                        ints.append(val)
                        idx.extend([i, j, k, l])

    integrals = np.array(ints, dtype=float)
    indices = np.array(idx, dtype=np.int32)
    nIntegrals = integrals.size

    print(f"[hdf5_exporter] Non-zero integrals stored: {nIntegrals}")

    # === Geometry (prepGeometry) ===
    nAtoms = mol.natom()
    coords = np.zeros(3 * nAtoms, dtype=float)
    Zs = np.zeros(nAtoms, dtype=float)
    for ia in range(nAtoms):
        Z = mol.Z(ia)
        x = mol.x(ia)
        y = mol.y(ia)
        z = mol.z(ia)
        Zs[ia] = float(Z)
        coords[3 * ia + 0] = x
        coords[3 * ia + 1] = y
        coords[3 * ia + 2] = z
    Enuc = mol.nuclear_repulsion_energy()

    # === sAO and dipole matrices (like getOverlapIntegrals/getDipoleLengths) ===
    S_ao = np.array(mints.ao_overlap())  # (nbf, nbf)
    Dx, Dy, Dz = mints.ao_dipole()
    Dx_arr = np.array(Dx)
    Dy_arr = np.array(Dy)
    Dz_arr = np.array(Dz)

    # === AO basis block (prepAOBasis analogue) ===
    print("[hdf5_exporter] Building AO basis block (nAO, ncc, cc, alpha, center, polynomial)...")
    nAO, maxNmbCC, ncc, cc, alpha, center, polynomial = _build_ao_basis_serenity_style(basis, mol)
    if nAO != nbf:
        raise RuntimeError(f"AO mismatch: nAO={nAO}, nbf={nbf}")

    # === Write HDF5 exactly in Serenity style ===
    print(f"[hdf5_exporter] Writing HDF5 file: {filename}")
    with h5py.File(filename, "w") as f:
        # 1) energies, occs, irreps, coefficients
        f.create_dataset("energies", data=energies)
        f.create_dataset("occs", data=occs)
        f.create_dataset("irreps", data=irreps)
        f.create_dataset("coefficients", data=coeffs_flat)

        # 2) integrals and indices
        f.create_dataset("integrals", data=integrals)
        f.create_dataset("indices", data=indices)

        # 3) AO basis block
        f.attrs["nAO"] = int(nAO)
        f.attrs["maxNmbCC"] = int(maxNmbCC)
        f.create_dataset("ncc", data=ncc)
        f.create_dataset("cc", data=cc)
        f.create_dataset("alpha", data=alpha)
        f.create_dataset("center", data=center)
        f.create_dataset("polynomial", data=polynomial)

        # 4) AO overlap and dipole matrices
        f.create_dataset("sAO", data=S_ao.flatten())
        f.create_dataset("dip_x", data=Dx_arr.flatten())
        f.create_dataset("dip_y", data=Dy_arr.flatten())
        f.create_dataset("dip_z", data=Dz_arr.flatten())

        # 5) Geometry
        f.attrs["nAtoms"] = int(nAtoms)
        f.attrs["Enuc"] = float(Enuc)
        f.create_dataset("geometry", data=coords)
        f.create_dataset("Zs", data=Zs)

        # 6) Attributes (C1 symmetry, Serenity-style)
        f.attrs["nSym"] = 1
        f.attrs["nCenters"] = int(nAtoms)
        f.attrs["nBasisFunc"] = int(nbf)
        f.attrs["nActOrbs"] = int(nActOrbs)
        f.attrs["nOccOrb"] = int(nOccActOrbs)
        f.attrs["energy"] = float(energy)
        f.attrs["nIntegrals"] = int(nIntegrals)

    print("[hdf5_exporter] Done.")
