# Symmetry and Orbital Representations

This document describes the scientific method, conventions, validation set, and
current limits of QEPlotter's Γ-point symmetry-representation engine.

## Scope

The engine starts from a periodic crystal structure and analyses a complete
atomic symmetry orbit using one of these bases:

- one scalar `s` orbital per selected atom;
- three real `p` orbitals (`px`, `py`, `pz`) per selected atom;
- five real `d` orbitals (`dz²`, `dx²-y²`, `dxy`, `dxz`, `dyz`) per atom;
- three Cartesian atomic displacements (`ux`, `uy`, `uz`) per atom.

It calculates the reducible characters, irreducible decomposition,
symmetry-adapted linear combinations (SALCs), and the orbital sets that share an
irrep and are therefore allowed to interact by symmetry.

## Scientific method and conventions

### 1. Primitive standardisation

The uploaded lattice, fractional positions, and atomic numbers are passed to:

```python
spglib.standardize_cell(
    cell,
    to_primitive=True,
    no_idealize=False,
    symprec=symprec,
)
```

Symmetry is redetected on the resulting primitive cell using the tolerance shown
in the GUI. Complete Wyckoff orbits are used instead of arbitrary atom subsets,
so the selected basis remains closed under every operation.

### 2. Γ-point factor group

At Γ, the translational Bloch phase is one. Operations with the same rotational
part can therefore be represented by one point-group operation.

ASE stores real-space lattice vectors as the rows of $A$. Fractional
rotations are converted to Cartesian form using:

$$
R_\mathrm{cart}
=A^\mathrm{T}R_\mathrm{frac}A^{-\mathrm{T}}.
$$

A polar decomposition removes insignificant floating-point loss of
orthogonality.

### 3. Orbital representation

A symmetry operation permutes the selected atoms and transforms the local
orbital block. QEPlotter uses:

$$
D_s(g)=1,\qquad
D_p(g)=R_\mathrm{cart}(g).
$$

Atomic displacements use the same polar-vector block as the `p` basis.

The five real `d` functions are represented as normalised symmetric traceless
tensors. Each tensor transforms as:

$$
Q'=RQR^\mathrm{T}.
$$

The component order is:

```text
dz², dx²-y², dxy, dxz, dyz
```

### 4. Reducible characters

For a selected orbital basis, the character is:

$$
\chi_\Gamma(g)=\operatorname{Tr}[D_\Gamma(g)].
$$

An atom moved to another site contributes zero to the trace. An atom left on
its site contributes the trace of its orbital transformation block. Characters
are averaged only within exact conjugacy classes.

### 5. Irreducible characters

QEPlotter does not use a partial hard-coded character table. It constructs:

1. the finite point-group multiplication table;
2. the conjugacy classes;
3. the regular representation;
4. the central class-sum operators.

The commuting class-sum operators are resolved numerically. Their eigenspaces
give the complex irreducible characters. Complex-conjugate pairs are then
combined into real crystallographic irreps.

The finite-group engine is regression-tested against all 32 crystallographic
point groups.

### 6. Irreducible decomposition

The multiplicity of irrep $\alpha$ is:

$$
n_\alpha
=\frac{1}{|G|}
\sum_C |C|\,\chi_\Gamma(C)\chi_\alpha(C)^*.
$$

A result is accepted only if every multiplicity is an integer and:

$$
\sum_\alpha n_\alpha d_\alpha
=\dim(\Gamma).
$$

This dimension check prevents a missing symmetry class from producing a
plausible-looking but incomplete result.

### 7. Symmetry-adapted linear combinations

SALCs are generated using the projection operator:

$$
P_\alpha
=\frac{d_\alpha}{|G|}
\sum_g\chi_\alpha(g)^*D_\Gamma(g).
$$

The projected subspace is orthonormalised before its coefficients are shown.
Different bases within a degenerate irrep are equivalent rotations of the same
symmetry subspace.

Large projection matrices are generated only when requested, keeping the
character and decomposition results interactive.

### 8. Orbital compatibility

Two orbital sets are reported as symmetry-compatible only when their
decompositions share at least one irrep.

This is a necessary selection rule, not proof of strong hybridisation. Energy
alignment and spatial overlap are not inferred from symmetry alone.

## Irrep labels

Conventional Mulliken labels are emitted where their mapping is unambiguous,
including the cubic `T`, `T_d`, `O`, and `O_h` families.

For other groups, the GUI keeps deterministic `Γ` labels together with
dimensional A/E/T-like hints. It does not present an orientation-dependent
Mulliken subscript as certain.

## Automated validation

| Structure | Detected space group | Point group | Validation |
|---|---:|---:|---|
| tetrahedral CH₄ test crystal | `P-43m` (#215) | `-43m` (`T_d`) | H `1s`: `A1 + T2`, SALC coefficients |
| diamond Si | `Fd-3m` (#227) | `m-3m` (`O_h`) | all primitive orbits, `s/p/d` dimensions |
| rocksalt NaCl | `Fm-3m` (#225) | `m-3m` (`O_h`) | both ionic sites, `s/p/d` dimensions |
| hcp Mg | `P6_3/mmc` (#194) | `6/mmm` (`D_6h`) | nonsymmorphic primitive-orbit action |
| wurtzite ZnO | `P6_3mc` (#186) | `6mm` (`C_6v`) | polar, non-centrosymmetric structure |
| 2H monolayer MoS₂ | `P-6m2` (#187) | `-6m2` (`D_3h`) | 2D slab with vacuum and two site orbits |

The tests verify:

- the full operation count and presence of every conjugacy class;
- integer irrep multiplicities;
- reconstruction of the original basis dimension;
- known CH₄/Td decomposition and its fully symmetric SALC;
- orbital compatibility between ligand and central-site channels;
- the finite-group character engine for all 32 crystallographic point groups.

## Current scientific boundary

The current page analyses single-valued, spinless, structure-derived
representations at Γ. It does not yet calculate:

- spinor or double-group representations;
- magnetic space-group representations;
- non-Γ little-group irreps;
- symmetry labels of individual Quantum ESPRESSO bands;
- wavefunction-derived orbital characters;
- `f`-orbital representations;
- disordered or partially occupied crystallographic sites.

Individual band labels require wavefunction-level symmetry data and cannot be
inferred reliably from geometry alone.
