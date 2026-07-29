# Irreducible k-grid engine

QEPlotter keeps uniform integration meshes separate from conventional
high-symmetry band paths. The **K-grid & IBZ** workspace accepts a periodic
structure and a uniform `nk1 × nk2 × nk3` mesh, then produces the
symmetry-inequivalent representatives, multiplicities, normalized weights,
full-grid mapping, a first-Brillouin-zone view, and Quantum ESPRESSO cards.

## Method

1. The uploaded lattice, fractional atomic positions, and atomic numbers are
   passed to spglib to determine the spatial symmetry operations.
2. Each direct-space rotation `W` is converted to its reciprocal action
   `R = (W⁻¹)ᵀ`.
3. The uniform mesh and its half-step shift are represented by integer grid
   addresses on a common exact denominator.
4. Operations that do not map the selected shifted mesh onto itself are
   excluded. The interface reports this explicitly; it never silently applies
   an incompatible crystal rotation.
5. QEPlotter applies every compatible reciprocal rotation, and optionally
   time reversal `k → −k`, to build exact symmetry orbits.
6. The smallest full-grid index is selected deterministically as each orbit's
   representative. Its multiplicity is the orbit size and its normalized
   integration weight is `multiplicity / (nk1 nk2 nk3)`.

The reduction code does **not** call spglib's reciprocal-mesh reducer. spglib is
the symmetry source only; QEPlotter owns grid generation, compatibility,
mapping, representative selection, and weighting.

## Coordinates and Quantum ESPRESSO

Fractional points use the reciprocal basis of the uploaded cell:

`B = 2π(A⁻¹)ᵀ`, with row-vector convention `k_cart = k_frac · B`.

The page exports:

- `K_POINTS automatic`, recommended for ordinary SCF/NSCF calculations because
  Quantum ESPRESSO performs its own symmetry handling; and
- explicit weighted `K_POINTS crystal`, with either normalized weights or
  integer multiplicities.

Normalized weights sum to one. Multiplicity weights sum to the number of
points in the complete grid; the two choices differ only by a common factor.

## Time reversal and magnetic calculations

Time reversal makes `k` and `−k` equivalent. It is normally valid for
non-magnetic calculations. Spin-orbit coupling alone does not necessarily
break time reversal, but magnetic order can. For magnetic or non-collinear
states, the user must select the option according to the physical calculation.
When a QE input contains magnetic indicators, the page displays a warning
instead of silently assuming that time reversal is valid.

The engine is deliberately structure-based: it does not infer a magnetic space
group from a geometry-only CIF/POSCAR. If time reversal is disabled, exported
QE cards include a reminder to set `noinv=.true.` in `&SYSTEM`; the `K_POINTS`
card alone cannot encode that QE setting.

## Validation invariants

Every generated result is checked before it reaches the interface:

- every full-grid point maps to exactly one irreducible representative;
- the sum of multiplicities equals the complete mesh size;
- normalized weights sum to one;
- only symmetry operations that preserve the chosen mesh and shift are used;
- reciprocal-space equivalence is computed with integer arithmetic rather
  than floating-point coordinate matching.
