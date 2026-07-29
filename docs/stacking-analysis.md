# Bilayer stacking analysis

QEPlotter first verifies that the uploaded structure contains two
vacuum-separated slabs. It does not run stacking classification for a
monolayer or for a bulk structure without a clear bilayer separation.

## Supported six-site convention

For commensurate trigonal-prismatic MX2-like layers, including ordinary and
Janus TMD monolayers, labels are ordered from the lower layer to the upper
layer:

| Orientation family | Registry | Geometric condition |
|---|---|---|
| R-type, parallel | AA | core over core and interface over interface |
| R-type, parallel | AB | upper interface over lower core |
| R-type, parallel | BA | upper core over lower interface |
| H-type, antiparallel | AA′ | both core/interface pairs cross-aligned |
| H-type, antiparallel | AB′ | core over core; interface sites staggered |
| H-type, antiparallel | A′B | interface over interface; core sites staggered |

Keeping the lower-to-upper order is important. AB and BA can be related by
layer exchange in an ideal symmetric homobilayer, but they are not generally
equivalent in heterobilayers, Janus structures, polar domains, or moiré
systems.

The code identifies the core and two surface atomic planes geometrically. It
does not guess the metal by choosing the least frequent element. A canonical
label is accepted only when at least 80% of both compared site sets coincide
within the in-plane tolerance. This bidirectional test prevents one accidental
coincidence in a supercell from defining the whole structure.

## Janus and heterobilayer output

A registry label alone is ambiguous for a Janus material because its two
surfaces have different species. QEPlotter therefore reports:

- lower and upper layer formula separately;
- homobilayer, heterobilayer, Janus homobilayer, or Janus heterobilayer;
- facing interface termination, such as `Se | S`;
- the ordered registry and orientation family.

For example:

```text
AA′ · H-type (antiparallel)
lower → upper: MoSSe / WSSe
facing interface: Se | S
```

## Safe fallback

The six labels are not universal names for every layered crystal. If the
structure is twisted, incommensurate, translated away from a high-symmetry
site, octahedral, buckled, or otherwise outside the supported
trigonal-prismatic geometry, QEPlotter returns:

```text
General registry
fractional core shift: (u, v)
```

It never chooses the nearest AA label merely to produce a result. Multilayer
sequences such as ABA/ABC and polytype names such as 2H/3R are different
descriptions and are not assigned to a verified two-layer structure by this
classifier.

## Naming references

- *Interlayer Registry Index of Layered Transition Metal
  Dichalcogenides*, Journal of Physical Chemistry Letters (parallel
  AA/AB/BA and antiparallel AA′/AB1/AB2 configurations):
  <https://doi.org/10.1021/acs.jpclett.1c04202>
- *Stacking stability of MoS2 bilayer: An ab initio study*
  (AA, AB, AA′, AB′, and A′B convention):
  <https://doi.org/10.1088/1674-1056/23/10/106801>
- *Relaxation effects in transition metal dichalcogenide
  bilayer heterostructures* (R/H local registry notation):
  <https://doi.org/10.1038/s41699-024-00477-6>

Different publications sometimes exchange AB1/AB2 with AB′/A′B or reverse
top/bottom order. QEPlotter avoids this ambiguity by publishing its geometric
condition and lower-to-upper convention beside the label.
