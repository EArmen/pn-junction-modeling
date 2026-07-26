# Domain and mesh guidance

For total depletion width `d`, use `xp=Nd*d/(Na+Nd)` and
`xn=Na*d/(Na+Nd)`, so `Na*xp=Nd*xn`. The minimum depletion interval is
`[-xp,xn]`; do not centre a strongly asymmetric junction in `[-d/2,d/2]`.
Place each numerical boundary in its quasi-neutral bulk. A practical starting
point is five majority-carrier Debye lengths beyond each depletion edge.

Run `notebooks/domain_mesh_sensitivity.ipynb`. Enlarge both buffers and refine
the mesh until (i) relative peak-field change and (ii) the field-profile RMS
change on the common interval are each below 1%. The notebook calculates both
quantities and prints PASS or REFINE.
