# Scaffolded fill-in-the-code worksheet

Complete only the marked expressions; hints are below each task.

```python
d, xp, xn = depletion_edges(Na, Nd, Vbi)
assert np.isclose(Na*___, Nd*___)       # Hint: charge neutrality.
peak = np.max(np.abs(___))              # Hint: electric-field array.
ratio = xp/xn
```

Explain why `delta_num` is not a Debye or diffusion length. Then add a loop
over three doping ratios and tabulate `xp/xn`.
