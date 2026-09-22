# Built-up upright torsion and buckling studies (Ferrite.jl)

Shell finite element studies of a cold-formed steel upright made of two 3 × 3 × 0.074 in lipped C sections
welded back-to-lip with intermittent 3 in welds, built with [Ferrite.jl](https://ferrite-fem.github.io/Ferrite.jl/)
and [QuadShellFiniteElement.jl](https://github.com/runtosolve/QuadShellFiniteElement.jl).

- `calculate_J_single_member/` — Saint-Venant torsion constant J of one C by the Moen (2008) §4.2.7.3.2.3 static twist
  method, checked against an exact 2D Saint-Venant solution; documents the drilling-dof finding that led to the
  Hughes–Brezzi option in QuadShellFiniteElement.jl and TriShellFiniteElement.jl.
- `calculate_J_built_up_member/` — effective J of the welded pair, weld spacing study vs. the Tlumak equation, and the
  global flexural-torsional eigenbuckling of the 44 in member.
- `docs/` — GitHub Pages site: interactive WGLMakie buckling mode shape (`index.html`) and summary figures.

**Hosted mode shape:** https://runtosolve.github.io/BuiltUpUprightTorsion/

Each study folder has its own README with the model description, results and file list. Run the scripts from their
folder with `julia --project=.` (the Project/Manifest reference packages developed under `~/.julia/dev`).
