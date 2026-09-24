# Built-up upright torsion and buckling studies (Ferrite.jl)

Shell finite element studies of a cold-formed steel upright made of two 3 × 3 × 0.074 in lipped C sections
welded back-to-lip with intermittent 3 in welds, built with [Ferrite.jl](https://ferrite-fem.github.io/Ferrite.jl/)
and [QuadShellFiniteElement.jl](https://github.com/runtosolve/QuadShellFiniteElement.jl).

- `rmi_education_summit_2026/calculate_J_single_member/` — Saint-Venant torsion constant J of one C by the Moen (2008) §4.2.7.3.2.3 static twist
  method, checked against an exact 2D Saint-Venant solution; documents the drilling-dof finding that led to the
  Hughes–Brezzi option in QuadShellFiniteElement.jl and TriShellFiniteElement.jl.
- `rmi_education_summit_2026/calculate_J_built_up_member/` — effective J of the welded pair, weld spacing study vs. the Tlumak equation, and the
  global flexural-torsional eigenbuckling of the 44 in member.
- `docs/` — GitHub Pages site: browser-interactive 3D buckling mode shape (`index.html`, plotly.js; rotate/zoom in the browser),
  the static WGLMakie export (`mode_global_1_wglmakie_static.html`, snapshot only — WGLMakie's camera needs a live Julia
  session, see `show_mode_wglmakie.jl`), and summary figures.

**Hosted pages**

- Global flexural-torsional buckling mode, L = 44 in (rigid-section model, P_cre = 298.9 kips): https://runtosolve.github.io/UprightLTB/
- Lowest global buckling mode, L = 120 in (unconstrained shell, P_cre = 49.2 kips): https://runtosolve.github.io/UprightLTB/L120.html
- L = 120 in with frame bracing at z = 6, 54, 102 in (L_x = 120, L_y = L_t = 48 in; P_cre = 51.2 kips): https://runtosolve.github.io/UprightLTB/L120_braced.html
- Static twist study for J_eff (L = 111 in, 3 in welds at 18 in) is hosted separately at
  https://runtosolve.github.io/Upright_torsion_J_study/ (repository runtosolve/Upright_torsion_J_study).

Each study folder has its own README with the model description, results and file list. Run the scripts from their
folder with `julia --project=.` (the Project/Manifest reference packages developed under `~/.julia/dev`).
