# buckling_braced_L120.jl — eigenbuckling of the welded pair, L = 120 in, pinned warping-free ends, 3 in welds at 18 in,
# with frame bracing at z = 6, 54, 102 in (AISI Lx = 120 in: deflection in Y unbraced; Ly = 48 in and Lt = 48 in:
# X translation and twist restrained at the braces). Unconstrained shell model (global mode is mode 1 at this length).
include(joinpath(@__DIR__, "buckling_built_up.jl"))
L = 120.0; braces = [6.0, 54.0, 102.0]
println("Braced, L = 120: braces at z = $braces restrain u_x of every node (Ly = Lt = 48 in), u_y free (Lx = 120 in)")
a = built_up_buckling(; L, mode = :all, nev = 4, braces, brace_components = [1])
# (the rigid-section model is not used here: brace Dirichlet conditions would fall on affine master dofs, which Ferrite
#  does not support; at this length the unconstrained shell's first mode is the global mode)
println("\nUnbraced reference (from buckling_results_L120.csv): unconstrained mode 1 = 49.2 kips")
# beam-theory estimate: P_ey over Lx = 120 (deflection in Y), P_t over Lt = 48 with G J_eff (C_w = 0), coupled FT formula
X1, Y1 = centerline(shape); sp = section_properties(X1, Y1, shape.t)
A2 = 2sp.A; Ix2 = 2sp.Ixx; Iy2 = 2 * (sp.Iyy + sp.A * (shape.B / 2)^2); J_eff = 0.7678; x_o = -2.355
r_o2 = (Ix2 + Iy2) / A2 + x_o^2; β = 1 - x_o^2 / r_o2
Pey = π^2 * E * Ix2 / 120^2; Pt = G * J_eff / r_o2; Pex48 = π^2 * E * Iy2 / 48^2
timo = ((Pey + Pt) - sqrt((Pey + Pt)^2 - 4β * Pey * Pt)) / (2β)
@printf("\nTimoshenko: P_ey(Lx = 120) = %.1f kips, P_t(Lt = 48, G J_eff, C_w = 0) = %.1f kips, P_ex(Ly = 48) = %.1f kips → coupled FT = %.1f kips\n", Pey/1e3, Pt/1e3, Pex48/1e3, timo/1e3)
open(joinpath(@__DIR__, "buckling_results_L120_braced.csv"), "w") do io
    println(io, "analysis,mode,P_cr_lbf,P_cr_kips,sigma_cr_ksi,rigid_section_participation,U_mid,V_mid,theta_mid,twist_to_translation")
    for (name, r) in (("all", a),), (k, P) in enumerate(r.P_cr)
        c = r.info[k]; println(io, join([name, k, P, P / 1000, P / r.A_total / 1000, c.participation, c.U, c.V, c.θ, c.twist_to_translation], ","))
    end
    println(io, "timoshenko,FT,$timo,$(timo/1000),,,,,,")
end
write_mode(joinpath(@__DIR__, "mode_all_1_L120_braced.csv"), a, 1)
println("Wrote buckling_results_L120_braced.csv, mode_all_1_L120_braced.csv")
