# buckling_braced_L120_perforated.jl
#
# Elastic eigenbuckling (P_cre) of the PERFORATED two-C welded upright, L = 120 in, pinned warping-free ends,
# 3 in welds at 18 in spacing (weld centers z = 6, 24, …, 114 in), with frame bracing at z = 6, 54, 102 in
# (AISI Lx = 120 in: deflection in Y unbraced; Ly = Lt = 48 in: X translation and twist restrained at the braces),
# i.e. the same problem as buckling_braced_L120.jl but on the Gmsh mixed quad/triangle mesh of perforated_J_gmsh.jl
# (teardrop web holes, square flange holes, same pattern as the torsion J study).
#
# Model: Ferrite.jl + QuadShellFiniteElement.jl (quads) + TriShellFiniteElement.jl (triangles at the hole
# boundaries), Hughes–Brezzi drilling. Reference load P_ref = 1 kip applied as tributary nodal forces at both
# (unperforated) end sections; the prebuckling membrane stresses of the unconstrained shell (non-uniform around
# the holes) feed the mixed geometric stiffness; K φ = P (−K_g) φ on Cᵀ K C by shift-invert Arnoldi.
# Braces: Dirichlet u_x = 0 on every node of the mesh cross line at z = 6, 54, 102 (the Gmsh strip is fragmented
# by these lines so nodes exist exactly there). The weld at z = 6 overlaps the brace at z = 6: its rigid-body master
# is taken off the brace line and the slaves on the brace line are tied in u_y, u_z, θx, θy, θz only (u_x = 0).
#
# Run:  julia --project=. buckling_braced_L120_perforated.jl            (perforated, writes *_perforated.csv)
#       julia --project=. buckling_braced_L120_perforated.jl gross      (same pipeline without holes: check against
#                                                                        buckling_braced_L120.jl, 51.2 kips)

include(joinpath(@__DIR__, "perforated_J_gmsh.jl"))
using ArnoldiMethod, LinearMaps, SparseArrays

const P_ref = 1000.0     # lbf (1 kip) total reference compression
const MAP_15_TO_18 = [1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17]

# rigid tie of the components in `comps` (1:3 translations, 4:6 rotations) of `slaves` to `master`
function add_rigid_tie_partial!(ch, nd, grid, master, slaves, comps)
    xm = grid.nodes[master].x; dm = nd[master]
    for s in slaves
        s == master && continue
        r = grid.nodes[s].x - xm; ds = nd[s]
        1 in comps && add!(ch, AffineConstraint(ds[1], [dm[1] => 1.0, dm[5] => r[3], dm[6] => -r[2]], 0.0))
        2 in comps && add!(ch, AffineConstraint(ds[2], [dm[2] => 1.0, dm[6] => r[1], dm[4] => -r[3]], 0.0))
        3 in comps && add!(ch, AffineConstraint(ds[3], [dm[3] => 1.0, dm[4] => r[2], dm[5] => -r[1]], 0.0))
        for k in 4:6
            k in comps && add!(ch, AffineConstraint(ds[k], [dm[k] => 1.0], 0.0))
        end
    end
end

symmetric_weld_locations(L, wl, ws) = (n = floor(Int, (L - wl) / ws) + 1; collect(range(L / 2 - (n - 1) * ws / 2, L / 2 + (n - 1) * ws / 2, n)))

"""
Mixed quad/triangle geometric stiffness from the membrane stresses of the prebuckling displacement `u`
(plane stress, element local frame). Returns K_g and the area-weighted mean global axial stress σ_zz.
"""
function assemble_mixed_Kg!(Kg, dh, u, t)
    Dσ = E / (1 - ν^2) * [1.0 ν 0.0; ν 1.0 0.0; 0.0 0.0 (1 - ν) / 2]
    cvq = CellValues(QuadratureRule{RefQuadrilateral}(2), IP4(), IP4())
    cvt = CellValues(QuadratureRule{RefTriangle}(1), TS.IP3(), TS.IP3())
    assembler = start_assemble(Kg)
    σz_sum = 0.0; area = 0.0
    for sdh in dh.subdofhandlers, cell in CellIterator(sdh)
        x = getcoordinates(cell); cd = celldofs(cell); k = length(x)
        if k == 4
            T = Q.calculation_rotation_matrix(x); xl = Q.global_nodal_coords_to_planar_coords(x, T); reinit!(cvq, xl); cv = cvq
        else
            T = TS.calculation_rotation_matrix(x); xl = TS.global_nodal_coords_to_planar_coords(x, T); reinit!(cvt, xl); cv = cvt
        end
        ul = zeros(2k)
        for i in 1:k
            ulo = T' * [u[cd[3(i-1)+m]] for m in 1:3]                 # :u dofs come first in the cell dof vector
            ul[2i-1] = ulo[1]; ul[2i] = ulo[2]
        end
        nq = getnquadpoints(cv); sx = zeros(nq); sy = zeros(nq); sxy = zeros(nq); Bm = zeros(3, 2k)
        for q in 1:nq
            for i in 1:k
                dN = shape_gradient(cv, q, i)
                Bm[1, 2i-1] = dN[1]; Bm[2, 2i] = dN[2]; Bm[3, 2i-1] = dN[2]; Bm[3, 2i] = dN[1]
            end
            σ = Dσ * (Bm * ul); sx[q] = σ[1]; sy[q] = σ[2]; sxy[q] = σ[3]
            dV = getdetJdV(cv, q); area += dV
            e1z = T[3, 1]; e2z = T[3, 2]                                     # global σ_zz from the local tensor
            σz_sum += (σ[1] * e1z^2 + σ[2] * e2z^2 + 2σ[3] * e1z * e2z) * dV
        end
        if k == 4
            kg = Q.calculate_element_geometric_stiffness_matrix(cv, sx .* t, sy .* t, sxy .* t)
            kg24 = zeros(24, 24); kg24[Q.MAP_20_TO_24, Q.MAP_20_TO_24] = kg
            Te = Q.rotation_matrix_for_element_stiffness_drilling(T)
            kgg = (Te * kg24 * Te')[Q.FIELD_ORDER_24, Q.FIELD_ORDER_24]
        else
            kg = TS.calculate_element_geometric_stiffness_matrix(cv, sx .* t, sy .* t, sxy .* t, T)
            kg18 = zeros(18, 18); kg18[MAP_15_TO_18, MAP_15_TO_18] = kg
            Te = TS.rotation_matrix_for_element_stiffness_drilling(T)
            kgg = (Te * kg18 * Te')[IND18, IND18]
        end
        assemble!(assembler, cd, kgg)
    end
    return Kg, σz_sum / area
end

"""
    braced_buckling_gmsh(; L, perforated, braces, weld_length, weld_spacing, mesh_size, nev)

Eigenbuckling of the welded pair on the Gmsh mesh with pinned warping-free ends and X-translation braces at `braces`.
"""
function braced_buckling_gmsh(; L = 120.0, perforated = true, braces = [6.0, 54.0, 102.0], weld_length = 3.0, weld_spacing = 18.0,
                                mesh_size = 0.2, n_flat = 4, n_corner = 4, nev = 4, verbose = true)
    t = shape.t; B = shape.B; D = shape.D; R = shape.R
    weld_locations = symmetric_weld_locations(L, weld_length, weld_spacing)
    t0 = time()
    grid2d, fold, S = strip_mesh(L; perforated, mesh_size, n_flat, n_corner, zlines = braces)
    grid = folded_grid(grid2d, fold, 2)
    dh, nd = mixed_dofhandler(grid)
    P = [n.x for n in grid.nodes]
    hasd(i) = haskey(nd, i)
    near(x, x0, tol) = abs(x - x0) <= tol + 1e-9
    part(i) = (P[i][1] - (P[i][1] > B + 0.5 ? B : 0.0), P[i][2]); shapeof(i) = P[i][1] > B + 0.5 ? 2 : 1
    end0 = Set(i for i in eachindex(P) if hasd(i) && near(P[i][3], 0.0, 1e-6)); endL = Set(i for i in eachindex(P) if hasd(i) && near(P[i][3], L, 1e-6))
    onbrace(i) = any(near(P[i][3], zb, 1e-6) for zb in braces)
    verbose && @printf("  mesh: %d quads + %d tris, %d nodes, %d dofs (%.0f s)\n", length(getcellset(grid, "quad")), length(getcellset(grid, "tri")), getnnodes(grid), ndofs(dh), time() - t0)

    # ---- constraints
    ch = ConstraintHandler(dh)
    add!(ch, Dirichlet(:u, union(end0, endL), (x, t) -> [0.0, 0.0], [1, 2]))                       # pinned, warping free
    datum = argmin([hasd(i) ? (part(i)[1])^2 + (part(i)[2] - D / 2)^2 + (P[i][3] - L / 2)^2 + (shapeof(i) == 2 ? 0.0 : 1e3) : Inf for i in eachindex(P)])
    add!(ch, Dirichlet(:u, Set([datum]), (x, t) -> [0.0], [3]))                                    # axial datum
    n_brace_nodes = 0
    for zb in braces
        bn = Set(i for i in eachindex(P) if hasd(i) && near(P[i][3], zb, 1e-6))
        n_brace_nodes += length(bn)
        add!(ch, Dirichlet(:u, bn, (x, t) -> [0.0], [1]))                                          # brace: u_x = 0 (Ly, Lt)
    end
    claimed = falses(length(P)); n_ties = 0
    for zw in weld_locations, yw in (D, 0.0)
        slaves = [i for i in eachindex(P) if hasd(i) && !claimed[i] && abs(P[i][3] - zw) <= weld_length / 2 + 1e-9 && P[i][3] > 1e-6 && P[i][3] < L - 1e-6 &&
                  near(part(i)[2], yw, R) && ((shapeof(i) == 1 && near(part(i)[1], B, R)) || (shapeof(i) == 2 && near(part(i)[1], 0.0, R)))]
        isempty(slaves) && continue
        claimed[slaves] .= true
        master = slaves[argmax([minimum(abs(P[i][3] - zb) for zb in vcat(braces, -1e3)) for i in slaves])]   # master off the brace lines
        add_rigid_tie_partial!(ch, nd, grid, master, [i for i in slaves if !onbrace(i)], 1:6)
        add_rigid_tie_partial!(ch, nd, grid, master, [i for i in slaves if onbrace(i)], 2:6)
        n_ties += length(slaves) - 1
    end
    close!(ch); update!(ch, 0.0)

    # ---- stiffness, reference load (tributary nodal forces on the end edges), prebuckling solve
    K = allocate_matrix(dh, ch); K = assemble_mixed_Ke!(K, dh, t)
    f = zeros(ndofs(dh)); A_end = Dict(0 => 0.0, 1 => 0.0); trib = Dict{Int,Float64}()
    for c in grid.cells, k in 1:length(c.nodes)
        a, b = c.nodes[k], c.nodes[mod1(k + 1, length(c.nodes))]
        for (e, zs) in ((0, 0.0), (1, L))
            if near(P[a][3], zs, 1e-6) && near(P[b][3], zs, 1e-6)
                len = hypot(P[a][1] - P[b][1], P[a][2] - P[b][2])
                trib[a] = get(trib, a, 0.0) + len / 2; trib[b] = get(trib, b, 0.0) + len / 2; A_end[e] += len * t
            end
        end
    end
    A_gross = 2 * S * t
    for (n, l) in trib
        f[nd[n][3]] += (near(P[n][3], 0.0, 1e-6) ? 1.0 : -1.0) * P_ref / A_gross * l * t
    end
    K0 = copy(K); apply!(K, f, ch); u = K \ f; apply!(u, ch)
    Kg = allocate_matrix(dh); Kg, σz_mean = assemble_mixed_Kg!(Kg, dh, u, t)
    verbose && @printf("  %d brace nodes, %d weld ties, %d prescribed dofs; A_gross = %.4f in² (end edge areas %.4f / %.4f), σ_ref = %.2f psi, mean prebuckling σ_z = %.2f psi (%.0f s)\n",
                       n_brace_nodes, n_ties, length(ch.prescribed_dofs), A_gross, A_end[0], A_end[1], P_ref / A_gross, -σz_mean, time() - t0)

    # ---- condensed eigenproblem
    C, _ = Ferrite.create_constraint_matrix(ch)
    Kc = sparse(C' * K0 * C); Kgc = sparse(C' * (-Kg) * C)
    F = lu(Kc)
    op = LinearMap{Float64}(x -> F \ (Kgc * x), size(Kc, 1); ismutating = false)
    decomp, hist = partialschur(op; nev = nev + 4, tol = 1e-8, which = isdefined(ArnoldiMethod, :LR) ? ArnoldiMethod.LR() : :LR)
    μ, Ψ = partialeigen(decomp)
    keep = findall(m -> real(m) > 1e-12, μ)
    order = sortperm(real.(μ[keep]); rev = true)[1:min(nev, length(keep))]
    P_cr = [P_ref / real(μ[keep[k]]) for k in order]
    modes = [C * real.(Ψ[:, keep[k]]) for k in order]

    # ---- classification: rigid-section share of the in-plane displacements (1 in bins along z), mid-length motion
    nodes_d = [i for i in eachindex(P) if hasd(i)]
    bins = [[i for i in nodes_d if z0 < P[i][3] <= z0 + 1.0] for z0 in 0.0:1.0:L-1.0]
    mid = [i for i in nodes_d if abs(P[i][3] - L / 2) <= 0.5]
    rigid_stats(a, ids) = begin
        xs = [P[i][1] for i in ids]; ys = [P[i][2] for i in ids]
        ux = [a[nd[i][1]] for i in ids]; uy = [a[nd[i][2]] for i in ids]
        U = mean(ux); V = mean(uy); θ = ls_rotation(xs, ys, ux, uy)
        uxr = U .- θ .* (ys .- mean(ys)); uyr = V .+ θ .* (xs .- mean(xs))
        (sum(uxr .^ 2 .+ uyr .^ 2), sum(ux .^ 2 .+ uy .^ 2), U, V, θ, maximum(hypot.(xs .- mean(xs), ys .- mean(ys))))
    end
    info = map(modes) do a
        rig = 0.0; tot = 0.0
        for b in bins
            length(b) < 6 && continue
            r_, t_, = rigid_stats(a, b); rig += r_; tot += t_
        end
        _, _, U, V, θ, rmax = rigid_stats(a, mid)
        (; participation = rig / tot, U, V, θ, twist_to_translation = abs(θ) * rmax / max(hypot(U, V), 1e-30))
    end
    if verbose
        for (k, Pk) in enumerate(P_cr)
            c = info[k]
            @printf("  mode %d: P_cr = %10.1f lbf = %8.2f kips, σ_cr(gross) = %7.2f ksi;  rigid-section participation %.3f;  mid-length U = %+.3f V = %+.3f θ·r_max/|U,V| = %.2f\n",
                    k, Pk, Pk / 1000, Pk / A_gross / 1000, c.participation, c.U / maximum(abs, [c.U, c.V, 1e-30]), c.V / maximum(abs, [c.U, c.V, 1e-30]), c.twist_to_translation)
        end
        @printf("  (%.0f s total)\n", time() - t0)
    end
    return (; P_cr, modes, info, grid, nd, dh, A_gross, weld_locations, braces, L, perforated,
            nquad = length(getcellset(grid, "quad")), ntri = length(getcellset(grid, "tri")), ndofs = ndofs(dh))
end

function write_mode_mesh(prefix, r, k)
    a = r.modes[k]; g = r.grid; nd = r.nd
    open(prefix * "_nodes.csv", "w") do io
        println(io, "x,y,z,ux,uy,uz")
        for (n, node) in enumerate(g.nodes)
            p = node.x
            if haskey(nd, n); d = nd[n]; println(io, join([p[1], p[2], p[3], a[d[1]], a[d[2]], a[d[3]]], ","))
            else; println(io, join([p[1], p[2], p[3], 0.0, 0.0, 0.0], ",")); end
        end
    end
    open(prefix * "_cells.csv", "w") do io
        for c in g.cells; println(io, join(c.nodes, ",")); end
    end
    open(prefix * "_summary.csv", "w") do io
        println(io, "L_in,P_cre_kips,mode,quads,tris,ndofs,A_gross_in2,braces,welds")
        println(io, join([r.L, r.P_cr[k] / 1000, k, r.nquad, r.ntri, r.ndofs, r.A_gross, "\"" * join(r.braces, ";") * "\"", "\"" * join(r.weld_locations, ";") * "\""], ","))
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    perforated = !(length(ARGS) >= 1 && ARGS[1] == "gross")
    tag = perforated ? "L120_braced_perforated" : "L120_braced_gmsh_gross"
    println("Two-C welded upright, L = 120 in, pinned warping-free, 3 in welds at 18 in, braces (u_x = 0) at z = 6, 54, 102 in; ",
            perforated ? "PERFORATED (teardrop web holes, square flange holes), " : "no holes (pipeline check), ", "Gmsh mixed mesh\n")
    r = braced_buckling_gmsh(; perforated, nev = 4)
    println("\nReference (structured quad mesh, no holes, buckling_results_L120_braced.csv): mode 1 = 51.2 kips (global, Y-flexure over Lx = 120 in)")
    open(joinpath(@__DIR__, "buckling_results_$(tag).csv"), "w") do io
        println(io, "analysis,mode,P_cr_lbf,P_cr_kips,sigma_cr_gross_ksi,rigid_section_participation,U_mid,V_mid,theta_mid,twist_to_translation")
        for (k, Pk) in enumerate(r.P_cr)
            c = r.info[k]; println(io, join(["all", k, Pk, Pk / 1000, Pk / r.A_gross / 1000, c.participation, c.U, c.V, c.θ, c.twist_to_translation], ","))
        end
    end
    write_mode_mesh(joinpath(@__DIR__, "mode_$(tag)"), r, 1)
    println("Wrote buckling_results_$(tag).csv, mode_$(tag)_{nodes,cells,summary}.csv")
end
