# make_mode_html_mesh.jl — browser-interactive 3D buckling mode page (plotly.js) from the Gmsh mixed mesh output of
# buckling_braced_L120_perforated.jl (mode_<tag>_{nodes,cells,summary}.csv). Same content and style as
# make_mode_html.jl: deformed shell coloured by displacement magnitude, welds (red), brace outlines (green), end
# loads (orange cones), mid-length section inset. Only boundary edges (hole outlines, ends, lip tips) are drawn to
# keep the file size down.
# Run:  julia --project=. make_mode_html_mesh.jl [tag = L120_braced_perforated]
using DelimitedFiles, Statistics, Printf
include(joinpath(@__DIR__, "..", "calculate_J_single_member", "calculate_J_single_member.jl"))

tag = length(ARGS) >= 1 ? ARGS[1] : "L120_braced_perforated"
prefix = joinpath(@__DIR__, "mode_" * tag)
nodes = readdlm(prefix * "_nodes.csv", ','; skipstart = 1)
xyz = Float64.(nodes[:, 1:3]); u = Float64.(nodes[:, 4:6])
cells = [parse.(Int, split(strip(l), ',')) for l in readlines(prefix * "_cells.csv")]
sm = readdlm(prefix * "_summary.csv", ','; skipstart = 1)
L = Float64(sm[1, 1]); P = Float64(sm[1, 2]); nq = Int(sm[1, 4]); nt = Int(sm[1, 5])
braces = parse.(Float64, split(strip(String(sm[1, 8])), ';')); centers = parse.(Float64, split(strip(String(sm[1, 9])), ';'))
perforated = occursin("perforated", tag)
weld_length = 3.0; weld_spacing = 18.0; B = 3.0; D = 3.0 - 0.074; R = 0.199; scale_frac = 0.035
umag = vec(sqrt.(sum(u .^ 2; dims = 2))); sc = scale_frac * L / maximum(umag)
def = xyz .+ sc .* u
I = Int[]; J = Int[]; K = Int[]
for c in cells
    push!(I, c[1] - 1); push!(J, c[2] - 1); push!(K, c[3] - 1)
    length(c) == 4 && (push!(I, c[1] - 1); push!(J, c[3] - 1); push!(K, c[4] - 1))
end
edge_count = Dict{Tuple{Int,Int},Int}()
for c in cells, k in 1:length(c)
    a, b = c[k], c[mod1(k + 1, length(c))]; e = (min(a, b), max(a, b)); edge_count[e] = get(edge_count, e, 0) + 1
end
edges = collect(keys(edge_count))
bedges = [e for (e, n) in edge_count if n == 1]
ex = Any[]; ey = Any[]; ez = Any[]
for (a, b) in bedges
    push!(ex, def[a, 1]); push!(ey, def[a, 2]); push!(ez, def[a, 3]); push!(ex, def[b, 1]); push!(ey, def[b, 2]); push!(ez, def[b, 3])
    push!(ex, nothing); push!(ey, nothing); push!(ez, nothing)
end
# welds: tie nodes (C1 corners, C2 web corners) within the weld windows
xp = xyz[:, 1] .- (xyz[:, 1] .> B + 0.5) .* B; s1 = xyz[:, 1] .<= B + 0.5
near(x, x0, tol) = abs(x - x0) <= tol + 1e-9
corner = (s1 .& near.(xp, B, R) .& (near.(xyz[:, 2], D, R) .| near.(xyz[:, 2], 0.0, R))) .| (.!s1 .& near.(xp, 0.0, R) .& (near.(xyz[:, 2], D, R) .| near.(xyz[:, 2], 0.0, R)))
inweld = [any(near(z, c, weld_length / 2) for c in centers) for z in xyz[:, 3]]
w = findall(corner .& inweld .& (xyz[:, 3] .> 0) .& (xyz[:, 3] .< L))
# brace outlines: mesh edges with both nodes on a brace line (undeformed)
bx = Any[]; by = Any[]; bz = Any[]
for zb in braces, (a, b) in edges
    (near(xyz[a, 3], zb, 1e-6) && near(xyz[b, 3], zb, 1e-6)) || continue
    push!(bx, xyz[a, 1]); push!(by, xyz[a, 2]); push!(bz, zb); push!(bx, xyz[b, 1]); push!(by, xyz[b, 2]); push!(bz, zb); push!(bx, nothing); push!(by, nothing); push!(bz, nothing)
end
# loads: cones on a subset of the end nodes (every ~0.45 in along the section), pointing into the member
function end_subset(z0)
    ids = findall(near.(xyz[:, 3], z0, 1e-6))
    keep = Int[]
    for i in sort(ids; by = k -> (xyz[k, 1] > B + 0.5, atan(xyz[k, 2] - D / 2, xp[k] - B / 2)))
        (isempty(keep) || hypot(xyz[i, 1] - xyz[keep[end], 1], xyz[i, 2] - xyz[keep[end], 2]) > 0.45) && push!(keep, i)
    end
    keep
end
ends = vcat(end_subset(0.0), end_subset(L))
cx = [xyz[k, 1] for k in ends]; cy = [xyz[k, 2] for k in ends]; cz = [xyz[k, 3] for k in ends]
cw = [xyz[k, 3] > L / 2 ? -1.0 : 1.0 for k in ends]; shaft = 2.0
ax_ = Any[]; ay_ = Any[]; az_ = Any[]
for k in eachindex(ends)
    push!(ax_, cx[k]); push!(ay_, cy[k]); push!(az_, cz[k]); push!(ax_, cx[k]); push!(ay_, cy[k]); push!(az_, cz[k] - cw[k] * shaft)
    push!(ax_, nothing); push!(ay_, nothing); push!(az_, nothing)
end
# mid-length section: nearest mesh node to each centerline point at z = L/2
X1, Y1 = centerline(shape); sec_traces = String[]
for s in 1:2
    xs0 = X1 .+ (s - 1) * B; ys0 = Y1
    rows = [argmin((xyz[:, 1] .- xs0[i]) .^ 2 .+ (xyz[:, 2] .- ys0[i]) .^ 2 .+ (xyz[:, 3] .- L / 2) .^ 2) for i in eachindex(xs0)]
    push!(sec_traces, "{type:'scatter',x:[$(join(round.(xs0; digits = 4), ','))],y:[$(join(round.(ys0; digits = 4), ','))],mode:'lines',line:{color:'#8a8a8a',width:4},showlegend:false,hoverinfo:'skip',xaxis:'x',yaxis:'y'}")
    push!(sec_traces, "{type:'scatter',x:[$(join(round.(xs0 .+ sc .* u[rows, 1]; digits = 4), ','))],y:[$(join(round.(ys0 .+ sc .* u[rows, 2]; digits = 4), ','))],mode:'lines',line:{color:'#2a78d6',width:4},showlegend:false,hoverinfo:'skip',xaxis:'x',yaxis:'y'}")
end
rx = maximum(def[:, 1]) - minimum(def[:, 1]); ry = maximum(def[:, 2]) - minimum(def[:, 2]); rz = L + 2 * shaft
ar = round.([rx, ry, rz] ./ max(rx, ry, rz) .* 1.9; digits = 3)
js(v) = join((x === nothing ? "null" : string(round(x, digits = 4)) for x in v), ',')
holes_note = perforated ? "PERFORATED (teardrop web holes, square flange holes; Gmsh mesh, $(nq) quads + $(nt) triangles)" : "no holes (Gmsh mesh, $(nq) quads)"
title = "Lowest buckling mode, two-C welded upright, $(holes_note), L = $(Int(L)) in, pinned warping-free, 3 in welds at 18 in:  P<sub>cre</sub> = $(round(P, digits = 1)) kips (unconstrained shell)<br>" *
        "frame bracing (green outlines) at z = $(join(Int.(braces), ", ")) in restrains cross-aisle translation and twist (L<sub>y</sub> = L<sub>t</sub> = 48 in); downaisle flexure unbraced (L<sub>x</sub> = $(Int(L)) in)"

html = """
<!doctype html>
<html><head><meta charset="utf-8"><title>Perforated built-up upright buckling mode</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>html,body{margin:0;height:100%;font-family:Helvetica,Arial,sans-serif;background:#fff}#plot{width:min(100vw,1400px);height:min(90vh,760px);margin:0 auto;overflow:hidden}</style></head>
<body><div id="plot"></div>
<script>
const data = [
 {type:'mesh3d',x:[$(js(def[:, 1]))],y:[$(js(def[:, 2]))],z:[$(js(def[:, 3]))],i:[$(join(I, ','))],j:[$(join(J, ','))],k:[$(join(K, ','))],
  intensity:[$(js(umag ./ maximum(umag)))],colorscale:'Viridis',cmin:0,cmax:1,flatshading:true,
  lighting:{ambient:0.9,diffuse:0.2,specular:0.0},colorbar:{title:{text:'normalized<br>displacement'},len:0.5,x:0.5},
  hoverinfo:'skip',name:'deformed shell',scene:'scene'},
 {type:'scatter3d',mode:'lines',x:[$(js(ex))],y:[$(js(ey))],z:[$(js(ez))],line:{color:'rgba(0,0,0,0.5)',width:1.5},hoverinfo:'skip',showlegend:false,scene:'scene'},
 {type:'scatter3d',mode:'markers',x:[$(js(def[w, 1]))],y:[$(js(def[w, 2]))],z:[$(js(def[w, 3]))],marker:{color:'#e34948',size:3.5},showlegend:false,hoverinfo:'skip',scene:'scene'},
 {type:'scatter3d',mode:'lines',x:[$(js(bx))],y:[$(js(by))],z:[$(js(bz))],line:{color:'#008300',width:7},showlegend:false,hoverinfo:'skip',scene:'scene'},
 {type:'scatter3d',mode:'lines',x:[$(js(ax_))],y:[$(js(ay_))],z:[$(js(az_))],line:{color:'#eb6834',width:4},showlegend:false,hoverinfo:'skip',scene:'scene'},
 {type:'cone',x:[$(js(cx))],y:[$(js(cy))],z:[$(js(cz))],u:[$(js(zeros(length(cx))))],v:[$(js(zeros(length(cx))))],w:[$(js(cw))],
  sizemode:'absolute',sizeref:0.55,anchor:'tip',colorscale:[[0,'#eb6834'],[1,'#eb6834']],showscale:false,showlegend:false,hoverinfo:'skip',scene:'scene'},
 $(join(sec_traces, ",\n "))
];
const layout = {
 title:{text:'$(title)',x:0.02,xanchor:'left',font:{size:15}},
 scene:{domain:{x:[0,0.58],y:[0,1]},aspectmode:'manual',aspectratio:{x:$(ar[1]),y:$(ar[2]),z:$(ar[3])},xaxis:{visible:false},yaxis:{visible:false},zaxis:{visible:false},
        camera:{projection:{type:'orthographic'},eye:{x:-2.0,y:-2.4,z:0.75},center:{x:0,y:0,z:0},up:{x:0,y:0,z:1}},dragmode:'orbit'},
 xaxis:{domain:[0.62,0.98],title:{text:'Y (in), cross-aisle'},scaleanchor:'y',scaleratio:1,zeroline:false},
 yaxis:{domain:[0.2,0.8],title:{text:'X (in), downaisle'},zeroline:false},
 showlegend:false,margin:{l:30,r:20,t:60,b:30},autosize:true,paper_bgcolor:'#fff'};
Plotly.newPlot('plot', data, layout, {responsive:true, displaylogo:false});
</script></body></html>
"""
tmp = joinpath(mktempdir(), "mode_" * tag * "_plotly.html"); write(tmp, html)
out = joinpath(@__DIR__, "mode_" * tag * "_plotly.html"); cp(tmp, out; force = true)
println("wrote ", out, "  (", round(filesize(out) / 1e6, digits = 1), " MB)")
