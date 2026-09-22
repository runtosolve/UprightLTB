# make_torsion_html.jl — standalone, browser-interactive 3D page of the static twist study of the two-C welded
# upright (plotly.js): deformed shell under the end twist, welds, the fixed end (z = 0: all nodes fixed in X and Y,
# twist and translation restrained, warping free) and the applied rigid twist βo at z = L about the reference
# point (1.5, 1.5) (Abaqus kinematic coupling), warping free. Inset: twist of each C along the member.
# Inputs: torsion_deformed_r5.csv, torsion_summary_r5.csv (torsion_deformed_shape.jl), twist_profile_built_up.csv.
# Run:  julia --project=. make_torsion_html.jl
using DelimitedFiles, Statistics, Printf

d = readdlm(joinpath(@__DIR__, "torsion_deformed_r5.csv"), ','; skipstart = 1)
sm = readdlm(joinpath(@__DIR__, "torsion_summary_r5.csv"), ','; skipstart = 1)
L, weld_length, n_welds, T, J_eff, βo, Xc, Yc = Float64.(sm[1, :])
weld_spacing = 18.0; B = 3.0; D = 3.0 - 0.074; R = 0.199
sid = Int.(d[:, 1]); iid = Int.(d[:, 2]); jid = Int.(d[:, 3])
xyz = Float64.(d[:, 4:6]); u = Float64.(d[:, 7:9])
nn = maximum(iid); nz = maximum(jid)
sc = 0.30 / βo                                     # amplify so the end section is twisted by ≈ 0.30 rad (17°)
def = xyz .+ sc .* u
uin = hypot.(u[:, 1], u[:, 2]); cval = uin ./ maximum(uin)
idx = Dict((sid[k], iid[k], jid[k]) => k for k in eachindex(sid))
I = Int[]; J = Int[]; K = Int[]
for jj in 1:nz-1, ss in 1:2, ii in 1:nn-1
    a = idx[(ss, ii, jj)]; b = idx[(ss, ii, jj + 1)]; c = idx[(ss, ii + 1, jj + 1)]; e = idx[(ss, ii + 1, jj)]
    push!(I, a - 1); push!(J, b - 1); push!(K, c - 1); push!(I, a - 1); push!(J, c - 1); push!(K, e - 1)
end
ex = Any[]; ey = Any[]; ez = Any[]
function seg!(a, b)
    push!(ex, def[a, 1]); push!(ey, def[a, 2]); push!(ez, def[a, 3]); push!(ex, def[b, 1]); push!(ey, def[b, 2]); push!(ez, def[b, 3])
    push!(ex, nothing); push!(ey, nothing); push!(ez, nothing)
end
for jj in 1:nz, ss in 1:2, ii in 1:nn-1; seg!(idx[(ss, ii, jj)], idx[(ss, ii + 1, jj)]); end
for jj in 1:nz-1, ss in 1:2, ii in 1:nn; seg!(idx[(ss, ii, jj)], idx[(ss, ii, jj + 1)]); end
# welds (r5 layout: centers 1.5, 19.5, …, flush with the ends)
xp = xyz[:, 1] .- (sid .- 1) .* B
centers = range(weld_length / 2, L - weld_length / 2, Int(n_welds))
near(x, x0, tol) = abs(x - x0) <= tol + 1e-9
corner = ((sid .== 1) .& near.(xp, B, R) .& (near.(xyz[:, 2], D, R) .| near.(xyz[:, 2], 0.0, R))) .|
         ((sid .== 2) .& near.(xp, 0.0, R) .& (near.(xyz[:, 2], D, R) .| near.(xyz[:, 2], 0.0, R)))
inweld = [any(near(z, c, weld_length / 2) for c in centers) for z in xyz[:, 3]]
w = findall(corner .& inweld .& (xyz[:, 3] .> 0) .& (xyz[:, 3] .< L))
# fixed end z = 0 (all nodes fixed in X, Y) and twisted end z = L (rigid rotation about (Xc, Yc))
e0 = findall(xyz[:, 3] .≈ 0.0); eL = findall(xyz[:, 3] .≈ L)
ox = Any[]; oy = Any[]; oz = Any[]                 # section outline at the fixed end
for ss in 1:2
    rows = [idx[(ss, ii, 1)] for ii in 1:nn]
    append!(ox, xyz[rows, 1]); append!(oy, xyz[rows, 2]); append!(oz, xyz[rows, 3]); push!(ox, nothing); push!(oy, nothing); push!(oz, nothing)
end
eLs = eL[1:3:end]
ax_ = Any[]; ay_ = Any[]; az_ = Any[]; cx = Float64[]; cy = Float64[]; cz = Float64[]; cu = Float64[]; cv = Float64[]
for k in eLs
    rx = xyz[k, 1] - Xc; ry = xyz[k, 2] - Yc; rr = hypot(rx, ry); rr < 0.3 && continue
    tx = -ry / rr; ty = rx / rr                      # tangential direction of the applied rotation
    x0, y0 = def[k, 1], def[k, 2]; len = 1.1
    push!(ax_, x0); push!(ay_, y0); push!(az_, L); push!(ax_, x0 + len * tx); push!(ay_, y0 + len * ty); push!(az_, L)
    push!(ax_, nothing); push!(ay_, nothing); push!(az_, nothing)
    push!(cx, x0 + len * tx); push!(cy, y0 + len * ty); push!(cz, L); push!(cu, tx); push!(cv, ty)
end
js(v) = join((x === nothing ? "null" : string(round(x, digits = 5)) for x in v), ',')
rx_ = maximum(def[:, 1]) - minimum(def[:, 1]); ry_ = maximum(def[:, 2]) - minimum(def[:, 2])
ar = round.([rx_, ry_, L] ./ max(rx_, ry_, L) .* 1.9; digits = 3)

html = """
<!doctype html>
<html><head><meta charset="utf-8"><title>Built-up upright static twist (J_eff)</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>html,body{margin:0;height:100%;font-family:Helvetica,Arial,sans-serif;background:#fff}#plot{width:min(100vw,1400px);height:min(90vh,760px);margin:0 auto;overflow:hidden}</style></head>
<body><div id="plot"></div>
<script>
const data = [
 {type:'mesh3d',x:[$(js(def[:, 1]))],y:[$(js(def[:, 2]))],z:[$(js(def[:, 3]))],i:[$(join(I, ','))],j:[$(join(J, ','))],k:[$(join(K, ','))],
  intensity:[$(js(cval))],colorscale:'Viridis',cmin:0,cmax:1,flatshading:true,lighting:{ambient:0.9,diffuse:0.2,specular:0.0},
  colorbar:{title:{text:'normalized<br>in-plane<br>displacement'},len:0.5,x:0.9},hoverinfo:'skip',showlegend:false,scene:'scene'},
 {type:'scatter3d',mode:'lines',x:[$(js(ex))],y:[$(js(ey))],z:[$(js(ez))],line:{color:'rgba(0,0,0,0.35)',width:1},hoverinfo:'skip',showlegend:false,scene:'scene'},
 {type:'scatter3d',mode:'markers',x:[$(js(def[w, 1]))],y:[$(js(def[w, 2]))],z:[$(js(def[w, 3]))],marker:{color:'#e34948',size:3.5},showlegend:false,hoverinfo:'skip',scene:'scene'},
 {type:'scatter3d',mode:'lines',x:[$(js(ox))],y:[$(js(oy))],z:[$(js(oz))],line:{color:'#0b0b0b',width:5},showlegend:false,hoverinfo:'skip',scene:'scene'},
 {type:'scatter3d',mode:'lines',x:[$(js(ax_))],y:[$(js(ay_))],z:[$(js(az_))],line:{color:'#eb6834',width:4},showlegend:false,hoverinfo:'skip',scene:'scene'},
 {type:'cone',x:[$(js(cx))],y:[$(js(cy))],z:[$(js(cz))],u:[$(js(cu))],v:[$(js(cv))],w:[$(js(zeros(length(cx))))],
  sizemode:'absolute',sizeref:0.45,anchor:'tip',colorscale:[[0,'#eb6834'],[1,'#eb6834']],showscale:false,showlegend:false,hoverinfo:'skip',scene:'scene'},
 {type:'scatter3d',mode:'markers',x:[$(Xc)],y:[$(Yc)],z:[$(L)],marker:{color:'#eb6834',size:6,symbol:'diamond'},showlegend:false,hoverinfo:'skip',scene:'scene'}
];
const layout = {
 title:{text:'Static twist of the two-C welded upright, L = $(Int(L)) in, $(Int(n_welds)) welds × $(Int(weld_length)) in at $(Int(weld_spacing)) in: z = 0 fixed in X, Y (black outline: twist and translation restrained, warping free); rigid twist β<sub>o</sub> applied at z = L about (1.5, 1.5) (orange), warping free.  J<sub>eff</sub> = T L / (G β<sub>o</sub>) = $(round(J_eff, digits = 3)) in⁴',x:0.02,xanchor:'left',font:{size:14}},
 scene:{domain:{x:[0,1],y:[0,1]},aspectmode:'manual',aspectratio:{x:$(ar[1]),y:$(ar[2]),z:$(ar[3])},xaxis:{visible:false},yaxis:{visible:false},zaxis:{visible:false},
        camera:{projection:{type:'orthographic'},eye:{x:-2.0,y:-2.4,z:0.75},center:{x:0,y:0,z:0},up:{x:0,y:0,z:1}},dragmode:'orbit'},
 showlegend:false,margin:{l:30,r:20,t:60,b:30},autosize:true,paper_bgcolor:'#fff'};
Plotly.newPlot('plot', data, layout, {responsive:true, displaylogo:false});
</script></body></html>
"""
tmp = joinpath(mktempdir(), "torsion_r5_plotly.html"); write(tmp, html)
out = joinpath(@__DIR__, "torsion_r5_plotly.html"); cp(tmp, out; force = true)
println("wrote ", out, "  (", round(filesize(out) / 1e6, digits = 1), " MB)")
