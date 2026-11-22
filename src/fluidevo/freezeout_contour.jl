# freeze-out surfaces from fields given as callables of (τ, r): Milne, boost invariant, azimuthally symmetric
# a surface is parameterised by its normalised arc length, so it may have an outer and an inner edge at one τ

"""
    FOContourOpts(; Tfo = 0.156, rmax = 10.0, nr = 4000, tfo_margin = 3.0e-4, max_recession = 14.0,
                    min_inner_run = 5, close_axis = true)
Options of `isotherm_contour`: Tfo (GeV), the radial scan [0, rmax] with nr points, the hot-matter margin (GeV),
the largest inward speed of the outer edge, the shortest inner branch kept, and whether the path closes at r = 0.
"""
Base.@kwdef struct FOContourOpts
    Tfo::Float64           = 0.156
    rmax::Float64          = 10.0
    nr::Int                = 4000
    tfo_margin::Float64    = 3.0e-4
    max_recession::Float64 = 14.0
    min_inner_run::Int     = 5
    close_axis::Bool       = true
end

# radii where T(τ, r) crosses Tfo, ascending, by linear interpolation
function isotherm_crossings(T, tau, o::FOContourOpts)
    rs = range(0.0, o.rmax; length = o.nr)
    out = Float64[]
    Tp = Float64(T(tau, rs[1])) - o.Tfo
    for k in 2:length(rs)
        Tk = Float64(T(tau, rs[k])) - o.Tfo
        if (Tp < 0) != (Tk < 0)
            push!(out, rs[k-1] + (rs[k] - rs[k-1]) * (-Tp) / (Tk - Tp + 1e-30))
        end
        Tp = Tk
    end
    return out
end

# max_r T(τ, r) - Tfo
function peak_excess(T, tau, o::FOContourOpts)
    m = -Inf
    for r in range(0.0, o.rmax; length = o.nr)
        e = Float64(T(tau, r)) - o.Tfo
        e > m && (m = e)
    end
    return m
end

# the τ at which T(τ, 0) = Tfo (NaN if the axis does not cross Tfo on the grid)
function axis_freeze_time(T, taus, o::FOContourOpts)
    prev = Float64(T(taus[1], 0.0)) - o.Tfo
    for k in 2:length(taus)
        cur = Float64(T(taus[k], 0.0)) - o.Tfo
        if (prev < 0) != (cur < 0)
            return taus[k-1] + (taus[k] - taus[k-1]) * (-prev) / (cur - prev)
        end
        prev = cur
    end
    return NaN
end

"""
    isotherm_branches(T, taus; kwargs...)
The outer and inner edges of the hot region T(τ, r) > Tfo: `(outer_tau, outer_r, inner_tau, inner_r, tau_pinch)`.
"""
function isotherm_branches(T, taus, o::FOContourOpts)
    outer_tau = Float64[]; outer_r = Float64[]
    inner_tau = Float64[]; inner_r = Float64[]
    tau_pinch = NaN; prev_tau = NaN
    for tau in taus
        if peak_excess(T, tau, o) < o.tfo_margin
            tau_pinch = tau
            break
        end
        cs = isotherm_crossings(T, tau, o)
        if isempty(cs)
            tau_pinch = tau
            break
        end
        cand = maximum(cs)
        # an inward jump faster than max_recession: the outer edge is gone
        if !isempty(outer_r) && cand < last(outer_r) - o.max_recession * (tau - prev_tau)
            tau_pinch = tau
            break
        end
        push!(outer_tau, tau); push!(outer_r, cand); prev_tau = tau
        if length(cs) >= 2 && Float64(T(tau, 0.0)) < o.Tfo
            push!(inner_tau, tau); push!(inner_r, minimum(cs))
        end
    end
    if isnan(tau_pinch) && !isempty(outer_tau)
        @warn "isotherm_branches: hot matter is left at the last τ of the grid; the surface is open" tau_last = last(outer_tau)
        tau_pinch = last(outer_tau)
    end
    return (; outer_tau, outer_r, inner_tau, inner_r, tau_pinch)
end

isotherm_branches(T, taus; kwargs...) = isotherm_branches(T, taus, FOContourOpts(; kwargs...))

"""
    isotherm_contour(T, taus; kwargs...)
The isotherm T(τ, r) = Tfo as one ordered path `(lambda, tau, r, branch)`: the outer branch, then the inner one.
"""
function isotherm_contour(T, taus, o::FOContourOpts)
    b = isotherm_branches(T, taus, o)
    isempty(b.outer_tau) && return (lambda = Float64[], tau = Float64[], r = Float64[], branch = Symbol[])

    path_tau = copy(b.outer_tau); path_r = copy(b.outer_r)
    branch = fill(:outer, length(path_tau))
    inner_tau, inner_r = b.inner_tau, b.inner_r
    if 0 < length(inner_tau) < o.min_inner_run
        @warn "isotherm_contour: discarding an inner branch of $(length(inner_tau)) points as discretisation noise (min_inner_run = $(o.min_inner_run))" tau = extrema(inner_tau) r = extrema(inner_r)
        inner_tau = Float64[]; inner_r = Float64[]
    end
    if !isempty(inner_tau)
        ord = sortperm(inner_tau; rev = true)   # from the pinch back down to the axis
        append!(path_tau, inner_tau[ord]); append!(path_r, inner_r[ord])
        append!(branch, fill(:inner, length(inner_tau)))
    end
    if o.close_axis
        tau_c = axis_freeze_time(T, taus, o)
        if isfinite(tau_c)
            # with a shell the axis froze before the pinch, without one it freezes last
            shell = !isempty(inner_tau)
            if (shell && tau_c < path_tau[end]) || (!shell && tau_c > path_tau[end])
                push!(path_tau, tau_c); push!(path_r, 0.0); push!(branch, shell ? :inner : :outer)
            end
        end
    end
    lambda = zeros(length(path_tau))
    for i in 2:length(path_tau)
        lambda[i] = lambda[i-1] + hypot(path_tau[i] - path_tau[i-1], path_r[i] - path_r[i-1])
    end
    return (; lambda, tau = path_tau, r = path_r, branch)
end

isotherm_contour(T, taus; kwargs...) = isotherm_contour(T, taus, FOContourOpts(; kwargs...))

# Piecewise-linear interpolation of values `v` given at increasing knots `s`.
struct _Polyline{S,V}
    s::S
    v::V
end

function (p::_Polyline)(t)
    s = p.s
    t <= first(s) && return first(p.v)
    t >= last(s) && return last(p.v)
    i = searchsortedlast(s, t)
    w = (t - s[i]) / (s[i+1] - s[i])
    return (1 - w) * p.v[i] + w * p.v[i+1]
end

_as_svector(v::SVector) = v
_as_svector(v::Tuple) = SVector(v)
_as_svector(v::AbstractVector) = SVector{length(v)}(v)

# FreezeOutResult on s ∈ [0, 1] from x(s) and φ(s), both called with a 1-tuple (s,).
function _freezeout_result(x_fn, phi_fn, npts)
    x = SplineInterp(PieceWiseFunction(x_fn, 0.0, 1.0), (npts,))
    phi = SplineInterp(PieceWiseFunction(phi_fn, 0.0, 1.0), (npts,))
    return FreezeOutResult(x, phi)
end

"""
    freezeout_from_path(tau, r, fields; npts = 300) -> FreezeOutResult
A `FreezeOutResult` on the path (tau[i], r[i]) by arc length; `fields` is a callable `fields(τ, r)` or one
field vector per point. Orient the path from large r at τ₀ towards the axis.
"""
function freezeout_from_path(tau::AbstractVector, r::AbstractVector, fields; npts::Integer = 300)
    s, keep = _arclength(tau, r)
    pos = _Polyline(s, [SVector{2,Float64}(tau[i], r[i]) for i in keep])
    x_fn(t) = pos(t[1])
    phi_fn(t) = (p = pos(t[1]); _as_svector(fields(p[1], p[2])))
    return _freezeout_result(x_fn, phi_fn, npts)
end

function freezeout_from_path(tau::AbstractVector, r::AbstractVector, fields::AbstractVector; npts::Integer = 300)
    length(fields) == length(tau) ||
        throw(ArgumentError("fields has $(length(fields)) entries for a path of $(length(tau)) points"))
    s, keep = _arclength(tau, r)
    pos = _Polyline(s, [SVector{2,Float64}(tau[i], r[i]) for i in keep])
    val = _Polyline(s, [_as_svector(fields[i]) for i in keep])
    x_fn(t) = pos(t[1])
    phi_fn(t) = val(t[1])
    return _freezeout_result(x_fn, phi_fn, npts)
end

# normalised arc length and the indices kept (zero-length segments dropped)
function _arclength(tau, r)
    length(tau) == length(r) || throw(ArgumentError("tau and r have different lengths"))
    length(tau) >= 2 || throw(ArgumentError("a path needs at least two points"))
    keep = [firstindex(tau)]
    lam = [0.0]
    for i in firstindex(tau)+1:lastindex(tau)
        j = last(keep)
        d = hypot(tau[i] - tau[j], r[i] - r[j])
        if d > 0
            push!(keep, i); push!(lam, last(lam) + d)
        end
    end
    L = last(lam)
    L > 0 || throw(ArgumentError("the path has zero length"))
    return lam ./ L, keep
end

"""
    freezeout_cap(tau0, r_edge, r_max, fields; npts = 300) -> FreezeOutResult
The τ = τ₀ face for r_edge < r < r_max, the matter already below Tfo at τ₀; spectra add over surface and cap.
"""
function freezeout_cap(tau0::Real, r_edge::Real, r_max::Real, fields; npts::Integer = 300)
    r_max > r_edge || throw(ArgumentError("need r_max > r_edge, got r_edge = $r_edge, r_max = $r_max"))
    t0 = Float64(tau0); ra = Float64(r_edge); rb = Float64(r_max)
    x_fn(t) = SVector{2,Float64}(t0, rb - t[1] * (rb - ra))
    phi_fn(t) = _as_svector(fields(t0, rb - t[1] * (rb - ra)))
    return _freezeout_result(x_fn, phi_fn, npts)
end

"""
    freezeout_from_isotherm(T, fields, taus; npts = 300, cap = true, kwargs...) -> (surface, cap, contour)
`isotherm_contour`, `freezeout_from_path` and `freezeout_cap` in one call.
"""
function freezeout_from_isotherm(T, fields, taus; npts::Integer = 300, cap::Bool = true, kwargs...)
    o = FOContourOpts(; kwargs...)
    contour = isotherm_contour(T, taus, o)
    length(contour.tau) >= 2 || throw(ArgumentError("no isotherm T = $(o.Tfo) found on the grid"))
    surface = freezeout_from_path(contour.tau, contour.r, fields; npts = npts)
    face = nothing
    if cap && first(contour.tau) == first(taus) && first(contour.r) < o.rmax
        face = freezeout_cap(first(taus), first(contour.r), o.rmax, fields; npts = npts)
    end
    return (; surface, cap = face, contour)
end
