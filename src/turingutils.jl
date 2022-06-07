q025(x) = quantile(x, 0.025)
q975(x) = quantile(x, 0.975)

function compute_stats(chain::Chains, name::Symbol)
    # x = Array(x)
    # x = DataFrame(x)
    x = Array(group(chain, name))
    # mean_x = mean(x, dims=1)
    # xlb = mean_x - mapslices(x -> quantile(x, 0.025), x, dims=1)
    # xub = mean_x - mapslices(x -> quantile(x, 0.975), x, dims=1)

    mean_x = mean.(eachslice(x, dims=2))
    xlb = mean_x - q025.(eachslice(x, dims=2))
    xub = mean_x - q975.(eachslice(x, dims=2))

    return (μ=mean_x, lb=xlb, ub=xub)
end

function compute_stats(x::Array)
    # mean_x = mean(x, dims=1)
    # xlb = mean_x - mapslices(x -> quantile(x, 0.025), x, dims=1)
    # xub = mean_x - mapslices(x -> quantile(x, 0.975), x, dims=1)

    mean_x = mean.(eachslice(x, dims=2))
    xlb = mean_x .- q025.(eachslice(x, dims=2))
    xub = mean_x .- q975.(eachslice(x, dims=2))

    return (μ=mean_x, lb=xlb, ub=xub)
end

function plot_estimate(s, x::Array)
    μ, lb, ub = compute_stats(x)
    lines(s.timesteps, μ, color=(:steelblue3))
    band!(s.timesteps, μ - lb, μ - ub, color=(:grey68, 0.5))
    current_figure()
end

function get_estimates(chain::Chains)
    wu = Array(group(chain, :ωu))
    we = Array(group(chain, :ωe))
    w = wu .+ we
    return (w=w, wu=wu, we=we)
end


function get_omega(chain::Chains)
    wu = Array(group(chain, :ωu))
    we = Array(group(chain, :ωe))
    w = wu .+ we
    return (w=w, wu=wu, we=we)
end


function summarize_omega(chain::Chains)
    wu_array = Array(group(chain, :ωu))
    we_array = Array(group(chain, :ωe))
    w_array = Array(group(chain, :ω))

    w = compute_stats(w_array)
    wu = compute_stats(wu_array)
    we = compute_stats(we_array)
    return (w=w, wu=wu, we=we)
end