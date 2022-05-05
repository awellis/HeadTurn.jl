function compute_stats(chain::Chains, name::Symbol)
    # x = Array(x)
    # x = DataFrame(x)
    x = Array(group(chain, name))
    mean_x = mean(x, dims=1)
    xlb = mean_x - mapslices(x -> quantile(x, 0.025), x, dims=1)
    xub = mean_x - mapslices(x -> quantile(x, 0.975), x, dims=1)
    return (μ = mean_x, lb = xlb, ub = xub)
end

function compute_stats(x::Array)
    mean_x = mean(x, dims=1)
    xlb = mean_x - mapslices(x -> quantile(x, 0.025), x, dims=1)
    xub = mean_x - mapslices(x -> quantile(x, 0.975), x, dims=1)
    return (μ = mean_x, lb = xlb, ub = xub)
end

function get_omega(chain::Chains)   
    wu = Array(group(chain, :ωu))
    we = Array(group(chain, :ωe))
    w = wu .+ we
    return (w = w, wu = wu, we = we)
end


function summarize_omega(chain::Chains)
    wu_array = Array(group(chain, :ωu))
    we_array = Array(group(chain, :ωe))
    w_array = Array(group(chain, :ω))

    w = compute_stats(w_array)
    wu = compute_stats(wu_array)
    we = compute_stats(we_array)
    return (w = w, wu = wu, we = we)
end