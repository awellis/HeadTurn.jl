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
    w1 = Array(group(chain, :ω1))
    w2 = Array(group(chain, :ω2))
    w = w1 .+ w2
    return (w = w, w1 = w1, w2 = w2)
end


function summarize_omega(chain::Chains)
    w1_array = Array(group(chain, :ω1))
    w2_array = Array(group(chain, :ω2))
    w_array = Array(group(chain, :ω))

    w = compute_stats(w_array)
    w1 = compute_stats(w1_array)
    w2 = compute_stats(w2_array)
    return (w = w, w1 = w1, w2 = w2)
end