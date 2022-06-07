## as anonymous function
@model StateSpaceModel(y) = begin
    N = length(y)
    x = tzeros(Real, N)
    σy = 1.0
    σx = 1.0

    x[1] ~ Normal(0, 1)
    y[1] ~ Normal(x[1], σy)

    # Iterators.drop(1:N, 1)
    for i ∈ 2:N
        x[i] ~ Normal(x[i-1], σx)
        y[i] ~ Normal(x[i], σy)
    end
end

## as function
    # @model function StateSpaceModel(y)
    #     N = length(y)
    #     x = tzeros(Real, N)
    #     σy = 1.0
    #     σx = 1.0

    #     x[1] ~ Normal(0, 1)
    #     y[1] ~ Normal(x[1], σy)
    #     for i ∈ 2:N
    #         x[i] ~ Normal(x[i-1], σx)
    #         y[i] ~ Normal(x[i-1], σy)
    #     end
    # end



@model random_acceleration_model(y; Δt = 0.01,  σω = 0.5, σα = 1.0, σy = 1.0, κ1, κ2) = begin
    N = length(y)
    α = tzeros(Real, N)
    ω = tzeros(Real, N)
    c = zeros(Real, N)
    θ = zeros(Real, N)
    # κ1, κ2 = cupula_dynamics(Δt)

    σα = σω
    # σω ~ truncated(Cauchy(0, 1), 0, Inf)
    # σy ~ truncated(Cauchy(0, 2), 0, Inf)

    α[1] ~ Normal(0, σα)
    ω[1] ~ Normal(0, σω)
    c[1] = κ1 * 0 + κ2 * ω[1]
    y[1] ~ Normal(ω[1] - c[1], σy)

    for i ∈ 2:N
        α[i] ~ Normal(α[i-1], σα)
        ω[i] ~ Normal(ω[i-1] + Δt * α[i], σω)
        c[i] = κ1 * c[i-1] + κ2 * ω[i]
        y[i] ~ Normal(ω[i] - c[i], σy)
    end
    return c
end


@model random_acceleration_model_2(y; Δt = 0.01,  σω = 0.5, σα = 1.0, σy = 1.0, κ1, κ2) = begin
    N = length(y)
    α = tzeros(Real, N)
    ω = tzeros(Real, N)
    c = zeros(Real, N)
    θ = zeros(Real, N)
    # κ1, κ2 = cupula_dynamics(Δt)

    σα = σω
    # σω ~ truncated(Cauchy(0, 1), 0, Inf)
    # σy ~ truncated(Cauchy(0, 2), 0, Inf)

    α[1] ~ Normal(0, σα)
    ω[1] ~ Normal(0, σω)
    c[1] ~ Determin(κ1 * 0 + κ2 * ω[1]) 
    y[1] ~ Normal(ω[1] - c[1], σy)

    for i ∈ 2:N
        α[i] ~ Normal(α[i-1], σα)
        ω[i] ~ Normal(ω[i-1] + Δt * α[i], σω)
        c[i] ~ Determin(κ1 * c[i-1] + κ2 * ω[i])
        y[i] ~ Normal(ω[i] - c[i], σy)
    end
end

@model RandomWalkAngularAccelerationModel(y, s; σω = 1.0) = begin
    N = length(y)
    α = tzeros(Real, N)
    ω = tzeros(Real, N)
    c = zeros(Real, N)
    θ = zeros(Real, N)
    Δt = s.Δt
    κ1, κ2 = cupula_dynamics(Δt)

    σα = σω
    # σω = 0.1
    σy = s.sensor.σ

    α[1] ~ Normal(0, σα)
    ω[1] ~ Normal(0, σω)
    c[1] ~ Determin(κ1 * 0 + κ2 * ω[1])
    y[1] ~ Normal(ω[1] - c[1], σy)

    for i ∈ Iterators.drop(1:N, 1)
        α[i] ~ Normal(α[i-1], σα)
        ω[i] ~ Normal(ω[i-1] + Δt * α[i], σω)
        c[i] ~ Determin(κ1 * c[i-1] + κ2 * ω[i])
        y[i] ~ Normal(ω[i] - c[i], σy)
    end
end


@model RandomWalkAngularVelocityModel(y, s; σω = 1.0) = begin
    N = length(y)
    ω = tzeros(Real, N)
    p = tzeros(Real, N)
    c = zeros(Real, N)
    Δt = s.Δt
    κ1, κ2 = cupula_dynamics(Δt)

    # σω = 1
    σy = s.sensor.σ
    ω[1] ~ Normal(0, σω)
    p[1] ~ Determin(ω[1] > 0 ? 1 : 0)
    c[1] ~ Determin(κ1 * 0 + κ2 * ω[1])
    y[1] ~ Normal(ω[1] - c[1], σy)

    for i ∈ Iterators.drop(1:N, 1)
        ω[i] ~ Normal(ω[i-1], σω)
        p[i] ~ Determin(ω[i] > 0 ? 1 : 0)
        c[i] ~ Determin(κ1 * c[i-1] + κ2 * ω[i])
        y[i] ~ Normal(ω[i] - c[i], σy)
    end
end

@model function headturnmodel(y, N; u, σωu, σωe, σy, timesteps, onsetᵤ, Du, Au, fu, duration_u, Δt, κ1=κ1, κ2=κ2) where {T}
    ωu = tzeros(Real, N)
    ωe = tzeros(Real, N)
    ω = tzeros(Real, N)
    c = tzeros(Real, N)

    ωu[1] ~ Normal(0, σωu)
    ωe[1] ~ Normal(0, σωe)
    ω[1] = ωu[1] + ωe[1]
    c[1] = κ1 * 0 + κ2 * ω[1]
    y[1] ~ Normal(ω[1] - c[1], σy)

    # for i ∈ Iterators.drop(1:N, 1)
    for i ∈ 2:N
        ωu[i] ~ Normal(ωu[i-1] + Δt * u[i], σωu)
        ωe[i] ~ Normal(ωe[i-1], σωe)
        ω[i] = ωu[i] + ωe[i]
        c[i] = κ1 * c[i-1] + κ2 * ω[i]
        y[i] ~ Normal(ω[i]-c[i], σy)
    end
    return (ω, c)
end

function makeheadturnmodel(s; σωu=0.05, σωe=0.1)
    y = s.y
    N = length(y)
    Du = s.mᵤ.D == "left" ? -1 : 1
    Au = s.mᵤ.A
    σy = 1.0 * s.sensor.σ # inflate the measurement noise	
    Δt = s.Δt
    duration_u = s.mᵤ.duration
    # onsetᵤ = hasproperty(s.mᵤ, :onset) ? s.mᵤ.onset : 1.0
    onsetᵤ = s.mᵤ.onset
    timesteps = range(0, stop=N * Δt, step=Δt)
    fu = 1 / (duration_u) # frequency: single sinusoidal head turn

    u = zeros(N)
    for i ∈ 2:N
        t = timesteps[i]
        u[i] = (t >= onsetᵤ) && (t <= onsetᵤ + duration_u) ? acceleration(Du, Au, fu, t, onsetᵤ) : 0.0
    end

    κ1, κ2 = cupula_dynamics(Δt)

    model = headturnmodel(y, N, u=u, σωu=σωu, σωe=σωe, σy=σy, timesteps=timesteps, onsetᵤ=onsetᵤ, Du=Du, Au=Au, fu=fu, duration_u=duration_u, Δt=Δt, κ1=κ1, κ2=κ2)
end





@model HeadTurnOnlyParamModel(y, s, m; σωᵤ = 0.01, σₐ = 1.0) = begin
    N = length(y)
    
    # intialize vectors
    α1 = tzeros(Real, N)
    ω = tzeros(Real, N)
    p = tzeros(Real, N)
    c = zeros(Real, N)

    # parameters
    D1 = m.D == "left" ? -1 : 1
    A1 ~ Normal(m.A, σₐ)

    Δt = s.Δt
    onsetᵤ = hasproperty(m, :onset) ? m.onset : s.onsetᵤ
    duration1 = m.duration
    f1 = 1/(duration1) # frequency: single sinusiodal head turn
    timesteps = range(0, stop = N, step = Δt)
    κ1, κ2 = cupula_dynamics(Δt)

    # priors
    σy = s.sensor.σ
    
    # initial state at t_0
    ω[1] ~ Normal(0, σωᵤ)
    # p[1] ~ Determin(ω[1] > 0 ? 1 : 0)
    c[1] ~ Determin(κ1 * 0 + κ2 * ω[1])
    y[1] ~ Normal(ω[1] - c[1], σy)

    for i ∈ Iterators.drop(1:N, 1)
        t = round(timesteps[i]; digits = 3)
        α1[i] = (t > onsetᵤ) & (t < onsetᵤ + duration1) ? acceleration(D1, A1, f1, t, onsetᵤ, noisy = false) : 0.0

        ω[i] ~ Normal(ω[i-1] + Δt * α1[i], σωᵤ)
        # p[i] ~ Determin(ω[i] > 0 ? 1 : 0)
        c[i] ~ Determin(κ1 * c[i-1] + κ2 * ω[i])
        y[i] ~ Normal(ω[i] - c[i], σy)
    end
end

@model headturn_simulation_model(y, s, m; Δt = 0.01, σωᵤ = 0.01, 
                                 σₐ = 1.0, onset = 1, κ1, κ2) = begin
    N = length(y)
    
    # intialize vectors
    α1 = tzeros(Real, N)
    ω = tzeros(Real, N)
    p = tzeros(Real, N)
    c = zeros(Real, N)

    # parameters
    D = m.D == "left" ? -1 : 1
    A ~ Normal(m.A, σₐ)

    # Δt = s.Δt
    # on = hasproperty(m, :onset) ? m.onset : s.onsetᵤ
    # on = shape * scale
    onsetᵤ ~ Gamma(20 * onset, 0.05)
    duration1 = m.duration
    f = 1/(duration1) # frequency: single sinusiodal head turn
    timesteps = range(0, stop = N, step = Δt)
    # κ1, κ2 = cupula_dynamics(Δt)

    # priors
    σy = s.sensor.σ
    
    # initial state at t_0
    ω[1] ~ Normal(0, σωᵤ)
    # p[1] ~ Determin(ω[1] > 0 ? 1 : 0)
    c[1] = κ1 * 0 + κ2 * ω[1]
    y[1] ~ Normal(ω[1] - c[1], σy)

    for i ∈ Iterators.drop(1:N, 1)
        t = round(timesteps[i]; digits = 3)
        α1[i] = (t > onsetᵤ) & (t < onsetᵤ + duration1) ? acceleration(D, A, f, t, onsetᵤ, noisy = false) : 0.0

        ω[i] ~ Normal(ω[i-1] + Δt * α1[i], σωᵤ)
        # p[i] ~ Determin(ω[i] > 0 ? 1 : 0)
        c[i] = κ1 * c[i-1] + κ2 * ω[i]
        y[i] ~ Normal(ω[i] - c[i], σy)
    end
    return c
end

@model HeadTurnOnlySimulationModel(y, s, m; σωᵤ = 0.01, σₐ = 1.0) = begin
    N = length(y)
    
    # intialize vectors
    α1 = tzeros(Real, N)
    ω = tzeros(Real, N)
    p = tzeros(Real, N)
    c = zeros(Real, N)

    # parameters
    D = m.D == "left" ? -1 : 1
    A ~ Normal(m.A, σₐ)

    Δt = s.Δt
    on = hasproperty(m, :onset) ? m.onset : s.onsetᵤ
    # on = shape * scale
    onsetᵤ ~ Gamma(20 * on, 0.05)
    duration1 = m.duration
    f = 1/(duration1) # frequency: single sinusiodal head turn
    timesteps = range(0, stop = N, step = Δt)
    κ1, κ2 = cupula_dynamics(Δt)

    # priors
    σy = s.sensor.σ
    
    # initial state at t_0
    ω[1] ~ Normal(0, σωᵤ)
    # p[1] ~ Determin(ω[1] > 0 ? 1 : 0)
    c[1] ~ Determin(κ1 * 0 + κ2 * ω[1])
    y[1] ~ Normal(ω[1] - c[1], σy)

    for i ∈ Iterators.drop(1:N, 1)
        t = round(timesteps[i]; digits = 3)
        α1[i] = (t > onsetᵤ) & (t < onsetᵤ + duration1) ? acceleration(D, A, f, t, onsetᵤ, noisy = false) : 0.0

        ω[i] ~ Normal(ω[i-1] + Δt * α1[i], σωᵤ)
        # p[i] ~ Determin(ω[i] > 0 ? 1 : 0)
        c[i] ~ Determin(κ1 * c[i-1] + κ2 * ω[i])
        y[i] ~ Normal(ω[i] - c[i], σy)
    end
end


# TODO: when moving, set σωᵤ much lower, when not moving higher
@model HeadTurnComplexModel(y, s, m; σωᵤ = 0.01, σωₑ = 0.1) = begin
    N = length(y)
    stochastic = s.stochastic
    # intialize vectors
    α1 = tzeros(Real, N)
    ω1 = tzeros(Real, N)
    p1 = tzeros(Real, N)
    ω2 = tzeros(Real, N)
    p2 = tzeros(Real, N)
    ω = tzeros(Real, N)
    move = tzeros(Real, N)
    p = tzeros(Real, N)
    c = zeros(Real, N)

    # parameters
    D1 = m.D == "left" ? -1 : 1
    A1 = stochastic == 0 ? m.A : rand(Normal(m1.A, 0.5))

    Δt = s.Δt
    onsetᵤ = hasproperty(m, :onset) ? m.onset : s.onsetᵤ
    duration1 = m.duration
    f1 = 1/(duration1) # frequency: single sinusiodal head turn
    timesteps = range(0, stop = N, step = Δt)
    κ1, κ2 = cupula_dynamics(Δt)

    # priors
    σy = s.sensor.σ
    
    # initial state at t_0
    ω1[1] ~ Normal(0, σωᵤ)
    p1[1] ~ Determin(ω1[1] > 0 ? 1 : 0)
    ω2[1] ~ Normal(0, σωₑ)
    p2[1] ~ Determin(ω2[1] > 0 ? 1 : 0)
    ω[1] ~ Determin(ω1[1] + ω2[1])
    move[1] ~ Determin(0)
    p[1] ~ Determin(ω[1] > 0 ? 1 : 0)
    c[1] ~ Determin(κ1 * 0 + κ2 * ω[1])
    y[1] ~ Normal(ω[1] - c[1], σy)

    for i ∈ Iterators.drop(1:N, 1)
        t = round(timesteps[i]; digits = 3)
        α1[i] = (t > onsetᵤ) & (t < onsetᵤ + duration1) ? acceleration(D1, A1, f1, t, onsetᵤ, noisy = stochastic) : 0.0

        ω1[i] ~ Normal(ω1[i-1] + Δt * α1[i], σωᵤ)
        move[i] ~ Determin((t > onsetᵤ) & (t < onsetᵤ + duration1) ? 1 : 0)
        p1[i] ~ Determin(ω1[i] > 0 ? 1 : 0)
        ω2[i] ~ Normal(ω2[i-1], σωₑ)
        p2[i] ~ Determin(ω2[i] > 0 ? 1 : 0)
        ω[i] ~ Determin(ω1[i] + ω2[i])
        p[i] ~ Determin(ω[i] > 0 ? 1 : 0)

        c[i] ~ Determin(κ1 * c[i-1] + κ2 * ω[i])
        y[i] ~ Normal(ω[i] - c[i], σy)
    end
end


@model CoupledModel(y, s, m; σωᵤ = 0.01, σωₑ = 0.1) = begin
    N = length(y)
    stochastic = s.stochastic
    # intialize vectors
    α1 = tzeros(Real, N)
    ω1 = tzeros(Real, N)
    # p1 = tzeros(Real, N)
    ω2 = tzeros(Real, N)
    p2 = tzeros(Real, N)
    ω = tzeros(Real, N)
    # p = tzeros(Real, N)
    c = zeros(Real, N)

    # parameters
    D1 = m.D == "left" ? -1 : 1
    A1 = stochastic == 0 ? m.A : rand(Normal(m1.A, 0.5))

    Δt = s.Δt
    onsetᵤ = hasproperty(m, :onset) ? m.onset : s.onsetᵤ
    duration1 = m.duration
    f1 = 1/(duration1) # frequency: single sinusiodal head turn
    timesteps = range(0, stop = N, step = Δt)
    κ1, κ2 = cupula_dynamics(Δt)

    # priors
    σy = s.sensor.σ
    
    # initial state at t_0
    ω1[1] ~ Normal(0, σωᵤ)
    # p1[1] ~ Determin(ω1[1] > 0 ? 1 : 0)
    ω2[1] ~ Normal(0, σωₑ)
    # p2[1] ~ Determin(ω2[1] > 0 ? 1 : 0)
    ω[1] ~ Determin(ω1[1] - ω2[1])
    # p[1] ~ Determin(ω[1] > 0 ? 1 : 0)
    c[1] ~ Determin(κ1 * 0 + κ2 * ω[1])
    y[1] ~ Normal(ω1[1] - c[1], σy)
    y[1] ~ Normal(ω2[1] - c[1], σy)

    for i ∈ Iterators.drop(1:N, 1)
        t = round(timesteps[i]; digits = 3)
        α1[i] = (t > onsetᵤ) & (t < onsetᵤ + duration1) ? acceleration(D1, A1, f1, t, onsetᵤ, noisy = stochastic) : 0.0

        ω1[i] ~ Normal(ω1[i-1] + Δt * α1[i], σωᵤ)
        # p1[i] ~ Determin(ω1[i] > 0 ? 1 : 0)
        ω2[i] ~ Normal(ω2[i-1], σωₑ)
        # p2[i] ~ Determin(ω2[i] > 0 ? 1 : 0)
        ω[i] ~ Determin(ω1[i] - ω2[i])
        # p[i] ~ Determin(ω[i] > 0 ? 1 : 0)

        c[i] ~ Determin(κ1 * c[i-1] + κ2 * ω[i])
        y[i] ~ Normal(ω1[i] - c[i], σy)
        y[i] ~ Normal(ω2[i] - c[i], σy)
    end
end
