using Parameters, Distributions, StatsPlots

@with_kw struct Sensor
    noise::Distribution = Normal(0, 1.0)
    # κ1, κ2 = cupula_dynamics(Δt)
end

@with_kw struct NigmatullinaTrial
    D::String # direction
    onset::Real = 0.0
end

@with_kw struct HeadTurn
    A::Real = 20 # amplitude
    D::String = "right" # direction
    onset::Real = 1
    duration::Real = 1
end

@with_kw struct PlannedMovement
    A::Real = 20 # amplitude
    D::String = "right" # direction
    onset::Real = 1
    duration::Real = 1
end

@with_kw struct ExternalDisturbance
    A::Real = 10 # amplitude
    D::String = "right" # direction
    onset::Real = 1
    duration::Real = 2
	@assert A > zero(A)
end


@with_kw struct Translation
    A::Real = 20 # amplitude
    D::String = "right" # direction
    onset::Real = 1
    duration::Real = 1
end


function cupula_dynamics(Δt::Real; τ::Real = 4)
    # τ = 4 same as Laurens & Angelaki (2017)
    κ₁ = τ / (τ + Δt)
    κ₂ = Δt / (τ + Δt)
    return κ₁, κ₂
end

function velocity(D::Number, A::Number, f::Number, t::Number, start::Number)
    D * A * 1/(2*π*f) * (1-cos.(2*π*f*(t-start)))
end

function acceleration(D::Number, A::Number, f::Number, t::Number, start::Number;
    noisy::Bool = false)
    noise = noisy == 0 ? 0 : rand(Normal(0, 2.0))
    D * A * sin.(2*π*f*(t-start)) + noise 
end

function displacement(A::Number; T::Number = 1)
    A * T^2 /(2π)
end 

function amplitude(θ::Number; T::Number = 1)
    (2π/T^2) * θ
end

function peak_velocity(A::Number; T::Number = 1)
    (T/π) * A
end


function simulate(m::NigmatullinaTrial, sensor::Sensor;
    onset::Real = 0, trial_duration::Real = 15, Δt = 0.01)
    #= An alternative “limit-finding” approach was adopted with each trial applying
    progressively increasing acceleration in the yaw plane, increasing by 0.5 deg/s^2
    every three seconds, reaching 82 deg/sec

    Cutfield, N. J., Cousins, S., Seemungal, B. M., Gresty, M. A., &
    Bronstein, A. M. (2011). Vestibular perceptual thresholds to angular
    rotation in acute unilateral vestibular paresis and with galvanic
    stimulation. Annals of the New York Academy of Sciences, 1233(1), 256–262.
    =#

    @assert onset >= 0
    @assert trial_duration > 0

    m.D ∈ ["left", "right"] || error("Direction not specified.")
    D = m.D == "left" ? -1 : 1
    noise = sensor.noise
    total_duration = onset + trial_duration

    timesteps = range(0, stop = total_duration, step = Δt)
    T = length(timesteps)

    κ1, κ2 = cupula_dynamics(Δt)

    # either use undef intialiser, or intialise with zeros
    # ω = Array{Float64}(undef, nsteps) # velocity
    # θ = Array{Float64}(undef, nsteps) # position
    # c = Array{Float64}(undef, nsteps) # cupula dynamics
    # y = Array{Float64}(undef, nsteps) # sensory observations
    α = zeros(T)
    ω = zeros(T)
    θ = zeros(T)
    c = zeros(T)
    y = zeros(T)

    increment = D * 0.3
    tstep = 3

    # start until motion onset
    j = 1
    t = round(timesteps[j]; digits = 3)

    while t < onset
        y[j] = ω[j] - c[j] + rand(noise)
        j += 1
        t = round(timesteps[j]; digits = 3)
    end

    k = 1 # counts the number of times that α is incremented
    # start moving at time j
    α[j] = α[j] + increment

    for i ∈ (j+1):T
        t = round(timesteps[i]; digits = 3)

        if t >= (onset + k * tstep)
            α[i] = α[i-1] + increment
            k += 1
        else
            α[i] = α[i-1]
        end
        ω[i] = ω[i-1] + Δt * α[i-1]
        θ[i] = θ[i-1] + Δt * ω[i-1] + 1/2 * Δt^2 * α[i-1]
        c[i] = κ1 * c[i-1] + κ2 * ω[i]
        y[i] = ω[i] - c[i] + rand(noise)
    end

    # return a named tuple
    out = (timesteps = collect(timesteps),
            y = y, α = α, ω = ω, θ = θ, c = c, 
            Δt = Δt, onset = onset)
    return out
end

function simulate(mᵤ::HeadTurn,
    mₑ::ExternalDisturbance,
    sensor::Sensor;
    Δt = 0.01,
    duration::Real = 5, 
    stochastic::Bool = false)

    mᵤ = @isdefined(mᵤ) ? mᵤ : HeadTurn(A = 0)
    mₑ = @isdefined(mₑ) ? mₑ : ExternalDisturbance(A = 0)

    onsetᵤ = mᵤ.onset
    onsetₑ = mₑ.onset
    @assert onsetᵤ >= 0
    @assert onsetₑ >= 0
    @assert duration > 0

    # onset = m.onset
    duration1 = mᵤ.duration
    duration2 = mₑ.duration
    total_duration = duration

    onsetᵤ + duration1 <= total_duration || error("Head turn extends beyond simulation event.")
    onsetₑ + duration2 <= total_duration || error("Disturbance extends beyond simulation event.")

    mᵤ.D ∈ ["left", "right"] || error("Direction not specified.")
    mₑ.D ∈ ["left", "right"] || error("Direction not specified.")


    D1 = mᵤ.D == "left" ? -1 : 1
    D2 = mₑ.D == "left" ? -1 : 1

    # motion amplitude is noisy if stochastic simulation
    A1 = stochastic == 0 ? mᵤ.A : rand(Normal(mᵤ.A, 2.0))
    A2 = stochastic == 0 ? mₑ.A : rand(Normal(mₑ.A, 2.0))

    noise = sensor.noise
    f1 = 1/(duration1) # frequency: single sinusiodal head turn
    f2 = 1/(duration2)

    timesteps = range(0, stop = total_duration, step = Δt)
    T = length(timesteps)

    κ1, κ2 = cupula_dynamics(Δt)

    α1 = zeros(T)
    α2 = zeros(T)
    α = zeros(T)
    ω1 = zeros(T)
    ω2 = zeros(T)
    ω = zeros(T) # ω = Array{Float64}(undef, T) # velocity
    θ = zeros(T)
    c = zeros(T)
    y = zeros(T)

    for i ∈ Iterators.drop(1:T, 1)

        t = round(timesteps[i]; digits = 3)

        α1[i] = (t > onsetᵤ) & (t < onsetᵤ + duration1) ? acceleration(D1, A1, f1, t, onsetᵤ, noisy = stochastic) : 0
        α2[i] = (t > onsetₑ) & (t < onsetₑ + duration2) ? acceleration(D2, A2, f2, t, onsetₑ, noisy = stochastic) : 0
        α[i] = α1[i] + α2[i]
        ω1[i] = ω1[i-1] + Δt * α1[i]
        ω2[i] = ω2[i-1] + Δt * α2[i]

        ω[i] =  ω1[i] + ω2[i]
        θ[i] = θ[i-1] + Δt * ω[i] + 1/2 * Δt^2 * α[i]
        c[i] = κ1 * c[i-1] + κ2 * ω[i]
        y[i] = ω[i] - c[i] + rand(noise)
    end
    
    out = (timesteps = collect(timesteps),
            y = y,
            αᵤ = α1, αₑ = α2, 
            α = α, 
            ωᵤ = ω1, ωₑ = ω2,
            ω = ω, 
            θ = θ, c = c, 
            Δt = Δt, onsetᵤ = onsetᵤ, onsetₑ = onsetₑ,
            sensor = noise, stochastic = stochastic)
    return out

end

simulate(mᵤ::HeadTurn,
    sensor::Sensor;
    Δt = 0.01,
    duration::Real = 5, 
    stochastic::Bool = false) = simulate(mᵤ, ExternalDisturbance(A = 0, duration = 0.1), sensor, 
        Δt = Δt, duration = duration, stochastic = stochastic)

simulate(mₑ::ExternalDisturbance,
    sensor::Sensor;
    Δt = 0.01,
    duration::Real = 5, 
    stochastic::Bool = false) = simulate(HeadTurn(A = 0, duration = 0.1), mₑ, sensor, 
        Δt = Δt, duration = duration, stochastic = stochastic)




