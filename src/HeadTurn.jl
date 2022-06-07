module HeadTurn

using Parameters, Turing, Random
using CairoMakie

export NigmatullinaTrial, PlannedHeadTurn, ExternalDisturbance, Sensor
export simulate
export cupula_dynamics, acceleration, velocity, displacement, amplitude, peak_velocity
export makeheadturnmodel, headturnmodel
export q025, q975, compute_stats, get_omega, summarize_omega, get_estimates, plot_estimate
export savefig, publication_theme

include("utils.jl")
include("selfmotion.jl")
include("turingmodels.jl")
include("turingutils.jl")

end
