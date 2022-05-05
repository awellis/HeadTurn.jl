module HeadTurn

using Parameters, Turing, Random

export NigmatullinaTrial, PlannedHeadTurn, ExternalDisturbance, Sensor
export simulate
export cupula_dynamics, acceleration, velocity, displacement, amplitude, peak_velocity
export compute_stats, get_omega, summarize_omega

include("selfmotion.jl")
include("turingmodels.jl")
include("turingutils.jl")

end
