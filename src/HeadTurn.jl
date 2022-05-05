module HeadTurn

using Parameters, Turing, Random

export NigmatullinaTrial, PlannedHeadTurn, ExternalDisturbance, Sensor
export simulate
export cupula_dynamics, acceleration, velocity, displacement, amplitude, peak_velocity


include("selfmotion.jl")
include("turingmodels.jl")

end
