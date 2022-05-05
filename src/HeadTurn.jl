module HeadTurn

using Parameters, Turing, Random

export NigmatullinaTrial, HeadTurn, PlannedMovement, ExternalDisturbance, Sensor
export simulate
export cupula_dynamics, acceleration, velocity, displacement, amplitude, peak_velocity


include("selfmotion.jl")
include("turingutils.jl")
include("turingmodels.jl"

end
