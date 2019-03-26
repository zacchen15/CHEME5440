# Sample script to run a model using GRNSimKit -
using GRNSimKit
using PyPlot

# time step -
time_step_size = 1.0*(1/60)

# what case are we doing? -> {:default, :broken}
simulation_case_flag = :default

# default -
path_to_model_file = "$(pwd())/Memory.json"
color_P1 = "red"
color_P2 = "green"
color_P3 = "blue"
if simulation_case_flag == :broken

    # point to broken file -
    path_to_model_file = "$(pwd())/Memory-Broken.json"

    # point to colors -
    color_P1 = "lightcoral"
    color_P2 = "lightgreen"
    color_P3 = "powderblue"
end

# Build a data dictionary from a model file -
ddd = build_discrete_dynamic_data_dictionary(time_step_size, path_to_model_file)

# Run the model to steady-state, before we do anything -
steady_state = GRNSteadyStateSolve(ddd)

# Run the model 100 time units *before* we add inducer -
ddd[:initial_condition_array] = steady_state
(T0, X0) = GRNDiscreteDynamicSolve((0.0,1.0,time_step_size), ddd)

# Add inducer T = 1 mM for 100 time units -
ddd[:initial_condition_array] = X0[end,:]
ddd[:initial_condition_array][7] = 10.0
tstart_1 = T0[end]
tstop_1 = tstart_1 + 1.0
(T1, X1) = GRNDiscreteDynamicSolve((tstart_1,tstop_1, time_step_size), ddd)

# Wash inducer out -
ddd[:initial_condition_array] = X1[end,:]
ddd[:initial_condition_array][7] = 0.0
tstart_2 = T1[end]
tstop_2 = tstart_2 + 10.0
(T2, X2) = GRNDiscreteDynamicSolve((tstart_2,tstop_2, time_step_size), ddd)

# Package -
T = [T0 ; T1; T2]
X = [X0 ; X1; X2]

# make a plot -
plot(T*(60),X[:,4],color_P1,linewidth=2, linestyle="--")
plot(T*(60),X[:,5],color_P2,linewidth=2, linestyle="--")
plot(T*(60),X[:,6],color_P3,linewidth=2, linestyle="--")

# axis -
xlabel("Time (min)", fontsize=16)
ylabel("Protein (nmol/gDW)", fontsize=16)
