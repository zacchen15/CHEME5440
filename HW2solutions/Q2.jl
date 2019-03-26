# Script for Q2 PS2-S19

# load the GRNSimKit and PyPlot packages -
using GRNSimKit
using PyPlot

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
dd = build_default_data_dictionary(path_to_model_file)

# Run the model to steady-state, before we do anything -
steady_state = GRNSteadyStateSolve(dd)

# Run the model 1 time unit *before* we add inducer -
dd[:initial_condition_array] = steady_state
(T0, X0) = GRNDynamicSolve((0.0,1.0), dd)

# Add inducer T = 1 mM for 100 time units -
dd[:initial_condition_array] = X0[end,:]
dd[:initial_condition_array][7] = 10.0 # mM
tstart_1 = T0[end]
tstop_1 = tstart_1 + 1.0
(T1, X1) = GRNDynamicSolve((tstart_1,tstop_1), dd)

# Wash inducer out -
dd[:initial_condition_array] = X1[end,:]
dd[:initial_condition_array][7] = 0.0
tstart_2 = T1[end]
tstop_2 = tstart_2 + 10.0
(T2, X2) = GRNDynamicSolve((tstart_2,tstop_2), dd)

# Package -
T = [T0 ; T1; T2]
X = [X0 ; X1; X2]

# plots -
if simulation_case_flag == :default

    # make a plot -
    plot(T*(60),X[:,4],color_P1,linewidth=2, label="P1")
    plot(T*(60),X[:,5],color_P2,linewidth=2, label="P2")
    plot(T*(60),X[:,6],color_P3,linewidth=2, label="P3")

    legend()

elseif simulation_case_flag == :broken

    # plot w/dash for broken
    plot(T*(60),X[:,4],color_P1, linewidth=2, linestyle="--")
    plot(T*(60),X[:,5],color_P2, linewidth=2, linestyle="--")
    plot(T*(60),X[:,6],color_P3, linewidth=2, linestyle="--")

end

# axis -
xlabel("Time (min)", fontsize=16)
ylabel("Protein (nmol/gDW)", fontsize=16)
