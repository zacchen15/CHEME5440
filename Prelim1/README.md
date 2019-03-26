# Prelim 1
## Question 1
### Part A
Please see ```Q1.jpg``` for written solution
### Part B
Please see ```Q1.jpg``` for written solution
## Question 2
### Part A
Use ```julia> include("Q2.jl")``` to run the file for Question 2.
Plot of mRNA concentration vs time is in Figure 1. Plot of protein concentration vs time is in Figure 2.
### Part B
To view the sensitivity array, sp<sub>ij</sub>, for each protein, p<sub>i</sub>, where i = {1,2,3} during Phase j, where j = {1,2,3} use
```julia> spij```

Note that the Phases are Phase I, Early Phase II, and Late Phase II, respectively.
### Part C
The **U** matrix can be called similarly to the sensitivity matrix using
```julia> Upij```

## Question 3
Note: Before running ```Q3.jl```, restart the Julia REPL. This file uses the variable, I, which will conflict with the LinearAlgebra package in ```Q2.jl``` if ```Q3.jl```is run after ```Q2.jl```
### Part A
To view stoichiometric matrix use
```julia> include("Q3.jl")```
```julia> print(S)```

Expressions for the rate of transcription and translation can be seen in the associated files ```RateTranscription.jl``` and ```RateTranslation.jl```

### Part B
Running ```Q3.jl``` will also create a figure for the protein concentration vs inducer concentration.

### Part C
To view the dual value/shadow price array at an inducer concentration use
```julia> print(dual_value_array[:,i])```
where i/10000 is the inducer concentration of interest. Ex: to look at I = 10 mM, use i = 100000
