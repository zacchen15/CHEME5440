## CHEME 5440 Problem Set 3
### Zac Chen
#### Part 1A
Using the Kyoto Encyclopedia of Genes and Genomes (KEGG) database, a stoichiometric matrix, **S**, was derived for arginine biosynthesis in human cells. Based on the boundaries cycle given, there are a total of 6 reaction fluxes, **v<sub>j</sub>**, and 12 transport fluxes, **b<sub>j</sub>**, which create the 18 rows of the stoichiometric matrix. Based on the each of the 6 reactions (shown below), there are a total of 20 metabolites that make up the 20 columns of the stoichiometric matrix.

*<center> **v<sub>1</sub>**: ATP + citrulline + aspartate &rightarrow; AMP + pyrophosphate + argininosccuinate</center>*

*<center> **v<sub>2</sub>**: argininosccuinate &rightarrow; fumarate + arginine </center>*

*<center> **v<sub>3</sub>**: arginine + H<sub>2</sub>O &rightarrow; ornithine + urea </center>*

*<center> **v<sub>4</sub>**: carbamoyl-phosphate + ornithine &rightarrow; phosphate + citrulline </center>*

*<center> **v<sub>5f</sub>/v<sub>5r</sub>**: arginine + 1.5 NADPH  + 1.5 H<sup>+</sup> + O<sub>2</sub> &leftarrow; &rightarrow;  citrulline + NO + 1.5 NADP<sup>+</sup> + 2 H<sub>2</sub>O </center>*

Transport fluxes were assigned a value of 1 if the metabolite was a reactant that needed to be brought into the system or -1 if they were a product that needed to be removed from the system. For example, ATP is a consumed metabolite in the system. In order to balance the consumption at reaction v<sub>1</sub>, the transport flux for ATP, b<sub>5</sub>, was set to 1.

Due to diffculties in typesetting matrices with more than 10 columns, it is suggested to view to the stoichiometric matrix in the command window after running the file:
```
julia> include("HW3.jl")
julia> S
```
<!--
$$
\left(\begin{array}{cccccccccccccccccccc}
-1 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0\\
1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\\
0 1 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0\\
0 1 -1 0 -1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0\\
0 0 1 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0\\
0 0 1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\\
0 0 0 -1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0\\
-1 0 0 1 1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0\\
-1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0\\
0 0 -1 0 2 -2 0 0 0 0 0 1 0 0 0 0 0 0 0 0\\
0 0 0 1 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0\\
1 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0\\
1 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0\\
0 0 0 0 1 -1 0 0 0 0 0 0 0 0 0 -1 0 0 0 0\\
0 0 0 0 -2 2 0 0 0 0 0 0 0 0 0 0 1 0 0 0\\
0 0 0 0 -1.5 1.5 0 0 0 0 0 0 0 0 0 0 0 1 0 0\\
0 0 0 0 -1.5 1.5 0 0 0 0 0 0 0 0 0 0 0 0 1 0\\
0 0 0 0 1.5 -1.5 0 0 0 0 0 0 0 0 0 0 0 0 0 -1\\
\end{array}\right)
\left(\begin{array}{cc}
v_1\\
v_2\\
v_3\\
v_4\\
v_{5f}\\
v_{5r}\\
b_1\\
b_2\\
b_3\\
b_4\\
b_5\\
b_6\\
b_7\\
b_8\\
b_9\\
b_{10}\\
b_{11}\\
b_{12}\\
b_{13}\\
b_{14}\\
\end{array}\right)
$$
-->

#### Part 1B
Next, a 6x18 element matrix, **A**, was created to determine if the model urea cycle was elementally balanced in C, H, N, O, P, and S. In this matrix, each row represented an element, while each column represented a metabolite. To determine if the cycle was balanced, the element matrix was multiplied by the stoichiometric matrix (18x20) to give a 6x20 elemental balance matrix, **E**, where the rows represented each element, and the columns represented each reaction.

$$
E = A*S
$$

In this matrix, each position represents the balance of each element for each flux. For example, a 1 in the row for H and column for v<sub>1</sub> represents an unbalanced hydrogen gained in the v<sub>1</sub> reaction. In the model, created, the product of these showed zeros for the balance at all the reactions, v<sub>j</sub>, whereas all the transport fluxes, b<sub>j</sub> returned values for C, H, N, O, P and S representative of the molecule transported. The element matrix and elemental balance matrix can be viewed using:
```
julia> include("HW3.jl")
julia> A
julia> E
```

#### Part 1C
The maximum urea production from the cycle was then determined using the Flux Balance Analysis file with the assumption that (1) the reaction fluxes, v<sub>n</sub>, were known, (2) 0 &le; b<sub>j</sub> &le; 10 mmol/gDW-hr, and (3) steady-state enzyme concentration, **E** = 0.01 umol/gDW.

The value of the reaction fluxes were determined using the equation:
$$
v_j = k_{cat, j}E\prod_{i=1}^{n}{\frac{X_i}{K_{M,i}+X_i}}
$$
where k<sub>cat</sub> is the specific rate of each enzyme (given), E is the steady state enzyme concentration (given), X<sub>i</sub> is the steady state metabolite concentration, and K<sub>M,i</sub> is the saturation constant for that metabolite. Values were determined using KEGG, Bionumbers, and Park, et al (2016). In cases that values could not be found, the assumption that X<sub>i</sub> >> K<sub>M,i</sub> was made so that the saturation term would be ~1.

Using this set of equations and constraints, the maximum urea concentration was determined to be 0.00278 mmol/gDW-s.
