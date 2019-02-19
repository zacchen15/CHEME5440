#CHEME5440 HW1
###Zac Chen

##Part 1
Using the model developed in Problem Set 1 the parameters in Table 1 were increased and decreased by 50%. Their observed effect was recorded below. Most of the parameters that affected mRNA concentration were part of the control function, $u$.

Table 1: Parameters changed in the transcription model and their observed effect on the graph of the steady state mRNA concentration
| Parameter |                Definition                |                                                                              Effects on mRNA Concentration                                                                               |
|:---------:|:----------------------------------------:|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|
|   $e_x$   |         Average Elongation Rate          |                                                                                        No Visible  Effect                                                                                         |
|    $L$    |           Average Gene Length            |                                                                                        No Visible Effect                                                                                         |
|  $K_xj$   |          Saturation Coefficient          |                                                                                        No Visible Effect                                                                                         |
|   $k_I$   |     Rate for Open Complex Formation      |                                                                                        No Visible Effect                                                                                         |
|  $R_X,T$  |          Concentration of RNAP           |                                                                                        No Visible Effect                                                                                         |
|   $K_1$   | Binding Constant of RNAP without Inducer | Increase in parameter makes basal transcription rate increase at low inducer concentration. Decrease in parameter makes basal transcription rate decrease at low inducer concentration.  |
|   $K_2$   |  Binding Constant of RNAP with Inducer   | Increase in parameter makes the slope of the curve steeper in the middle part of the diagram. Decrease in parameters decreases the slope of the curve in the middle part of the diagram. |
|    $n$    |          Cooperativity Constant          | Increase in parameter makes the slope of the curve steeper in the middle part of the diagram. Decrease in parameters decreases the slope of the curve in the middle part of the diagram. |
|   $K_d$   |       Binding Function for Inducer       |                                                                                    No Visible Effect                                                                                     |

![](HW2Q1_Original.png)
Figure 1: Concentration of steady state mRNA over a range of inducer concentrations.

##Part 2
###Part 2a
Mass Balances for mRNA, $m_i$ and protein, $p_i$, were derived based on transcription/translation rate, degradation rate, and dilution rate.
For the mRNA balance and protein balances:
$$ \begin{aligned}
\frac{dm_i}{dt} = r_{x,i}u_{x,i} -  (k_{x,d}-\mu)m_i \\
\frac{dp_i}{dt} = r_{L,i}u_{L,i} -  (k_{L,d}-\mu)p_i
\end{aligned} $$
Note that since the first protein is outside the cell, there is no dilution term for $m_1$ and $p_1$. The system of 6 ODEs can be simplified as
$ \frac{dx}{dt} = Ax -  Sr $
where
$$ A =
\left(\begin{array}{cccccc}
-k_{x,d} & 0 & 0 & 0 & 0 & 0\\
0 & -k_{x,d}-\mu & 0 & 0 & 0 & 0\\
0 & 0 & -k_{x,d}-\mu & 0 & 0 & 0\\
0 & 0 & 0 & -k_{L,d} & 0 & 0\\
0 & 0 & 0 & 0 & k_{L,d}-\mu & 0\\
0 & 0 & 0 & 0 & 0 & k_{L,d}-\mu\\
\end{array}\right)
$$

$$ x =
\left(\begin{array}{cccccc}
m_1\\
m_2\\
m_3\\
p_1\\
p_2\\
p_3\\
\end{array}\right)
$$

$$ S =
\left(\begin{array}{cccccc}
1 & 0 & 0 & 0 & 0 & 0\\
0 & 1 & 0 & 0 & 0 & 0\\
0 & 0 & 1 & 0 & 0 & 0\\
0 & 0 & 0 & 1 & 0 & 0\\
0 & 0 & 0 & 0 & 1 & 0\\
0 & 0 & 0 & 0 & 0 & 1\\
\end{array}\right)
$$

$$ x =
\left(\begin{array}{cccccc}
r_{x,1}u_{x,1}\\
r_{x,2}u_{x,2}\\
r_{x,3}u_{x,3}\\
r_{L,1}u_{L,1}\\
r_{L,2}u_{L,2}\\
r_{L,3}u_{L,3}\\
\end{array}\right)
$$

Where $r_i$ are the kinetic transcription and translation rates and $u_i$ are control functions for transcription and translation. The control function, $u_i$, is given as:

$$
u_i = \frac{K_x+\sum{K_{ij}f_{i,j}}}{1+K_x+\sum{K_{ij}f_{i,j}}}
$$

###Part 2b
Note: to run the code for ```HW2.jl``` the functions ```RateTranslation.jl``` and ```RateTransription.jl``` must be initialized by typing
```
include("RateTranslation.jl")
include("RateTransription.jl")
```
For all simulations of the memory circuit, the inducer concentration was set to 10 mM for 60 minutes, and set to zero after 50 minutes. These changes are reflected in the sharp drop in the mRNA charts at the 1 hour mark.

####Functional Memory Circuit
In the functional memory circuit after the inducer is removed, the mRNA concentration for protein 1 drops dramatically, however, the mRNA concentrations remains at steady state as shown in Figure 3. Due to the mutual induction between protein 2 and 3, once they are inducer by protein 1, they remain in the cell as shown in Figure 4.

![](HW1Q2b_MemoryCircuitmRNA.png)
Figure 3: mRNA concentration as a function of time in the functional memory circuit. $m_1$ is in black, $m_2$ is in blue, $m_3$ is in red.
![](HW1Q2b_MemoryCircuitProtein.png)
Figure 4: Protein concentration as a function of time in the functional memory circuit. $p_1$ is in black, $p_2$ is in blue, $p_3$ is in red.



####Broken Memory Circuit
In the broken memory circuit, protein 2 does not induce protein 3 production. In this circuit, after the inducer is removed, the protein 1 concentration decreases, followed by the protein 3 concentration. protein 2 remains longest of the three because protein 3 continues to induce the production of protein 2. In this circuit, the induction of protein 3 by protein 2 is removed by setting
$$ f_{2,3} = \frac{(p_2)^n}{(K_d)^n + (p_2)^n} = 0 $$

![](HW1Q2b_BrokenCircuitmRNA.png)
Figure 5: mRNA concentration as a function of time in the broken memory circuit. $m_1$ is in black, $m_2$ is in blue, $m_3$ is in red.
![](HW1Q2b_BrokenCircuitProtein.png)
Figure 6: Protein concentration as a function of time in the broken memory circuit. $p_1$ is in black, $p_2$ is in blue, $p_3$ is in red.

###Part 2c
In Part 2c, a discretization was used to approximate the solution achieved by the ODE solver in part 2b. The approximation is as follows:
$$ x_{k+1} = A_{hat}x_k + S_{hat}r_k $$
where
$$  A_{hat} = exp(A*\tau) $$
$$  S_{hat} = A^{-1}[A_{hat}-I]S $$

Using this approximation, a solution of the same order of magnitude was determined as shown in the mRNA concentration and protein concentration shown in Figure 7 and 8.

![](HW1Q2c_MemoryCircuitmRNA.png)
Figure 7: mRNA concentration as a function of time in the functional memory circuit. $m_1$ is in black, $m_2$ is in blue, $m_3$ is in red.
![](HW1Q2c_MemoryCircuitProtein.png)
Figure 8: Protein concentration as a function of time in the functional memory circuit. $p_1$ is in black, $p_2$ is in blue, $p_3$ is in red.
