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


##Part 2
###Part 2a
Mass Balances for mRNA, $m_i$ and protein, $p_i$, were derived based on transcription/translation rate, degradation rate, and dilution rate.
For the mRNA balance:
$ \frac{dm_i}{dt} = r_{x,i}u_{x,i} -  (k_{x,d}-\mu)m_i $
For the protein balance:
$ \frac{dp_i}{dt} = r_{L,i}u_{L,i} -  (k_{L,d}-\mu)p_i $
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


###Part 2b
####Functional Memory Circuit
####Broken Memory Circuit
###Part 2c
