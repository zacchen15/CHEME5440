#CHEME5440 HW1
###Zac Chen

##Part A

### Derivation of Transcription Rate Constant
The mass balances around the open and closed complex are:
$$ \frac{d}{dt}(G_j:R_x)_{C} = k_{+}(G_j)(R_x) - k_{-}(G_j:R_x)_C - k_I(G_j:R_x)_C $$
$$ \frac{d}{dt}(G_j:R_x)_O = k_I(G_j:R_x)_C -k_A(G_j:R_x)_O - k_{E,j}(G_j:R_x)_{O} $$

The total amount of RNA polymerase is:
$$ R_{X,T} = R_X + (G_j:R_x)_{C} + (G_j:R_x)_{O} $$

Assuming steady-state, the amount of open and closed complexes are:
$$ (G_j:R_x)_{C} = (\frac{k_+}{k_-+k_I})(G_j)(R_x) $$
$$ (G_j:R_x)_{O} = (\frac{k_I}{k_A+k_{E}})(G_j:R_x)_{c} $$

Let the ratio of parameters be:
$$ K_{X,j} = \frac{k_{-}+k_{I}}{k_{+}} $$
$$ \tau_{X,j} = \frac{k_A+k_E}{k_{I}} $$

Replacing for terms in the expression for open complex and for the total amount of RNA polymerase, we get:
$$ (G_j:R_x)_{O} = \frac{(G_j)(R_x)}{(K_{X,j})(\tau_{X,j})} $$
$$ R_{X,T} = R_X + (G_j)(R_x)(K_{X,j}^{-1}) + (G_j)(R_x)(K_{X,j}^{-1})(\tau_{X,j}^{-1}) $$

Rearranging, we get:
$$ R_{X} = \frac{R_{X,T}(K_{X,j})(\tau_{X,j})}{(K_{X,j})(\tau_{X,j})+(\tau_{X,j}+1)G_j} $$
$$ (G_j:R_x)_{O} = \frac{R_{X,T}G_j}{(K_{X,j})(\tau_{X,j})+(\tau_{X,j}+1)G_j} $$

Finally, we get an equation for the kinetic rate of transcription in terms of total RNA polymerase and other solvable parameters
$$ r_{X,j} = k_{E,j}R_{X,T}(\frac{R_{X,T}G_j}{(K_{X,j})(\tau_{X,j})+(\tau_{X,j}+1)G_j}) $$

### List of Parameters Used from BioNumbers
|    Parameter    |  Value  | Units |             Source             |
|:---------------:|:-------:|:-----:|:------------------------------:|
|      $e_x$      |   27    | nt/s  |   Proshkin S, et al. (2010)    |
|       $L$       |   924   |  nt   |      Xu L, et al. (2006)       |
|     $K_xj$      | 2.48e-8 |   M   |    McClure, et al. (1980)*     |
|      $k_I$      |  0.04   |  1/s  |     McClure, et al. (1980)     |
|     $R_X,T$     |   30    |  nM   |     Arkin A, et al. (1998)     |
| $Radius_{cell}$ |  0.306  |  um   | Rosenberger RF , et al. (1978) |
| $Length_{cell}$ |  2.62   |  um   | Rosenberger RF , et al. (1978) |
\* Note: This parameter was calculated based on the intercept and slope obtained from Figure 2 in McClure, et al. (1980).

## Part B

The constant $\tau_{X,j}$ is time constant which compares the initiation, elongation, and abortive constant as shown below.
$$ \tau_{X,j} = \frac{k_A+k_E}{k_{I}} $$

Assuming that the abortive constant is negligibly small, then the time constant is a ratio between the elongation constant and initiation constant:
$$ \tau_{X,j} = \frac{k_E}{k_{I}} $$

Using parameters determined from BioNumbers, this was determined to be:
$$ \tau_{X,j} = 0.220 $$

Since $\tau$ is less than one, transcription is elongation limited.


## Part C

The following mass balance for the mRNA concentration was given:
$$ \frac{dm_j}{dt} = r_{x,j}u_j(I) - ( k_{x,j}^d - \mu )m_j$$

where the utility function $u_j$ is a function of the inducer concentration:
$$ u_j = \frac{W_1+W_2f_I}{1+W_1+W_2f_I} $$
$$ f_I = \frac{I^n}{K+I^n} $$

Using the given parameters and the steady state assumption for mRNA concentrations, the following equation was derived and plotted against inducer concentration:
$$ m_j = \frac{r_{x,j}u_j}{k_{x,j}^d+\mu} $$

Please see HW1.jl file for actual plot. 
