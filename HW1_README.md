#CHEME5440 HW1
###Zac Chen

##Part A

The mass balances around the open and closed complex are:
$$ \frac{d}{dt}(G_j:R_x)_{C} = k_{+}(G_j)(R_x) - k_-(G_j:R_x)_C - k(G_j:R_x)_C $$
$$ \frac{d}{dt}(G_j:R_x)_O = k(G_j:R_x)_C -k_A(G_j:R_x)_O - k_{E,j}(G_j:R_x)_{O} $$

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
$$ r_{X,j} = k_{E,j}R_{X,T}(\frac{R_{X,T}G_j}{(K_{X,j})(\tau_{X,j})+(\tau_{X,j}+1)G_j})
