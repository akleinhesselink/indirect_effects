# Essential Resource Competition Demo 
# R-code supporting 'Indirect effects of environmental change in resource competition models' 
# Written on R version 3.2.1 (2015-06-18) 
# Author: Andrew R. Kleinhesselink
# Contact: arklein@aggiemail.usu.edu 
# Requires loading the essential_resource_functions.R for function definitions 
# Contents: 
#   Section 1: demonstrates two species essential resource model 
#   Section 2: demonstrates equilibrium change with change in resource supply point
#   Section 3: compares analytical and simulated indirect effects 

rm(list = ls())
require( deSolve ) # Package for performing numerical solutions to differential equations
source('essential_resource_functions.R') # load functions bundled with this script

# Define parameters -----------------------------------------------------------

N = c(0.01, 0.01)  # initial populations/biomasses 
R = c(20, 20) # initial resource concentrations 
NR <- c(N = N, R = R)

r = c(1, 1) # growth rates
k = matrix(c(35, 10, 10, 35), nrow = 2, byrow = T) # resource half-saturation for species one and two
m = c(0.21, 0.2) # mortality rates -- giving species two a slight advantage
a = c(0.2, 0.2)  # resource growth rates
Q = matrix(c(3, 1, 2, 4), nrow = 2, byrow = T) # species one and two per capita consumption of resources
S = as.numeric(NR[3:4]) # resource supply points 

parms <- list(r = r, k = k, m = m, a = a, Q = Q, S = S) 
# ------------------------------------------------------------------------------

# Section 1: Demonstrate two species essential resource competition and coexistence
Rstars <- calc_rstars(parms = parms)
Rstars

Rstars[1, 1] > Rstars[1, 2]  #### check that species one is limited by R1

beta <- calc_essential_consumption(parms = parms) # slopes of consumption vectors 
beta 

rho <- calc_essential_rho(parms = parms)  # calculate rho directly from parameters
rho

t <- 1:1000
RRout <- ode(y = NR, times = t, func = essential_mod, parms = parms)  # run model 

# plot simulation to equilibrium
ylab1 = "Population/Biomass"
xlab1 = "Time"
ylab2 = "Resources"

par(mfrow = c(1, 2), pty = "s", tcl = -0.2, mgp = c(2.2, 0.5, 0), mar = c(4, 5, 1, 1), 
    las = "1", cex.lab = 1.2, family = "")

clrs <- c("green", "blue")

matplot(RRout[, 1], RRout[, 2:3], type = "l", xlab = xlab1, ylab = ylab1, lty = 1:2, col = clrs)
legend("topright", c("Sp. 1", "Sp. 2"), lty = 1:2, bty = "n", col = clrs, cex = 1.2)

matplot(RRout[, 1], RRout[, 4:5], type = "l", xlab = xlab1, ylab = ylab1)
legend("topright", c("R1", "R2"), lty = 1:3, bty = "n", cex = 1.2)

par(mfrow = c(1, 1), pty = "s", tcl = -0.2, mgp = c(2.2, 0.5, 0), mar = c(4, 5, 1, 1), las = "1", cex.lab = 1.2)

plot_essential_zngis(parms, RRout)  # Resource isocline plot compare figure 2 in main text 
text(35, 25, bquote(rho ~ "=" ~ .(rho)))

# Section 2: Demonstrate how change in resource supply point effects focal species (N1) 

newT <- 5000
S1delta <- 1 # change resource 1 supply point 
S2delta <- 0 # hold resource 2 supply point constant

out <- simulate_essential_effects(S1delta, S2delta, NR, parms, t = 1:newT)

par(mfrow = c(1, 3), pty = "s", tcl = -0.2, mgp = c(2.2, 0.5, 0), mar = c(3, 5, 2, 1), 
    oma = c(3, 1, 3, 1), las = "1", cex.lab = 1.2)
plot_ifx_panel(out[[1]], out[[2]], out[[3]], t = newT)  ## plot perturbation in S1

# resource 2 sensitivity
S1delta2 <- 0 # hold resource 1 supply point constant
S2delta2 <- 1 # change resource 2 supply point 

out2 <- simulate_essential_effects(S1delta2, S2delta2, NR, parms, t = 1:newT)
plot_ifx_panel(out2[[1]], out2[[2]], out2[[3]], t = newT)  ## plot perturbation in S2


# Section 3: Check whether the simulation matches algebra for defining direct and indirect effects
t <- 1:1000

NR1 <- ode(y = NR, times = t, func = essential_mod, parms = parms)  ### calcs N1, N2, R1, R2
preEq <- as.numeric(tail(NR1[, 2:5], 1))  # first equilibrium 
names(preEq) <- names(NR)

newParms <- parms
newParms$S <- parms$S * c(1.01, 1)  # change S1 slightly
deltaS1 <- newParms$S[1] - parms$S[1]

NR2 <- ode(y = preEq, times = t, func = essential_mod, parms = newParms)
Eq2 <- as.numeric(tail(NR2[, 2:5], 1))  # net effects equilibrium 
names(Eq2) <- names(NR)

NR2direct <- ode(y = preEq, times = t, func = essential_mod_hold_comp, parms = newParms)
Eq2direct <- as.numeric(tail(NR2direct[, 2:5], 1))  # direct effects equilibrium 
names(Eq2direct) <- names(NR)

indFX <- Eq2[1] - Eq2direct[1]  # change due to indirect effects 
indFX

dirFX <- Eq2direct[1] - preEq[1]  # change due to direct effects
dirFX

# direct sensitivity (change in population/change in resource)
simDirSens <- dirFX/(deltaS1)  
simDirSens

simIndSens <- indFX/(deltaS1)  # Indirect sensitivity from simulation
simIndSens

# Calculates direct and indirect sensitivity analytically
rho <- calc_essential_rho(parms)

dirFXN1_calc <- with(parms, (a[1])/(m[1] * Q[1, 1])) # see equation 12 in main text
indFXN1_calc <- with(parms, ((a[1]/(m[1] * Q[1, 1])) * (rho^2/(1 - rho^2)))) # see equation 13 in main text

cbind(dirFXN1_calc, simDirSens)  # check that calculated and simulated values match
cbind(indFXN1_calc, simIndSens)  # check that calculated and simulated values match

# Effects of a change in the supply of the non-limiting resource
t <- 1:1000
NR1 <- ode(y = NR, times = t, func = essential_mod, parms = parms)  ### calcs N1, N2, R1, R2
preEq <- as.numeric(tail(NR1[, 2:5], 1))  # first equilibrium 
names(preEq) <- names(NR)

newParms <- parms
newParms$S <- parms$S * c(1, 1.01)  # change S2 slightly
deltaS2 <- newParms$S[2] - parms$S[2]

NR2 <- ode(y = preEq, times = t, func = essential_mod, parms = newParms)
Eq2 <- as.numeric(tail(NR2[, 2:5], 1))  # net effects equilibrium 
names(Eq2) <- names(NR)

indFX2 <- Eq2[1] - preEq[1]  # Indirect change is equal to net effects 
indFX2

simIndSens2 <- indFX2/(deltaS2)  # Indirect sensitivity (change in population/change in resource)
simIndSens2

# This section calculates direct and indirect effect analytically 
indFX2_calc <- with(parms, -((a[1]/(m[1] * Q[1, 2])) * (rho^2/(1 - rho^2)))) # see equation 15 in main text

# check that calculated and simulated values match
cbind(indFX2_calc, simIndSens2)

