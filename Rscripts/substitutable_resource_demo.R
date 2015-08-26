# Substitutable Resource Competition Demo 
# R-code supporting 'Indirect effects of environmental change in resource competition models' 
# Written on R version 3.2.1 (2015-06-18) 
# Author: Andrew R. Kleinhesselink
# Contact: arklein@aggiemail.usu.edu 
# Requires loading the substitutable_resource_functions.R for function definitions 
# Contents: 
#   Section 1: demonstrates two species substitutable resource model 
#   Section 2: demonstrates equilibrium change with change in resource supply point
#   Section 3: compares analytical and simulated indirect effects 

rm(list = ls())
require( deSolve ) # Package for performing numerical solutions to differential equations
source('substitutable_resource_functions.R') # load functions bundled with this script

# Define parameters -----------------------------------------------------------
N = c(1, 1) # initial population sizes
R = c(4000, 4000)  # initial resource levels

NR <- c(N = N, R = R)  # initial populations and resource levels

S <- as.numeric( c(1000, 800) ) # resource supply points 
r = c(0.81, 0.8) # growth rates
k = c(152, 150)  # resource half-saturation for species one and two
tau = c(0.61, 0.52)  # minimum resource requirement for each species -- giving species two a slight advantage
D = 0.734 # mortality rates
Q = matrix(c(2.01, 5.04, 5.1, 2.234), nrow = 2, byrow = T) # species one and two per capita consumption of resource one and two
W = matrix(c(1.12, 3.0145, 3.111, 1.091), nrow = 2, byrow = T) # value of resource one and two to each species

parms <- list(r = r, k = k, tau = tau, D = D, Q = Q, W = W, S = S)
# ------------------------------------------------------------------------------


# Section 1: Demonstrate two species essential resource competition and coexistence

B <- calc_intercepts(parms)  # Calculate Zero Net Growth Isocline y-intercepts

# Plot Zero-net growth isoclines and consumption vectors and supply point
# compare to figure A1 in Appendix A
par(mfrow = c(1,1))
plot_substitutable_zngis(parms, B1 = B[[1]], B2 = B[[2]]) 
Rstars <- calc_substitutable_rstars(parms, B1 = B[[1]], B2 = B[[2]])
plot_substitutable_consumption(parms, B1 = B[[1]], B2 = B[[2]], R1star = Rstars[[1]], R2star = Rstars[[2]])
points(parms$S[1], parms$S[2])

# Run model
NRout <- ode(y = NR, func = substitutable_mod, t = 1:2000, parms = parms, method = "rk4")
preEq <- as.matrix(tail(NRout, 1))[, -1]

# show that calculated Rstars match equilibrium resources in simulation:
cbind( R1 = Rstars[[1]], R2 = Rstars[[2]])
preEq[3:4]


# Section 2: Demonstrate how change in resource supply point effects focal species (N1) 
newParms <- parms
newParms$S <- parms$S * c(1.1, 1)  # increases supply point for resource one

# Re-run model to equilibrium starting at old equilibrium
NRout2 <- ode(y = preEq, func = substitutable_mod, t = 1:2000, parms = newParms, method = "rk4")
Eq2 <- as.matrix(tail(NRout2, 1))[, -1]
Eq2

# find net change in N1
Eq2[1] - preEq[1]

# Sensitivity of N1 to change in S1 (partial derivative)
simNetSensitivity <- (Eq2[1] - preEq[1])/(newParms$S[1] - parms$S[1])
simNetSensitivity

# Re-run model with constant competitor abundance starting from first equilibrium
NRoutDirect <- ode(y = preEq, func = substitutable_mod_hold_comp, t = 1:2000, parms = newParms, method = "rk4")
dirEq <- as.matrix(tail(NRoutDirect, 1))[, -1]
dirEq[1] - preEq[1]  # change in N1 due to direct effects

points(dirEq[3], dirEq[4], pch = 20)  # show that R1*, R2* shifts a tiny bit when competitor is held constant
text(dirEq[3], dirEq[4], "new Eq", cex = 0.8, pos = 4)

# plot population time series
par(mfrow = c(1, 3))
matplot(NRout[1:1000, 2:3], type = "l", xlab = "Time", ylim = c(0, 110), ylab = "Population Density", main = "first equilibrium")

# Supply point 2
matplot(NRout2[1:1000, 2:3], type = "l", xlab = "Time", ylab = "Population Density", ylim = c(0, 110), main = "net effects")
matplot(NRoutDirect[1:1000, 2:3], type = "l", xlab = "Time", ylab = "Population Density", ylim = c(0, 110), main = "direct effects")

# Sensitivity of N1 to direct effects of change in S1
simDirSensitivity <- (dirEq[1] - preEq[1])/(newParms$S[1] - parms$S[1])
simDirSensitivity

# Sensitivity of N1 to indirect effects is net minus direct sensitivity
simIndSensitivity <- simNetSensitivity - simDirSensitivity
simIndSensitivity


# Section 3: Check whether the simulation matches algebra for defining direct and indirect effects

rho <- calc_substitutable_rho(parms)  # find rho based on model parameters
rho

calcDirSens <- calc_direct_sensitivity(parms = parms)
calcDirSens

# check that calculated direct sensitivity matches simulated
cbind(calculated = calcDirSens[1], simulated = simDirSensitivity)

# calc net sensitivity
calcNetSens <- calc_net_sensitivity(parms)

# check that calculated net sensitivity matches simulated
cbind(calculated = calcNetSens, simulated = simNetSensitivity)

# calculate indirect sensitivity using equation 21
calcIndSensSimple <- calc_substitutable_ifx_simple(parms)  # simplified formula 
cbind(calculated = calcIndSensSimple, simulated = simIndSensitivity) # compare to simulation

# calculate indirect sensitivity using equation 22 
calcIndSens <- calc_substitutable_ifx_mechanistic(dirN1 = calcDirSens[1], dirN2 = calcDirSens[2], parms = parms, rho = rho)
cbind( calculated = calcIndSens, simulated = simIndSensitivity)  # compare with simulation

# Run across a range of parameters and rho's and calculate indirect effects
parmsList <- list(NA)
parmsList[[1]] <- parms
rhoGradient <- NA
rhoGradient[1] <- calc_substitutable_rho(parms)
outputDF <- data.frame(simNetFX = NA, simDirFX = NA, simIndFX = NA, calcNetSens = NA, calcDirN1 = NA, calcDirN2 = NA, calcIndSens = NA, rho = NA)

outputDF[1, 1:3] <- simulate_substitutuble_effects(parms)
outputDF[1, 4:8] <- calc_substitutable_effects(parms)

for (i in 2:10) {
    newParms <- parms
    newParms$Q <- parms$Q + matrix(c(0, i, i, 0), nrow = 2, byrow = TRUE)
    parmsList[[i]] <- newParms
    outputDF[i, 1:3] <- simulate_substitutuble_effects(newParms)
    outputDF[i, 4:8] <- calc_substitutable_effects(newParms)
}

# Check that simulated and calculated values match columns.  
cbind(simulated = outputDF$simNetFX, calculated = outputDF$calcNetSens)
cbind(simulated = outputDF$simDirFX, calculated = outputDF$calcDirN1)
cbind(simulated = outputDF$simIndFX, calculated = outputDF$calcIndSens)

# Should be ~ equal (some error due to rounding etc.)
round(outputDF$simNetFX, 5) == round(outputDF$calcNetSens, 5)
round(outputDF$simDirFX, 5) == round(outputDF$calcDirN1, 5)
round(outputDF$simIndFX, 5) == round(outputDF$calcIndSens, 5)

