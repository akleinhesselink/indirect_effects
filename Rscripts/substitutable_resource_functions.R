# substitutable resource model functions 
# R-functions for numerical simulation supporting 'Indirect
# effects of environmental change in resource competition models' 
# Written on R version 3.2.1 (2015-06-18) 
# Author: Andrew R. Kleinhesselink

substitutable_mod <- function(t, NR, parms) {

  # Substitutable resource model 
  # Refer to equations 16 and 17 in the main text.
  # Intended use is to be passed as a parameter to the ode() function 
  # in the deSolve package for ordinary differential equations.  
  
  with(as.list(c(NR, parms)), {
    
    dN1dt <- N1 * (r[1] * (W[1, 1] * R1 + W[1, 2] * R2 - tau[1])/
                     (k[1] + W[1,1] * R1 + W[1, 2] * R2 - tau[1]) - D)  # sp. 1 growth
    dN2dt <- N2 * (r[2] * (W[2, 1] * R1 + W[2, 2] * R2 - tau[2])/
                     (k[2] + W[2,1] * R1 + W[2, 2] * R2 - tau[2]) - D)  # sp. 2 growth
    dR1dt <- D * (S[1] - R1) - Q[1, 1] * N1 - Q[2, 1] * N2
    dR2dt <- D * (S[2] - R2) - Q[1, 2] * N1 - Q[2, 2] * N2
    
    list(c(dN1dt, dN2dt, dR1dt, dR2dt))
  })
}


substitutable_mod_hold_comp <- function(t, NR, parms, hold = 2) {
  
  # Continuous two species resource model with substitutable resources This model
  # finds equilibrium populatoin of focal species (N1) while holding population of
  # competitor (N2) constant.
  # Intended use is to be passed as a parameter to the ode() function 
  # in the deSolve package for ordinary differential equations.  
  
  with(as.list(c(NR, parms)), {
    if (hold == 2) {
      dN1dt <- N1 * (r[1] * (W[1, 1] * R1 + W[1, 2] * R2 - tau[1])/
                       (k[1] + W[1, 1] * R1 + W[1, 2] * R2 - tau[1]) - D)  # sp. 1 growth
      dN2dt <- 0  # set to zero to hold N2 constant 
    } else if (hold == 1) {
      dN1dt <- 0  # set N1 constant 
      dN2dt <- N2 * (r[2] * (W[2, 1] * R1 + W[2, 2] * R2 - m[2])/
                       (k[2] + W[2, 1] * R1 + W[2, 2] * R2 - m[2]) - D)  # sp. 1 growth
    }
    dR1dt <- D * (S[1] - R1) - Q[1, 1] * N1 - Q[2, 1] * N2
    dR2dt <- D * (S[2] - R2) - Q[1, 2] * N1 - Q[2, 2] * N2
    list(c(dN1dt, dN2dt, dR1dt, dR2dt))
  })
}

calc_substitutable_rho <- function(parms) {
  
  # calculate rho directly from mechanistic resource model parameters
  # refer to equation 18 in the main text 
  
  rho <- NA
  with(parms, {
    rho <- sqrt(((Q[2, 2] * W[1, 2] + Q[2, 1] * W[1, 1]) * (Q[1, 2] * W[2, 2] + Q[1, 1] * W[2, 1]))/
                  ((Q[1, 2] * W[1, 2] + Q[1, 1] * W[1, 1]) * (Q[2,2] * W[2, 2] + W[2, 1] * Q[2, 1])))
    return(rho)
  })
}


calc_intercepts <- function(parms) {
  
  # Returns values for ZNGI y-intercepts
  
  with(parms, {
    B1 <- ((D/r[1]) * (k[1] - tau[1]) + tau[1])/(W[1, 2] * (1 - D/r[1]))
    B2 <- ((D/r[2]) * (k[2] - tau[2]) + tau[2])/(W[2, 2] * (1 - D/r[2]))
    return(list(B1, B2))
  })
}

plot_substitutable_zngis <- function(parms, B1, B2, ymax = 1000) {
  
  # Plot Tilman-style zero net growth isoclines for species one and two
  # compare to figure A1 in Appendix A
  
  with(parms, {
    xend <- max((B1/(W[1, 1]/W[1, 2])), (B2/(W[2, 1]/W[2, 2])))
    yl <- max(B1, B2, ymax)
    curve(B1 - (W[1, 1]/W[1, 2]) * x, 0, xend, col = 1, ylim = c(0, yl), xlab = "Resource 1", 
          ylab = "Resource 2")
    curve(B2 - (W[2, 1]/W[2, 2]) * x, 0, xend, add = TRUE, col = 2)
  })
}

calc_substitutable_rstars <- function(parms, B1, B2) {
  
  # Calculate the point at which zero net growth isoclines for two 
  # competing species cross.  This gives the equilibrium resource 
  # availability in competition
  
  with(parms, {
    R1star <- (B2 - B1)/(W[2, 1]/W[2, 2] - W[1, 1]/W[1, 2])
    R2star <- (B2 - (W[2, 1]/W[2, 2]) * ((B2 - B1)/
                                           (W[2, 1]/W[2, 2] - W[1, 1]/W[1,2])))
    
    return(list(R1star, R2star))
    
  })
}


plot_substitutable_consumption <- function(parms, B1, B2, R1star, R2star, species2Color = 2) {
  
  # plot consumption vector on graph with resource one on x-axis and resource two
  # on the y-axis.  Intercept determined by resource equilibrium.
  # Compare to figure A1 in Appendix A
  
  with(parms, {
    cv1 <- Q[1, 2]/Q[1, 1]
    cv2 <- Q[2, 2]/Q[2, 1]
    
    # Draw on graph
    points(R1star, R2star)    
    curve(R2star - (R1star * cv1) + cv1 * x, add = TRUE, lty = 2)
    curve(R2star - (R1star * cv2) + cv2 * x, add = TRUE, lty = 2, col = species2Color)
    
  })
  
}

calc_direct_sensitivity <- function(parms) {
  
  # Calculate direct sensitivity ot a change in resource 
  # supply point for resource one.
  # Refer to equation 20 in the main text 
  
  with(parms, {
    dirSensN1 <- D * (W[1, 1])/(Q[1, 2] * W[1, 2] + W[1, 1] * Q[1, 1])
    dirSensN2 <- D * (W[2, 1])/(Q[2, 2] * W[2, 2] + W[2, 1] * Q[2, 1])
    return(c(dirSensN1, dirSensN2))
  })
}

calc_net_sensitivity <- function(parms) {
  
  # Calculate net sensitivity to a change in resource 
  # supply point for substitutable resource model.
  # Refer to equation 19 in the main text 
  
  with(parms, {
    Net <- D/Q[1, 1] * (1/(1 - Q[1, 2] * Q[2, 1]/(Q[1, 1] * Q[2, 2])))
    return(Net)
  })
}

calc_substitutable_ifx <- function(dirN1, dirN2, Beta, rho) {
  
  # Calculate indirect sensitivity to change in resource supply
  # point when provided Beta and rho
  
  indSensN1 <- (dirN1 - (1/Beta) * dirN2) * ((rho^2)/(1 - rho^2))
  return(calcIndSens = indSensN1)
}

calc_substitutable_ifx_mechanistic <- function(dirN1, dirN2, parms, rho) {
  
  # Calcluate indirect sensitivity in a substitutable resource model. 
  # Calculates Beta directly from mechanistic parameters. 
  # Refer to equation 22 in the main text 
  
  with(parms, {
    Beta <- (Q[2, 2] + Q[2, 1] * W[2, 1]/W[2, 2])/(Q[1, 2] + Q[1, 1] * W[2, 1]/W[2,2])
    IFX2 <- (dirN1 - Beta * dirN2) * (rho^2/(1 - rho^2))
    return(IFX2)
  })
}

calc_substitutable_ifx_simple <- function(parms) {
  
  # Calculates indirect sensitivity with simplified formula directly from model
  # parameters, not including rho. 
  # Refer to equation 21 in the main text 

  NetSenseN1 <- calc_net_sensitivity(parms)
  
  DirSenseN1 <- with(parms, { D * (W[1, 1]/W[1, 2])/(Q[1, 2] + (W[1, 1] * Q[1, 1])/W[1, 2])
  })
  IFX <- NetSenseN1 - DirSenseN1
  return(IFX)
}

simulate_substitutuble_effects <- function(parms) {
  
  # helper function to find simulated indirect effects 
  # in a substitutable resource model 
  
  NR <- c(N = c(1, 1), R = c(4000, 4000))  # initial populations and resource levels
  NRout <- ode(y = NR, func = substitutable_mod, t = 1:2000, parms = parms, method = "rk4")
  preEq <- as.matrix(tail(NRout, 1))[, -1]
  
  newParms <- parms
  newParms$S <- parms$S * c(1.01, 1)
  NR2 <- ode(y = preEq, func = substitutable_mod, t = 1:2000, parms = newParms, 
             method = "rk4")
  NR3 <- ode(y = preEq, func = substitutable_mod_hold_comp, t = 1:2000, parms = newParms, 
             method = "rk4")
  Eq2 <- as.matrix(tail(NR2, 1))[, -1]
  Eq3 <- as.matrix(tail(NR3, 1))[, -1]
  
  netFX <- (Eq2 - preEq)/(newParms$S[1] - parms$S[1])
  dirFX <- (Eq3 - preEq)/(newParms$S[1] - parms$S[1])
  indFX <- netFX[1] - dirFX[1]
  
  return(c(netFX[1], dirFX[1], indFX[1]))
  
}

calc_substitutable_effects <- function(parms) {
  
  # helper function to calculate net, direct and indirect 
  # effects for a substitutable resource model 
  
  B <- calc_intercepts(parms = parms)
  rho <- calc_substitutable_rho(parms = parms)
  dirSens <- calc_direct_sensitivity(parms = parms)
  netSens <- calc_net_sensitivity(parms = parms)
  indSens <- calc_substitutable_ifx_mechanistic(dirN1 = dirSens[1], dirN2 = dirSens[2], parms = parms, rho = rho)
  return(c(netSens, dirSens[1], dirSens[2], indSens, rho))
}
