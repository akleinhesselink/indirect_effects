# essential resource model functions 
# R-functions for numerical simulation supporting 'Indirect
# effects of environmental change in resource competition models' 
# Written on R version 3.2.1 (2015-06-18) 
# Author: Andrew R. Kleinhesselink

essential_mod <- function(t, NR, parms) {
  
  # Continuous two species resource model with essential resources. 
  # Refer to equations 6 and 7 in the main text. 
  # Intended use is to be passed as a parameter to the ode() function 
  # in the deSolve package for ordinary differential equations.  
  
  with(as.list(c(NR, parms)), {
    dN1dt <- N1 * min((r[1] * R1/(R1 + k[1, 1]) - m[1]), 
                      (r[1] * R2/(R2 + k[1, 2]) - m[1]))  # sp. 1 growth
    dN2dt <- N2 * min((r[2] * R1/(R1 + k[2, 1]) - m[2]), 
                      (r[2] * R2/(R2 + k[2, 2]) - m[2]))  # sp. 2 growth
    dR1dt <- a[1] * (S[1] - R1) - ((dN1dt + m[1] * N1) * Q[1, 1] + 
                                     (dN2dt + m[2] * N2) * Q[2, 1])  # res. 1 change
    dR2dt <- a[2] * (S[2] - R2) - ((dN1dt + m[1] * N1) * Q[1, 2] + 
                                     (dN2dt + m[2] * N2) * Q[2, 2])  # res. 2 change
    list(c(dN1dt, dN2dt, dR1dt, dR2dt))
  })
}

essential_mod_hold_comp <- function(t, NR, parms, hold = 2) {
  
  # Continuous two species resource model with essential resources This model finds
  # equilibrium populatoin of focal species (N1) while holding population of
  # competitor (N2) constant.
  # Intended use is to be passed as a parameter to the ode() function 
  # in the deSolve package for ordinary differential equations.  
  
  with(as.list(c(NR, parms)), {
    if (hold == 2) {
      dN1dt <- N1 * min((r[1] * R1/(R1 + k[1, 1]) - m[1]), 
                        (r[1] * R2/(R2 + k[1, 2]) - m[1]))  # sp. 1 growth
      dN2dt <- 0  # set to zero to hold N2 constant 
    } else if (hold == 1) {
      dN1dt <- 0  # hold N1 constant
      dN2dt <- N2 * min((r[1] * R1/(R1 + k[2, 1]) - m[2]), 
                        (r[2] * R2/(R2 +  k[2, 2]) - m[2]))  # sp. 2 growth
    }
    dR1dt <- a[1] * (S[1] - R1) - ((dN1dt + m[1] * N1) * Q[1, 1] + 
                                     (dN2dt + m[2] * N2) * Q[2, 1])  # res. 1 change
    dR2dt <- a[2] * (S[2] - R2) - ((dN1dt + m[1] * N1) * Q[1, 2] + 
                                     (dN2dt + m[2] *  N2) * Q[2, 2])  # res. 2 change
    list(c(dN1dt, dN2dt, dR1dt, dR2dt))
  })
}

calc_rstars <- function(parms) {

  # calculate R stars for each species for an essential 
  # resource model 
  
  with(parms, { out <- matrix(
                c(k[1, 1] * m[1]/(r[1] - m[1]), 
                  k[1, 2] * m[1]/(r[1] - m[1]), 
                  k[2, 1] * m[2]/(r[2] - m[2]), 
                  k[2, 2] * m[2]/(r[2] - m[2])), nrow = 2, byrow = T)  
                return(out) 
                } 
       ) 
}

calc_essential_consumption <- function(parms) {
  
  # calculate consumption rates for each species in an 
  # essential resource model 
  
  with(parms, { 
    out <- c(C1 = Q[1, 2]/Q[1, 1], C2 = Q[2, 2]/Q[2, 1])
    return(out) 
    }
  ) 
}


plot_essential_zngis <- function(parms, RRout) {
  
  # Make Tilman style plot of species zero-net growth isoclines 
  # in relation to resource availability for two essential resources 
  # and two competing consumers. 
  # Similar (but not identical) to figure 2 in the main text 
  
  Rstars <- calc_rstars(parms = parms)
  beta <- calc_essential_consumption(parms = parms)
  with(as.list(c(parms, RRout, Rstars, beta)), {
    
    plot(S[1] * 2, S[2] * 2, xlab = "R1", ylab = "R2", 
         xlim = c(0, max(S[1]) * 2.2), pch = 16, 
         ylim = c(0, max(S[2]) * 2.2), type = "n")
    
    points(RRout[1, 4], RRout[1, 5], xlab = "R1", ylab = "R2", 
           xlim = c(0, max(S[1]) * 2.2), 
           ylim = c(0, max(S[2]) * 2.2), type = "p")
    
    lines(x = c(Rstars[1, 1], Rstars[1, 1], max(S[1]) * 2.2), 
          y = c(max(S[2]) * 2.2, Rstars[1, 2], Rstars[1, 2]), lwd = 2)
    lines(x = c(Rstars[2, 1], Rstars[2, 1], max(S[1]) * 2.2), 
          y = c(max(S[2]) * 2.2, Rstars[2, 2], Rstars[2, 2]), lwd = 2, col = "red")
    
    CEQ <- c(max(Rstars[1, 1], Rstars[2, 1]), max(Rstars[1, 2], Rstars[2, 2]))  # Coexistence Equilibrium
    REQ <- c(RRout[nrow(RRout), 4], RRout[nrow(RRout), 5])  # Final Resource State (Resource Equilibrium)
    
    abline(CEQ[2] - C1 * CEQ[1], C1, lty = "dashed")  # species one consumption vector
    abline(CEQ[2] - C2 * CEQ[1], C2, lty = "dashed", col = "red")  # species two consumption vector
    
    lines(RRout[, 4], RRout[, 5])
    arrows(RRout[1, 4], RRout[1, 5], RRout[3, 4], RRout[3, 5], length = 0.05, 
           angle = 45, code = 2)
    points(REQ[1], REQ[2], pch = 17, col = "dark green")  # label final resource state

  })
  legend("topright", 45, c("N1 ZNGI", "N2 ZNGI"), lty = 1, bty = "n", cex = 1.2, 
         col = c(1, 2))
}


plot_ifx_panel <- function(RRout_net, RRout_direct, RRout_indirect, t) {
  
  # make three graphs representing how a change in resource supply affects the two
  # species equilibrium. Left graph shows overall effect; middle graph shows direct
  # effect--i.e. effect of change when competitor is held constant; and right graph
  # shows indirect effect--i.e. net effect minus direct effect
  
  start <- min(t)
  p1 <- max(t)
  end <- max(t) * 2
  
  lines <- c(1, 1, 2)
  clrs <- c("black", "red", "red")
  ylims <- c(0, 1.1 * (max(RRout_net[, 2], RRout_net[, 3])))
  matplot(RRout_net[, 1], RRout_net[, 2:3], type = "l", xlab = "time", ylab = "Population", 
          lty = lines, col = clrs, ylim = ylims, main = "Net effects")
  abline(v = p1, lty = 2)
  
  RRout_direct <- cbind(RRout_direct, RRout_direct[, 3])
  RRout_direct[p1:end, 3] <- NA  # don't graph second period of species 2 in direct effects graph
  RRout_direct[start:p1, 4] <- NA
  
  matplot(RRout_direct[, 1], RRout_direct[, 2:4], type = "l", xlab = "time", ylab = "", 
          lty = lines, col = clrs, ylim = ylims, main = "Direct effects")
  abline(v = p1, lty = 2)
  matplot(RRout_indirect[, 1], RRout_indirect[, 2:3], type = "l", xlab = "time", 
          ylab = "", lty = lines, col = clrs, ylim = ylims, main = "Indirect effects")
  abline(v = p1, lty = 2)
  legend("bottomright", c("N1", "N2"), lty = 1, bty = "n", col = clrs, cex = 1.2, 
         seg.len = 1)
}


simulate_essential_effects <- function(S1delta = 0, S2delta = 0, NR, parms, t = 1:500) {
  
  # utility function to change resource supply and return simulated net, 
  # direct and indirect effects. 
  
  parms2 <- parms
  parms2$S <- parms2$S + c(S1delta, S2delta)
  
  out1 <- ode(y = NR, times = t, func = essential_mod, parms = parms)  # calcs N1, N2, R1, R2
  eq1 <- out1[max(t), c(2:5)]
  out2 <- ode(y = eq1, times = t, func = essential_mod, parms = parms2)
  out2[, 1] <- seq(max(t) + 1, max(t) * 2)
  
  outNet <- rbind(out1[, 1:3], out2[, 1:3])
  out3 <- ode(y = eq1, times = t, func = essential_mod_hold_comp, parms = parms2)
  out3[, 1] <- out2[, 1]
  
  outDirect <- rbind(out1[, 1:3], out3[, 1:3])
  out4 <- out2[, 1:3]
  out4[, 2] <- out2[, 2] - out3[, 2]
  
  out4[, 2] <- eq1[1] + out4[, 2]
  outIndirect <- rbind(out1[, 1:3], out4[, 1:3])
  
  return(list(outNet = outNet, outDirect = outDirect, outIndirect = outIndirect))
}


calc_essential_rho <- function(parms, N1limitedbyR1 = TRUE) {
  
  # calculate rho directly from mechanistic resource uptake parameters
  # see equation 9 in the main text 
  
  if (N1limitedbyR1) {
    rho <- with(parms, sqrt((Q[1, 2] * Q[2, 1])/(Q[1, 1] * Q[2, 2])))
  }
  return(rho)
}

calc_essential_ifx1 <- function(rho, dfxN1 = 1) {
  
  # Indirect effect strength as a function of rho
  # see equation 13 in the main text 
  
  ifx1 <- dfxN1 * (rho^2/(1 - rho^2))
  return(ifx1)
}

calc_essential_ifx2 <- function(rho, focalSensitivity) {
  
  # Indirect effects when the non-limiting resource changes. Determined 
  # by rho and focal species sensitivity
  # see equation 15 in the main text 

  ifx <- -focalSensitivity * (rho^2/(1 - rho^2))
  return(ifx)
}
