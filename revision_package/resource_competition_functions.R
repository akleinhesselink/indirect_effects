#### Online Appendix B: R-code for numerical simulation supporting "Indirect effects of environmental change in resource competition models"
#### Written on R version 3.0.2 (2013-09-25)
#### Andrew R. Kleinhesselink

require(deSolve)

########################################################
####### Section 1: Essential resource model analysis 
########################################################


RRmod = function(t, NR, parms) {
  
  #### Continuous two species resource model with essential resources
  
  with(as.list(c(NR, parms)), { 
    dN1dt = N1*min((r[1]*R1/(R1 + k[1,1]) - m[1]), (r[1]*R2/(R2 + k[1,2]) - m[1])) ## sp. 1 growth
    dN2dt = N2*min((r[2]*R1/(R1 + k[2,1]) - m[2]), (r[2]*R2/(R2 + k[2,2]) - m[2])) ## sp. 2 growth
    dR1dt = a[1]*(S[1]-R1) - ((dN1dt + m[1]*N1)*Q[1,1] + (dN2dt + m[2]*N2)*Q[2,1]) ## res. 1 change
    dR2dt = a[2]*(S[2]-R2) - ((dN1dt + m[1]*N1)*Q[1,2] + (dN2dt + m[2]*N2)*Q[2,2]) ## res. 2 change
    list(c(dN1dt, dN2dt, dR1dt, dR2dt))
  }) 
}

RRmod_constant_N2 = function(t, NR, parms, hold = 2) {
  
  #### Continuous two species resource model with essential resources
  #### This model finds equilibrium populatoin of focal species (N1) 
  #### while holding population of competitor (N2) constant. 
  
  with(as.list(c(NR, parms)), { 
    if (hold == 2){ 
      dN1dt = N1*min((r[1]*R1/(R1 + k[1,1]) - m[1]), (r[1]*R2/(R2 + k[1,2]) - m[1])) ## sp. 1 growth
      dN2dt = 0 # set to zero to hold N2 constant N2*min((r[1]*R1/(R1 + k[2,1]) - m[2]), (r[2]*R2/(R2 + k[2,2]) - m[2])) 
    }
    else if (hold ==1){ 
      dN1dt = 0 # N1*min((r[1]*R1/(R1 + k[1,1]) - m[1]), (r[1]*R2/(R2 + k[1,2]) - m[1])) ## sp. 1 growth
      dN2dt = N2*min((r[1]*R1/(R1 + k[2,1]) - m[2]), (r[2]*R2/(R2 + k[2,2]) - m[2])) ## sp. 2 growth
    }
    dR1dt = a[1]*(S[1]-R1) - ((dN1dt + m[1]*N1)*Q[1,1] + (dN2dt + m[2]*N2)*Q[2,1]) ## res. 1 change
    dR2dt = a[2]*(S[2]-R2) - ((dN1dt + m[1]*N1)*Q[1,2] + (dN2dt + m[2]*N2)*Q[2,2]) ## res. 2 change
    list(c(dN1dt, dN2dt, dR1dt, dR2dt))
  }) 
}

get_Rstars = function(parms){
  attach(parms)
  
  out = matrix(c(k[1,1]*m[1]/(r[1]-m[1]), k[1,2]*m[1]/(r[1]-m[1]),    #### R* for species one on R1 and R2
                 k[2,1]*m[2]/(r[2]-m[2]), k[2,2]*m[2]/(r[2]-m[2])),   #### R* for species two on R1 and R2
               nrow = 2, byrow = T) 
  detach(parms)
  return(out)
}

get_consumption = function(parms){ 
  attach(parms)
  out = c(C1 = Q[1,2]/Q[1,1], C2 = Q[2,2]/Q[2,1])      
  detach(parms)
  return(out)
}

getKs = function( parms ,  t = t) { 
  initR1 = parms$S[1]
  initR2 = parms$S[2]  
  N2eq = ode(y = c(N1 = 0, N2 = 1, R1 = initR1, R2 = initR2), times = t, func = RRmod, parms = parms) ### calcs N1, N2, R1, R2
  N1eq = ode(y = c(N1 = 1, N2 = 0, R1 = initR1, R2 = initR2), times = t, func = RRmod, parms = parms) ### calcs N1, N2, R1, R2
  K1 = N1eq[nrow(N1eq), "N1"]
  K2 = N2eq[nrow(N2eq), "N2"]
  return(list(K1, K2))
}

get_lvParms = function(parms ){
  Rstars = get_Rstars(parms)
  alpha = matrix(NA, nrow = 2, ncol = 2)  
  Q = parms$Q
  
  if (Rstars[1,1] < Rstars[1,2]) {
  #### case where species one is limited by resource two
    alpha = matrix(NA, nrow = 2, ncol = 2)  
    Q = parms$Q
    alphaT = Q[2,2]/Q[1,2]
    betaT = Q[1,1]/Q[2,1]
    K1 = getKs(parms, t = 1:1000)[[1]]
    K2 = getKs(parms, t = 1:1000)[[2]]
    alpha[1,1] = 1/K1
    alpha[2,2] = 1/K2
    alpha[1,2] = alphaT*alpha[1,1]
    alpha[2,1] = betaT*alpha[2,2]  
  }

  else if (Rstars[1,1] > Rstars[1,2]) {
    #### case where species one is limited by resource one
    alpha = matrix(NA, nrow = 2, ncol = 2)  
    Q = parms$Q
    alphaT = Q[1,2]/Q[2,2]
    betaT = Q[2,1]/Q[1,1]
    K1 = getKs(parms, t = 1:1000)[[1]]
    K2 = getKs(parms, t = 1:1000)[[2]]
    alpha[1,1] = 1/K1
    alpha[2,2] = 1/K2
    alpha[1,2] = alphaT*alpha[1,1]
    alpha[2,1] = betaT*alpha[2,2]  
  }
  
  
  return(list(alpha = alpha, K1 = K1, K2 = K2))
  
}
  
  
findRho = function(alpha){
  #### rho is the niche overlap measure from Chesson 2013
  rho = sqrt((alpha[1,2]*alpha[2,1])/(alpha[1,1]*alpha[2,2]))
  return(rho)  
}

plot_isocline = function(parms, RRout){ 
  #### Make Tilman style plot of species growth isoclines in relation to resource availability 
  #### for two essential resources and two competitors.  
  Rstars = get_Rstars(parms = parms)
  beta = get_consumption(parms = parms)
  with(as.list(c(parms, RRout, Rstars, beta)), {
    #plot(S[1]*2,S[2]*2,xlab="R1",ylab="R2",xlim=c(0,40),pch=16,ylim=c(0,40),type="n")

    plot(S[1]*2,S[2]*2,xlab="R1",ylab="R2",xlim=c(0,max(S[1])*2.2),pch=16,ylim=c(0,max(S[2])*2.2),type="n")
    
    points(RRout[1,4],RRout[1,5],xlab="R1",ylab="R2",xlim=c(0,max(S[1])*2.2),ylim=c(0,max(S[2])*2.2),type="p")
    lines(x=c(Rstars[1,1],Rstars[1,1],max(S[1])*2.2),y=c(max(S[2])*2.2,Rstars[1,2],Rstars[1,2]),lwd=2)
    lines(x=c(Rstars[2,1],Rstars[2,1],max(S[1])*2.2),y=c(max(S[2])*2.2,Rstars[2,2],Rstars[2,2]),lwd=2,col="red")
    
    CEQ=c(max(Rstars[1,1],Rstars[2,1]),max(Rstars[1,2], Rstars[2,2]))  ## Coexistence Equilibrium
    REQ=c(RRout[nrow(RRout), 4], RRout[nrow(RRout), 5])         ## Final Resource State (Resource Equilibrium)
    
    abline(CEQ[2] - C1*CEQ[1],C1,lty = "dashed")  ## species one consumption vector
    abline(CEQ[2] - C2*CEQ[1],C2,lty = "dashed", col="red") ## species two consumption vector
    
    lines(RRout[,4],RRout[,5])
    arrows(RRout[1, 4],RRout[1,5],RRout[3,4],RRout[3,5],length=0.05,angle=45,code=2)
    points(REQ[1], REQ[2], pch = 17, col = "dark green") ## label final resource state
    #text(REQ[1] + 2, REQ[2] - 2, cex = 0.5, "REQ", col = "dark green") ## label final resource state
    #lines(lR2, c(parms$S[2],parms$S[2]))
    #lines(c(parms$S[1],parms$S[1]), lR1)
  })
  legend("topright", 45, c("N1 ZNGI", "N2 ZNGI"), lty = 1, bty = "n", cex = 1.2, col = c(1,2))  
}  


plot_ifx_panel = function(RRout_net, RRout_direct, RRout_indirect, t){
  
  #### make three graphs representing how a change in resource supply affects the two species equilibrium
  #### left graph shows overall effect
  #### middle graph shows direct effect--i.e. effect of change when competitor is held constant
  #### right graph shows indirect effect--i.e. net effect minus direct effect  
  start = min(t)
  p1 = max(t)
  end = max(t)*2
  
  lines = c(1,1, 2)
  clrs = c("black", "red", "red")
  ylims = c(0, 1.1*(max(RRout_net[, 2], RRout_net[,3])))
  matplot( RRout_net[, 1], RRout_net[,2:3], type = "l", xlab ="time", ylab = "Population", lty = lines, col = clrs, ylim = ylims, 
           main = "Net effects") 
  abline(v=p1, lty = 2)
  
  RRout_direct = cbind(RRout_direct, RRout_direct[ , 3])
  RRout_direct[p1:end, 3] = NA #### don't graph second period of species 2 in direct effects graph
  RRout_direct[start:p1, 4] = NA
  
  matplot( RRout_direct[, 1], RRout_direct[,2:4], type = "l", xlab ="time", ylab = "", lty = lines, col = clrs, ylim = ylims, 
           main = "Direct effects") 
  abline(v=p1, lty = 2)
  matplot( RRout_indirect[, 1], RRout_indirect[,2:3], type = "l", xlab ="time", ylab = "", lty = lines, col = clrs, ylim = ylims, 
           main = "Indirect effects") 
  abline(v=p1, lty = 2)
  legend("bottomright", c("N1", "N2"), lty = 1, bty = "n", col = clrs, cex = 1.2, seg.len = 1)
}


changeSupply = function( S1delta = 0, S2delta = 0, NR, parms, t = 1:500){
  
  parms2 = parms
  parms2$S = parms2$S + c(S1delta, S2delta)
  
  out1 = ode(y = NR, times = t, func = RRmod, parms = parms) ### calcs N1, N2, R1, R2
  eq1 = out1[max(t), c(2:5)]
  out2 = ode(y = eq1, times = t, func = RRmod, parms = parms2)
  out2[, 1] = seq(max(t)+1, max(t)*2) 
  
  outNet = rbind(out1[, 1:3], out2[, 1:3])
  out3 = ode(y = eq1, times = t, func = RRmod_constant_N2, parms = parms2)
  out3[, 1] = out2[, 1]
  
  outDirect = rbind(out1[, 1:3], out3[, 1:3])
  out4 = out2[, 1:3]
  out4[, 2] = out2[, 2] - out3[, 2]
  
  out4[, 2] = eq1[1] + out4[,2]
  outIndirect = rbind(out1[, 1:3], out4[, 1:3])
  
  return(list(outNet = outNet, outDirect = outDirect, outIndirect = outIndirect))
}


getRho2 = function( parms, N1limitedbyR1 = TRUE){
  if (N1limitedbyR1 ) { 
    rho = with(parms, sqrt((Q[1,2]*Q[2,1])/(Q[1,1]*Q[2,2])))    
  }
  return(rho)
}

calcIFX1 = function( rho, dfxN1 = 1) { 
  ##### Draw indirect effect strength as a function of rho
  ifx1 = dfxN1*(rho^2/ ( 1 - rho^2))
  return(ifx1)  
}

calcFX2 = function( rho, focalSensitivity){
  #### Draw second type of indirect effects as function of rho and 
  #### focal species sensitivity 
  ifxtotal = - focalSensitivity*(rho^2/(1 - rho^2))
  return(ifxtotal)
}




########################################################

########################################################
####### Section 2: Substitutable resource model analysis 
########################################################

#### Substitutable resource model Functions 

subMod = function(t, NR, parms) {  
  #### Substitutable resource model  
  with(as.list(c(NR, parms)), { 
    
    dN1dt = N1*(r[1]*(w[1,1]*R1 + w[1,2]*R2 - tau[1])/(k[1] + w[1,1]*R1 + w[1,2]*R2 - tau[1])  - D) ## sp. 1 growth
    dN2dt = N2*(r[2]*(w[2,1]*R1 + w[2,2]*R2 - tau[2])/(k[2] + w[2,1]*R1 + w[2,2]*R2 - tau[2])  - D) ## sp. 2 growth
    dR1dt = D*(S[1] - R1) - q[1,1]*N1 - q[2,1]*N2
    dR2dt = D*(S[2] - R2) - q[1,2]*N1 - q[2,2]*N2  
    
    list(c(dN1dt, dN2dt, dR1dt, dR2dt))
  }) 
}


submod_constant_N2 = function(t, NR, parms, hold = 2) {
  
  #### Continuous two species resource model with substitutable resources
  #### This model finds equilibrium populatoin of focal species (N1) 
  #### while holding population of competitor (N2) constant. 
  
  with(as.list(c(NR, parms)), { 
    if (hold == 2){ 
      dN1dt = N1*(r[1]*(w[1,1]*R1 + w[1,2]*R2 - tau[1])/(k[1] + w[1,1]*R1 + w[1,2]*R2 - tau[1])  - D) ## sp. 1 growth
      dN2dt = 0 # set to zero to hold N2 constant N2*min((r[1]*R1/(R1 + k[2,1]) - m[2]), (r[2]*R2/(R2 + k[2,2]) - m[2])) 
    }
    else if (hold ==1){ 
      dN1dt = 0 # N1*min((r[1]*R1/(R1 + k[1,1]) - m[1]), (r[1]*R2/(R2 + k[1,2]) - m[1])) ## sp. 1 growth
      dN2dt = N2*(r[2]*(w[2,1]*R1 + w[2,2]*R2 - m[2])/(k[2] + w[2,1]*R1 + w[2,2]*R2 - m[2])  - D) ## sp. 1 growth
    }
    dR1dt = D*(S[1] - R1) - q[1,1]*N1 - q[2,1]*N2
    dR2dt = D*(S[2] - R2) - q[1,2]*N1 - q[2,2]*N2 
    list(c(dN1dt, dN2dt, dR1dt, dR2dt))
  }) 
}



calcRhoSub = function(parms ) { 
  rho = NA
  with(parms, { 
    rho = sqrt( ( (q[2,2] + q[2,1]*w[1,1]/w[1,2])*( q[1,2] + q[1,1]*w[2,1]/w[2,2]) )
                /( (q[1,2] + q[1,1]*w[1,1]/w[1,2])*( q[2,2] + w[2,1]*q[2,1]/w[2,2] ) ) )
    return(rho)
  } )
}


calcIntercepts = function(parms){ 
  #### Returns values for B in Tilman 1981 Appendix 
  with(parms, { 
    B1 = ((D/r[1])*(k[1] - tau[1]) + tau[1])/(w[1,2]*(1 - D/r[1]))
    B2 = ((D/r[2])*(k[2] - tau[2]) + tau[2])/(w[2,2]*(1 - D/r[2]))
    return(list(B1, B2))
  })  
}

plotZingis = function(parms, B1, B2, ymax = 1000 ){
  #### Plot Zero-Net Growth Isoclines for species one and two
  with(parms, { 
    xend = max((B1/(w[1,1]/w[1,2])),(B2/(w[2,1]/w[2,2])))
    yl = max(B1, B2, ymax)
    curve(B1 - (w[1,1]/w[1,2])*x, 0, xend, col = 1, ylim = c(0, yl), xlab = 'Resource 1', 
          ylab = 'Resource 2')
    curve(B2 - (w[2,1]/w[2,2])*x, 0, xend, add= TRUE, col =2)     
  })  
}

calcRstars = function(parms, B1, B2) { 
  #### calculate the point at which two species Zero Net Growth Isoclines Cross
  #### gives the equilibrium resource availability in competition
  with(parms, { 
    R1star = (B2 - B1)/(w[2,1]/w[2,2] - w[1,1]/w[1,2])
    R2star = (B2 - (w[2,1]/w[2,2])*((B2 - B1)/(w[2,1]/w[2,2] - w[1,1]/w[1,2])))
    
    return(list(R1star, R2star))
    
  })
}


plotConsumption = function(parms, B1, B2, R1star, R2star, species2Color = 2 ){ 
  ##### plot consumption vector on graph with resource one on x-axis 
  ##### and resource two on the y-axis.  Intercept determined by 
  ##### resource equilibrium. 
  with(parms, { 
    cv1 = q[1,2]/q[1,1]
    cv2 = q[2,2]/q[2,1]
    
    points(R1star, R2star)      
    
    curve( R2star - (R1star*cv1) + cv1*x , add= TRUE , lty = 2)
    curve( R2star - (R1star*cv2) + cv2*x, add = TRUE, lty = 2, col = species2Color)
    
  })
  
}


calcLVparms = function(parms, B){ 
  ##### Calculate Lotka-Volterra equivalent competition coefficients 
  ##### based on substitutable resource model parameters 
  with(parms, { 
    
    K = c(  (D*(S[2] + S[1]*w[1,1]/w[1,2] - B[1]))/(q[1,2] + w[1,1]*q[1,1]/w[1,2]), 
            (D*(S[2] + S[1]*w[2,1]/w[2,2] - B[2]))/(q[2,2] + w[2,1]*q[2,1]/w[2,2])  )
    
    alphaT = (q[2,2] + q[2,1]*w[1,1]/w[1,2])/(q[1,2] + q[1,1]*w[1,1]/w[1,2])
    betaT =  (q[1,2] + q[1,1]*w[2,1]/w[2,2])/(q[2,2] + q[2,1]*w[2,1]/w[2,2])
    
    alpha = matrix(c(NA), nrow = 2, ncol = 2, byrow = 2)
    
    alpha[1,1] = 1/K[1]
    alpha[2,2] = 1/K[2]
    alpha[1,2] = alpha[1,1]*alphaT
    alpha[2,1] = alpha[2,2]*betaT
    
    LVP = list(K = K, alphaT = alphaT, betaT = betaT, alpha = alpha)
    return(LVP)
  })
}

calcDirSensitivity = function(parms){ 
  #### Calculate direct sensitivity ot a change in supply of resource one
  with(parms, {  
    dirSensN1 = D*(w[1,1]/w[1,2])/(q[1,2] + (w[1,1]*q[1,1])/w[1,2])
    dirSensN2 = D*(w[2,1]/w[2,2])/(q[2,2] + (w[2,1]*q[2,1])/w[2,2])
    return(c(dirSensN1, dirSensN2))
  })
}

calcNetSensitivity = function(dirN1, dirN2, lvp){
  #### dirN1 and dirN2 give the direct sensitivity 
  #### of N1 and N2 to change in S1 
  #### calculated in 'calcDirSensitivity' function 
  
  with(lvp, { 
    Net = ((1 - alphaT*betaT)^(-1))*(dirN1 - alphaT*dirN2)
    return(Net)
  }) 
}

calcNetSensitivity2 = function( parms) { 

  with(parms, { 
    Net = D/q[1,1]*( 1 / (1 - q[1,2]*q[2,1]/(q[1,1]*q[2,2]) )) 
    return(Net)  
  } ) 
}

calcIFX = function( dirN1, dirN2, Beta, rho){ 
  #### Uses LV parameters to calculate indirect sensitivity to 
  #### change in resource supply 
  indSensN1 = (dirN1 - (1/Beta)*dirN2)*((rho^2)/(1 - rho^2))
  return( calcIndSens = indSensN1 )  
}

calcIFX2 = function(dirN1, dirN2, parms, rho) { 
  #### like calcIFX but does not use LV parameters
  #### goes directly from mechanistic parameters 
  #### to indirect effects strength 
  with(parms, { 
    C = ( q[2,2] + q[2,1]*w[2,1]/w[2,2]) / ( q[1,2] + q[1,1]*w[2,1]/w[2,2] )
    IFX2 = (dirN1 - C*dirN2)*(rho^2/(1-rho^2))
    return(IFX2)
  })
}

simulatedEffects = function( parms ){ 
  
  NR = c(N = c(1, 1), R = c(4000,4000))   ### initial populations and resource levels
  NRout = ode(y = NR, func= subMod, t = 1:2000, parms = parms, method = 'rk4' )
  preEq = as.matrix(tail(NRout, 1))[,-1]
  
  newParms = parms
  newParms$S = parms$S*c(1.01, 1)
  NR2 = ode(y = preEq, func = subMod, t = 1:2000, parms = newParms, method = 'rk4')
  NR3 = ode(y = preEq, func = submod_constant_N2, t = 1:2000, parms = newParms, method = 'rk4')
  Eq2 = as.matrix(tail(NR2, 1))[, -1]
  Eq3 = as.matrix(tail(NR3, 1))[, -1]
  
  netFX = (Eq2 - preEq)/(newParms$S[1] - parms$S[1])
  dirFX = (Eq3 - preEq)/(newParms$S[1] - parms$S[1])
  indFX = netFX[1] - dirFX[1]
  
  return(c(netFX[1], dirFX[1], indFX[1]))
  
}

analyticalCalculations = function( parms ){ 
  B = calcIntercepts(parms=parms)
  LVP = calcLVparms(parms=parms, B=c(B[[1]], B[[2]]))
  rho = calcRhoSub( parms= parms)
  dirSens = calcDirSensitivity(parms=parms)
  netSens = calcNetSensitivity(dirN1= dirSens[1], dirN2 = dirSens[2], lvp= LVP)
  indSens = calcIFX(dirN1=dirSens[1], dirN2 = dirSens[2], Beta = LVP$betaT, rho=rho)  
  indSens2 = calcIFX2(dirN1 = dirSens[1], dirN2= dirSens[2], parms = parms, rho = rho)
  return(c(netSens, dirSens[1], dirSens[2], indSens, indSens2, rho, LVP$betaT))
}

############################################################
####### End Section 2: Substitutable resource model analysis 
############################################################

