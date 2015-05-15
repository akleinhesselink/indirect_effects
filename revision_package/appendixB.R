#### Online Appendix B: R-code for numerical simulation supporting "Indirect effects of environmental change in resource competition models"
#### Written on R version 3.0.2 (2013-09-25)
#### Andrew R. Kleinhesselink

require(deSolve)

########################################################
####### Section 1: Essential resource model analysis 
########################################################


RRmod = function(t, NR, parms) {
  
  #### RRmod "Resource Ratio Model" 
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
  
  #### RRmod "Resource Ratio Model" 
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

source('resource_competition_functions.R')

###### Demonstrates Two species two essential resource competition model of coexistence
NR = c(N = c(.01, .01), R = c(20,20))   ### initial populations and resource levels

#### set parameters
parms = list(r  = c(1, 1),                                         ### growth rates
             k = matrix( c(35, 10, 10, 35), nrow = 2, byrow = T),  ### resource half saturation for species one and two
             m  = c(0.21, 0.2),                                     ### mortality rates 
             a  = c(0.2, 0.2),                                     ### resource growth rate
             Q = matrix(c(3   , 1,
                          2   , 4   ), nrow = 2, byrow = T),       ### sp. one and two per capita use of resource one and two
             S = as.numeric(NR[3:4]))                              ### resource supply rates 

Rstars = get_Rstars(parms = parms)
Rstars

Rstars[1,1] > Rstars[1,2]

beta = get_consumption(parms = parms) #### beta is a consumption vector 

getKs(parms=parms, t = 1:1000)

lvParms = get_lvParms(parms = parms) #### calculate Lotka-Volterra params
lvParms

getRho2( parms = parms) #### calculate rho directly from parameters

rho = findRho(lvParms$alpha) #### calculate rho from Lotka-Volterra parms
rho

##### run model 
t = 1:1000
RRout = ode(y = NR, times = t, func = RRmod, parms = parms) ### calcs N1, N2, R1, R2

##### plot simulation to equilibrium
par(mfrow=c(2,2),pty="s",tcl=-0.2,mgp=c(2.2,0.5,0),mar=c(4,5,1,1), las = "1", cex.lab = 1.2)
clrs = c("green", "blue")

matplot( RRout[, 1], RRout[,2:3], type = "l", xlab ="time", ylab = "Population", 
         lty = 1:2, col = clrs) 
legend("topright", c("Sp. 1", "Sp. 2"), lty = 1:2, bty = "n", col = clrs, cex = 1.2)

matplot( RRout[, 1], RRout[,4:5], type = "l", xlab ="time", ylab = "Resources") 
legend("topright", c("R1", "R2"), lty = 1:3, bty = "n", cex = 1.2)

CEQ=c(max(Rstars[1,1],Rstars[2,1]),max(Rstars[1,2], Rstars[2,2]))  #### Coexistence Equilibrium
REQ=c(RRout[nrow(RRout), 4], RRout[nrow(RRout), 5])         #### Final Resource State (Resource Equilibrium)
lR2 = c(CEQ[2]-(1/beta[1])*CEQ[1] + (1/beta[1])*parms$S[2], CEQ[2]-(1/beta[2])*CEQ[2] + 
          (1/beta[2])*parms$S[2]) 
lR1 = c(CEQ[2]-beta[1]*CEQ[1] + beta[1]*parms$S[1], CEQ[2]-beta[2]*CEQ[2] + 
          beta[2]*parms$S[1])

par(mfrow=c(1,1),pty="s",tcl=-0.2,mgp=c(2.2,0.5,0),mar=c(4,5,1,1), las = "1", cex.lab = 1.2)
plot_isocline(parms, RRout) #### Resource isocline plot
text(35, 25, bquote( rho ~ '=' ~ .(rho) ) )



#### 
#### run the model with a change in resource supply and observe effects on N1* 
#### 

newT = 5000
S1delta = 1
S2delta = 0

out = changeSupply(S1delta, S2delta, NR, parms, t = 1:newT)
par(mfrow=c(1,3),pty="s",tcl=-0.2,mgp=c(2.2,0.5,0),mar=c(3,5,2,1), 
    oma = c(3,1,3,1), las = "1", cex.lab = 1.2)
plot_ifx_panel(out[[1]], out[[2]], out[[3]], t = newT)   ##### plot perturbation in S1


#### resource 2 sensitivity
S1delta2 = 0
S2delta2 = 1

out2 = changeSupply(S1delta2, S2delta2, NR, parms, t = 1:newT)
plot_ifx_panel(out2[[1]], out2[[2]], out2[[3]], t = newT) ##### plot perturbation in S2

#### Check whether the simulation matches algebra 
#### for defining direct and indirect effects 
t = 1:1000
NR1 = ode(y = NR, times = t, func = RRmod, parms = parms) ### calcs N1, N2, R1, R2
preEq = as.numeric(tail(NR1[, 2:5], 1)) #### first equilibrium 
names(preEq) <- names(NR)

newParms = parms
newParms$S = parms$S*c(1.01, 1) #### change S1 by 1%
deltaS1 = newParms$S[1] - parms$S[1]

NR2 = ode(y = preEq, times = t, func=RRmod, parms = newParms)
Eq2 = as.numeric(tail(NR2[ , 2:5], 1)) #### net effects equilibrium 
names(Eq2) <- names(NR)

NR2direct = ode(y = preEq, times = t, func = RRmod_constant_N2, parms = newParms)
Eq2direct = as.numeric(tail(NR2direct[ , 2:5], 1)) #### just direct effects equilibrium 
names(Eq2direct) <- names(NR)

indFX = Eq2[1] - Eq2direct[1] #### Indirect effects 
indFX 

dirFX = Eq2direct[1] - preEq[1] #### direct effects
dirFX
simDirSens = dirFX/(deltaS1) #### direct effects sensitivity

simIndSens = indFX/(deltaS1) #### Indirect effects sensitivity from simulation
simIndSens

#### This section calculates direct and indirect effect analytically 
rho = getRho2(parms)

dirFXN1_calc = with(parms, (a[1])/(m[1]*Q[1,1])) 
indFXN1_calc = with(parms, ((a[1]/(m[1]*Q[1,1]))*(rho^2/(1-rho^2)))) 

cbind( dirFXN1_calc, simDirSens)  #### check that calculated and simulated values match
cbind( indFXN1_calc, simIndSens ) ##### check that calculated and simulated values match


####### Section for calculating effects of a change in the supply of the non-limiting resource 
t = 1:1000
NR1 = ode(y = NR, times = t, func = RRmod, parms = parms) ### calcs N1, N2, R1, R2
preEq = as.numeric(tail(NR1[, 2:5], 1)) #### first equilibrium 
names(preEq) <- names(NR)

newParms = parms
newParms$S = parms$S*c(1, 1.01) #### change S2 by 1%
deltaS2 = newParms$S[2] - parms$S[2]

NR2 = ode(y = preEq, times = t, func=RRmod, parms = newParms)
Eq2 = as.numeric(tail(NR2[ , 2:5], 1)) #### net effects equilibrium 
names(Eq2) <- names(NR)

indFX2 = Eq2[1] - preEq[1] #### indirect effects are equal to net effects 
indFX2
simIndSens2 = indFX2/(deltaS2) #### Indirect effects sensitivity 
simIndSens2

#### This section calculates direct and indirect effect analytically 
rho = getRho2(parms)

indFX2_calc = with(parms, -((a[1]/(m[1]*Q[1,2]))*(rho^2/(1-rho^2)))) 

##### check that calculated and simulated values match
cbind( indFX2_calc, simIndSens2 ) 


#############################################################
#### Generate figures for main text of manuscript: 
#############################################################

######## Figure settings 
par(mfrow= c(1,1), oma = c(1,1,1,1), mar = c(5,5.5,1,1), las = 1, cex.axis = 0.8)

######## Figure labels and titles 
yl1 = expression("Indirect effects" ~ ~ ~ italic(frac('dN'['F indirect']^'*','dS'[1])))
yl2 = expression("Indirect effects" ~ ~ ~ italic(frac('dN'['F indirect']^'*', 'dS'[2])))
xl = expression( "Niche overlap" ~ "("*italic(rho)*")")
xldirect = expression('Direct effects' ~ ~ '('*italic(frac(delta*'N'['F']^'*', delta*'S'[1]))*')')
xldirect2 = expression('Focal species sensitivity' ~ '('*italic(frac('a'[1], 'm'[1]*'q'['F2']))*')')

title1 = expression(italic(frac(delta * 'N'['F']^'*',  delta * 'S'[1]) ~ '='))
title2 = expression(italic(frac('a'[1], 'm'[1]*'q'['F2'])) ~ '=')
title3 = expression( italic( rho )  ~ '=' )

#### Figure 4a revised showing indirect effect strength calculated from equation 17
curve(expr = calcIFX1(x, dfxN1 = 0.5), from = 0, to = 1, n = 100, ylim = c(0, 5), 
      xlab = '', ylab = yl1, lty = 3, lwd = 2, cex.lab = 1) 
mtext( xl, side=1, line=2, cex.lab = 1)
curve(expr = calcIFX1(x, dfxN1 = 2), add = TRUE, lty = 2, lwd = 2)
curve(expr = calcIFX1(x, dfxN1 = 5), add = TRUE, lty = 1)
legend(x = 0.1, y = 4.2, legend=c('0.5', '2', '5'), lty = c(3,2,1), lwd = c(2, 2, 1), cex = 0.9, bty='n', title= title1)

##### Figure 4b revised showing indirect effects as a function of direct effects 
curve(expr = calcIFX1( rho = 0.75, dfxN1 = x), from = 0, to = 5, n = 10, 
      xlab = '', ylab = yl1, lty = 3, lwd = 2, cex.lab = 1, ylim = c(0, 5))
mtext(xldirect, side=1, line=4, cex.lab = 1)
curve( expr = calcIFX1 ( rho = 0.5, dfxN1 = x), add = TRUE, lty = 2, lwd = 2)
curve( expr = calcIFX1 ( rho = 0.25, dfxN1 = x), add = TRUE, lty = 1)
legend(x = 0.4, y = 4.8, legend=c('0.75', '0.5', '0.25'), lty = c(3,2,1), lwd = c(2, 2, 1), cex = 0.9, bty='n', title= title3)

#### Figure 5 indirect effects of changing non-limiting resource in relation to rho
#### different lines show different amounts of focal species' sensitivity to resource 2
fSense = c(0.5, 2, 5)
curve(expr = calcFX2 (rho = x, focalSensitivity= fSense[1]) , from = 0, to = 1, n = 100, xlab = xl, 
      ylab = yl2,  ylim = c(-5, 0), lty = 3, lwd = 2, cex.lab = 1)

curve(expr = calcFX2 (rho = x, focalSensitivity = fSense[2]), add = TRUE, lty = 2, lwd = 2)

curve(expr = calcFX2 (rho = x, focalSensitivity= fSense[3]), add = TRUE, lty = 1, lwd = 1)
legend(x = 0, y= -3, legend= as.character(fSense), lty = c(3,2,1), lwd = c(2, 2, 1), cex = 0.9, bty='n', title= title2)

###### figure 5b 
rhos = c(0.75, 0.5, 0.25)
curve( expr = calcFX2 ( rho = rhos[1], focalSensitivity= x), from = 0, to = 5, 
       n = 10, xlab = '', ylab = yl2, lty = 3, lwd = 2, cex.lab = 1, ylim = c(-5, 0))

curve( expr = calcFX2 ( rho = rhos[2], focalSensitivity= x), add = TRUE, 
       lty = 2, lwd = 2)
curve( expr = calcFX2 ( rho = rhos[3], focalSensitivity= x), add = TRUE, 
       lty = 1)
legend(x = 0, y= -3, legend= as.character(rhos), lty = c(3,2,1), lwd = c(2, 2, 1), 
       cex = 0.9, bty='n', title= title3)
mtext(xldirect2, side=1, line=4)

########################################################
####### End Section 1: essential resource model anaylsis 
########################################################

########################################################

########################################################
####### Section 2: Substitutable resource model analysis 
########################################################

rm(list = ls()) #### remove previous parameters 

#### Substitutable resource model Functions 

source('resource_competition_functions.R')

par(mfrow= c(1, 1) ) 

#### set parameters
S = c(1000, 800)

parms = list(r  = c(0.8, 0.8),              #### growth rates
             k = c(152, 150),               #### resource half saturation for resources one and two
             tau  = c(0.6, 0.5),            #### minimum resource requirement 
             D  = 0.7,                      #### 
             q = matrix(c(2 , 5,
                          5  ,2), nrow = 2, byrow = T),    #### sp. one and two per capita use of resource one and two
             w = matrix(c(1, 3, 
                          3, 1) , nrow = 2, byrow = T),
             S = as.numeric(S))                              #### resource supply rates 

B = calcIntercepts(parms)  #### Calculate Zero Net Growth Isocline y-intercepts

##### Plot Zero-net growth isoclines and consumption vectors and supply point 
plotZingis(parms, B1 = B[[1]], B2 = B[[2]])
Rstars = calcRstars(parms, B1 = B[[1]], B2 = B[[2]])
plotConsumption(parms, B1 = B[[1]], B2 = B[[2]], R1star= Rstars[[1]], R2star = Rstars[[2]])
points(parms$S[1], parms$S[2])

##### Run model 
NR = c(N = c(1, 1), R = c(4000,4000))   ### initial populations and resource levels
NRout = ode(y = NR, func= subMod, t = 1:2000, parms = parms, method = 'rk4' )
preEq = as.matrix(tail(NRout, 1))[,-1]

#### show that calculated Rstars match equilibrium resources in simulation: 
Rstars
preEq[ 3:4]

#### change resource supply 
newParms = parms
newParms$S = parms$S*c(1.10, 1) #### increases supply of resource one by 1%

##### Re-run model to equilibrium starting at old equilibrium
NRout2 = ode(y = preEq, func = subMod, t = 1:2000, parms = newParms, method = 'rk4')
Eq2 = as.matrix(tail(NRout2, 1))[,-1]
Eq2

##### find net change in N1
Eq2[1] - preEq[1]

#### Sensitivity of N1 to change in S1 (partial derivative)
simNetSensitivity = (Eq2[1] - preEq[1])/(newParms$S[1] - parms$S[1])
simNetSensitivity

##### Re-run model with constant competitor abundance starting from first equilibrium
NRoutDirect = ode(y = preEq, func = submod_constant_N2, t = 1:2000, parms = newParms, method = 'rk4')
dirEq = as.matrix(tail(NRoutDirect, 1))[, -1]
dirEq[1] - preEq[1] #### change in N1 due to direct effects

points(dirEq[3], dirEq[4], pch = 20) #### show that R1*, R2* shifts a tiny bit when competitor is held constant
text( dirEq[3], dirEq[4], 'new Eq', cex= 0.8, pos = 4)

###### plot population time series 
par(mfrow = c(3,1))
matplot(NRout[ 1:1000, 2:3] , type = 'l', xlab = 'Time', ylim = c(0, 110), 
        ylab = 'Population Density', main = 'first equilibrium') 

##### Supply point 2 
matplot(NRout2[ 1:1000, 2:3], type = 'l', xlab = 'Time', ylab = 'Population Density', 
        ylim = c(0, 110), main = 'net effects')
matplot(NRoutDirect[ 1:1000, 2:3], type = 'l', xlab = 'Time', ylab = 'Population Density', 
        ylim = c(0, 110), main = 'direct effects')

##### Sensitivity of N1 to direct effects of change in S1
simDirSensitivity = (dirEq[1] - preEq[1])/(newParms$S[1] - parms$S[1])
simDirSensitivity

##### Sensitivity of N1 to indirect effects is net minus direct sensitivity
simIndSensitivity = simNetSensitivity - simDirSensitivity
simIndSensitivity

##### Use Lotka-Volterra Equivalents to calculate 
##### net, direct and indirect sensitivity to change in S1
rho = calcRhoSub( parms ) #### find rho based on model parameters
rho

LVP = calcLVparms(parms = parms, B = c(B[[1]], B[[2]]))
LVP

calcDirSens = calcDirSensitivity(parms=parms)

#### check that calculated direct sensitivity matches simulated 
cbind(calculated = calcDirSens[1], simulated = simDirSensitivity)

#### calc net sensitivity 
calcNetSens = calcNetSensitivity(dirN1 = calcDirSens[1], dirN2 = calcDirSens[2], lvp= LVP)
calcNetSens

calcNetSens2 = calcNetSensitivity2 ( parms)
calcNetSens2

#### check that calculated net sensitivity matches simulated
cbind(calculated = calcNetSens, simulated = simNetSensitivity)

#### calculate indirect sensitivity 
calcIndSens = calcNetSens - calcDirSens
calcIndSens[1]

#### check that calc indirect matches simulated indirect 
cbind(calculated = calcIndSens[1], simulated = simIndSensitivity)

#### calculated indirect sensitivity using a function with rho and LV parms
calcIFX(dirN1 = calcDirSens[1], dirN2 = calcDirSens[2], Beta=LVP$betaT, rho= rho)

#### calculate indirect sensitivity directly from mechanistic parameters
calcIFX2(dirN1= calcDirSens[1], dirN2 = calcDirSens[2], parms= parms, rho=rho)

###### Run across a range of parameters and rho's and calculate indirect effects 

parmsList = list(NA)
parmsList[[1]] = parms
rhoGradient = NA
rhoGradient[1] = calcRhoSub(parms)
outputDF = data.frame(simNetFX= NA, simDirFX= NA, simIndFX = NA, calcNetSens = NA, 
                      calcDirN1 = NA, calcDirN2 = NA, calcIndSens = NA, 
                      calcIndSens2 = NA, rho = NA, Beta = NA)

outputDF[1, 1:3] = simulatedEffects(parms)
outputDF[1, 4:10] = analyticalCalculations(parms) 

for (i in 2:10 ){ 
  newParms = parms
  newParms$q = parms$q + matrix( c(0, i, i, 0), nrow = 2, byrow= TRUE)
  parmsList[[i]] = newParms  
  outputDF[i, 1:3 ] = simulatedEffects(newParms)
  outputDF[i, 4:10 ] = analyticalCalculations(newParms)  
}

##### check that simulated and calculated values match
##### columns should be ~ equal (some error due to rounding etc.)
cbind( outputDF$simNetFX , outputDF$calcNetSens ) 
cbind( outputDF$simDirFX, outputDF$calcDirN1)
cbind( outputDF$simIndFX, outputDF$calcIndSens, outputDF$calcIndSens2)

round( outputDF$simNetFX, 5 )  == round(  outputDF$calcNetSens , 5)
round( outputDF$simDirFX, 5 )  == round(  outputDF$calcDirN1, 5 )
round( outputDF$simIndFX, 5 )  == round(  outputDF$calcIndSens, 5 )


#####
##### Generate figure 6 for main text of manuscript
#####

par( mfrow = c(1,1), las=1, oma = c(1,1,1,1), mar = c(5,6,1,1))
yl = 'Indirect effect strength'
yl2 = expression("Indirect effect" ~ "("*italic(frac(delta*'N'['F']^'*', delta*'S'[1])*")"))

curve(expr = calcIFX (dirN1 = 2, dirN2 = 0, rho = x, Beta = 1) , 
      from = 0, to = 1, n = 100, xlab = expression(  rho  ), ylab = yl2,  
      ylim = c(-10, 10), lty = 1, lwd = 1)
curve(expr = calcIFX (dirN1 = 2, dirN2 = 1.33333333, rho = x, Beta = 1) , 
      add = TRUE, lty = 2, lwd = 1) 
curve(expr = calcIFX (dirN1 = 2, dirN2 = 2, rho = x, Beta = 1) , 
      add = TRUE, lty = 1, lwd = 1) 
curve(expr = calcIFX (dirN1 = 1.33333333333, dirN2 = 2, rho = x, Beta = 1) , 
      add = TRUE, lty = 2, lwd = 1) 
curve(expr = calcIFX (dirN1 = 0, dirN2 = 2, rho = x, Beta = 1) , 
      add = TRUE, lty = 1, lwd = 1) 

text(0.85, 9, label = expression( infinity ), cex = 1.3)
text(0.85, -9, label = expression(0), cex = 1)
text(0.99, 6, label = expression( frac(3, 2)), cex = 0.9)
text(0.95, 1, label = expression( 1 ), cex = 1)
text(0.99, -6, label = expression( frac(2,3) ), cex = 0.9)

############################################################
####### End Section 2: Substitutable resource model analysis 
############################################################

