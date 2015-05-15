##### quick demonstration of how a single species resource equilibrium 
##### changes in a substitutable model 
rm(list = ls())

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
plotZingis(parms, B1 = B[[1]], B2 = NULL)
Rstars = calcRstars(parms, B1 = B[[1]], B2 = B[[2]])
points(parms$S[1], parms$S[2])
text( parms$S[1], parms$S[2], 'S1,S2', pos = 3, cex = 0.8)

##### Run model 
NR = c(N = c(1, 0), #### set compitor to zero
       R = c(4000,4000))   
NRout = ode(y = NR, func= subMod, t = 1:2000, parms = parms, method = 'rk4' )
preEq = as.matrix(tail(NRout, 1))[,-1]
points(preEq[3], preEq[4]) #### plot R1*, R2* for single species 
text( preEq[3], preEq[4], labels='Eq1', pos = 4)
plotConsumption(parms, B1 = B[[1]], B2 = B[[2]], R1star= preEq[[3]], R2star = preEq[[4]], species2Color=NULL)

parms2 = parms
parms2$S = parms$S*c(0.80, 1)  #### change S1
points(parms2$S[1], parms2$S[2])
text( parms2$S[1], parms2$S[2], 'New\nS1,S2', pos = 2, cex = 0.8)

NR = c(N = c(1, 0), R = c(4000,4000))   ### initial populations and resource levels
NRout = ode(y = NR, func= subMod, t = 1:2000, parms = parms2, method = 'rk4' )
preEq = as.matrix(tail(NRout, 1))[,-1]
points(preEq[3], preEq[4], pch = 19)
text( preEq[3], preEq[4], labels='Eq2', pos = 3)
plotConsumption(parms, B1 = B[[1]], B2 = B[[2]], R1star= preEq[[3]], R2star = preEq[[4]], species2Color=NULL)






