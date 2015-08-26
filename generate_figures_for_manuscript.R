# Generate figures for main text 
# R-code supporting 'Indirect effects of environmental change in resource competition models' 
# Written on R version 3.2.1 (2015-06-18) 
# Author: Andrew R. Kleinhesselink
# Contact: arklein@aggiemail.usu.edu 
# Requires loading the essential_resource_functions.R and the 
# substitutable_resource_functions.R scripts for function definitions 

rm(list = ls())
source('essential_resource_functions.R') 
source('substitutable_resource_functions.R')

# Basic figure settings and labels --------------------------------------------
yl1 <- expression(scriptstyle("Indirect effects") ~ "(" * group( "", 
                                                                scriptstyle( 
                                                                  frac( 
                                                                    italic("dN")["F"]^"*", italic("dS")[1]*' '
                                                                  )
                                                                ), 
                                                                "|")[scriptscriptstyle('indirect')] * ")")


yl2 <- expression(scriptstyle("Indirect effects") ~ "(" * group( "", 
                                                                 scriptstyle( 
                                                                   frac( 
                                                                     italic("dN")["F"]^"*", italic("dS")[2]*' '
                                                                   )
                                                                 ), 
                                                                 "|")[scriptscriptstyle('indirect')] * ")")


xlpdf <- "Niche overlap    "
xldirect <- expression (scriptstyle("Direct effects") ~ "(" * 
                           scriptstyle(
                             frac( italic(partialdiff * "N")["F"]^"*", 
                                   italic(partialdiff * "S")[1]
                               )
                             ) * ")")

xldirect2.1 <- paste ("Focal species'\n   sensitivity" ) 
xldirect2.2 <- expression( "("* frac( italic("a")[1], italic("m")[1]*italic("q")["F2"])*  ")"  )

title1 <- expression(italic(frac(partialdiff * "N"["F"]^"*", partialdiff * "S"[1]) ~ "="))
title2 <- expression( frac(italic("a")[1], italic("m")[1] * italic("q")["F2"]) ~ "=") 

# ------------------------------------------------------------------------------
xlimA <- 1
xlimB <- 5
ylim <- 5

# Figure 3 --------------------------------------------------------------------
pdf(file='figs/fig3.pdf', paper= 'letter', useDingbats=FALSE, pointsize= 12, family='ArialMT' )
par(mfrow = c(1, 1), oma = c(15, 1, 0.2, 0), mar = c(5, 3.8, 2, 0.5), las = 1, cex.axis = 0.8)
par(mfrow = c(1,2), xpd = FALSE ) 

curve(expr = calc_essential_ifx1(x, dfxN1 = 0.5), 
      from = 0, to = xlimA, n = 100, ylim = c(0, ylim), 
      lty = 3, lwd = 2,  
      xlab = "", ylab = yl1, cex.lab = 1.5, mgp = c(1.25,0.75,0) )

curve(expr = calc_essential_ifx1(x, dfxN1 = 2), add = TRUE, lty = 2, lwd = 2)
curve(expr = calc_essential_ifx1(x, dfxN1 = 5), add = TRUE, lty = 1)

legend(x = 0.01*xlimA, y = 0.9*ylim, legend = c("0.5", "2", "5"), lty = c(3, 2, 1), lwd = c(2, 2, 1), 
       cex = 0.85, bty = "n", title = title1)

mtext("A)", side =3, line = 0.1, adj = 0)
par(xpd = TRUE)
mtext("Figure 3", side = 1, line = 0, adj = 0, outer = TRUE)
text( 0.45, -1.51, xlpdf, cex = 1.05 ) 
text( 0.70, -1.51, "(  )")
text(0.70, -1.46, "\\*r", vfont = c('serif symbol', 'italic'), cex= 1)

# 3B --------------------------------------------------------------------------
par(xpd = FALSE, mar = c(5,0.5,2,3.8))

curve(expr = calc_essential_ifx1(rho = 0.75, dfxN1 = x), from = 0, to = xlimB, n = 10, 
      xlab = "", ylab = "", lty = 3, lwd = 2, 
      ylim = c(0, ylim), 
      cex.lab = 1.5, 
      mgp = c(1.25, 0.75, 0), 
      axes = FALSE, 
      frame.plot = TRUE )

Axis(side = 1, labels = TRUE, mgp = c(1.25, 0.75, 0))
Axis(side = 2, labels = FALSE)

curve(expr = calc_essential_ifx1(rho = 0.5, dfxN1 = x), add = TRUE, lty = 2, lwd = 2)
curve(expr = calc_essential_ifx1(rho = 0.25, dfxN1 = x), add = TRUE, lty = 1)

legend(x = 0.01*xlimB, y = 0.9*ylim, legend = c("0.75", "0.5", "0.25"), lty = c(3, 2, 1), lwd = c(2, 2, 1), 
       cex = 0.85, bty = "n", title = "=")

mtext(xldirect, side = 1, line = 3, cex = 1.5)
text(0.6, 4.25, "\\*r", vfont = c('serif symbol', 'italic'), cex= 1)
mtext("B)", side = 3, line = 0.1, adj = 0 )

dev.off()

# Figure 4  -------------------------------------------------------------------
pdf(file='figs/fig4.pdf', paper= 'letter', useDingbats=FALSE, pointsize= 12, family='ArialMT' )
par(mfrow = c(1, 1), oma = c(15, 1, 0.2, 0), mar = c(5, 3.8, 2, 0.5), las = 1, cex.axis = 0.8)
par(mfrow = c(1,2), xpd = FALSE ) 

curve(expr = calc_essential_ifx2(x, focalSensitivity = 0.5), 
      from = 0, to = xlimA, n = 100, ylim = c(-ylim, 0), 
      lty = 3, lwd = 2,  
      xlab = "", ylab = yl2, cex.lab = 1.5, mgp = c(1.7,0.75,0) )

curve(expr = calc_essential_ifx2(x, focalSensitivity= 2), add = TRUE, lty = 2, lwd = 2)
curve(expr = calc_essential_ifx2(x, focalSensitivity= 5), add = TRUE, lty = 1)

legend(x = 0.01*xlimA, y = -0.4*ylim, legend = c("0.5", "2", "5"), lty = c(3, 2, 1), lwd = c(2, 2, 1), 
       cex = 0.85, bty = "n", title = title2)

mtext("A)", side =3, line = 0.1, adj = 0)
par(xpd = TRUE)
mtext("Figure 4", side = 1, line = 0, adj = 0, outer = TRUE)
text( 0.45, -4.82 + -1.51, xlpdf, cex = 1.05 ) 
text( 0.70, -4.82 + -1.51, "(  )", cex = 1.05)
text(0.70, -4.82 + -1.46, "\\*r", vfont = c('serif symbol', 'italic'), cex= 1)

# 4B --------------------------------------------------------------------------
par(xpd = FALSE, mar = c(5,0.5,2,3.8))

curve(expr = calc_essential_ifx2(rho = 0.75,focalSensitivity = x), from = 0, to = xlimB, n = 10, 
      xlab = "", ylab = "", lty = 3, lwd = 2, 
      ylim = c(-ylim, 0), 
      cex.lab = 1.5, 
      mgp = c(1.25, 0.75, 0), 
      axes = FALSE, 
      frame.plot = TRUE )

Axis(side = 1, labels = TRUE, mgp = c(1.25, 0.75, 0))
Axis(side = 2, labels = FALSE)

curve(expr = calc_essential_ifx2(rho = 0.5, focalSensitivity = x), add = TRUE, lty = 2, lwd = 2)
curve(expr = calc_essential_ifx2(rho = 0.25, focalSensitivity = x), add = TRUE, lty = 1)

legend(x = 0.01*xlimB, y = -0.4*ylim, legend = c("0.75", "0.5", "0.25"), lty = c(3, 2, 1), lwd = c(2, 2, 1), 
       cex = 0.85, bty = "n", title = "=")

par(xpd = TRUE)
text(0.6, -2.25, "\\*r", vfont = c('serif symbol', 'italic'), cex= 1)
mtext("B)", side = 3, line = 0.1, adj = 0 )
text( 1.5, -4.8 + -1.71, xldirect2.1, cex = 1.05 ) 
text( 3.5, -4.8 + -1.71, xldirect2.2, cex = 0.9)

dev.off()



## Figure 5 Indirect effects as a function of rho in a substitutable resource model 
pdf(file='figs/fig5.pdf', paper= 'letter', useDingbats=FALSE, pointsize= 12, family='ArialMT' )
par(mfrow = c(1, 1), oma = c(15, 1, 0.2, 0), mar = c(5, 3.8, 2, 0.5), las = 1, cex.axis = 0.8)
par(mfrow = c(1,2), xpd = FALSE ) 



#par(mfrow = c(1, 2), oma = c(15, 2, 0.2, 13), mar = c(5, 3.8, 2, 0.5), las = 1, cex.axis = 0.8)
#par(xpd = FALSE ) 

curve(expr = calc_substitutable_ifx(dirN1 = 4, dirN2 = 0, rho = x, Beta = 1), 
      from = 0, to = 1, n = 100, xlab = '', 
      ylab = yl1, 
      cex.lab = 1.5, 
      mgp = c(1.5, 0.75, 0),
      ylim = c(-10, 10), 
      lty = 1, lwd = 1)

curve(expr = calc_substitutable_ifx(dirN1 = 2, dirN2 = 1, rho = x, Beta = 1), add = TRUE, lty = 2, lwd = 1)
curve(expr = calc_substitutable_ifx(dirN1 = 1, dirN2 = 1, rho = x, Beta = 1), add = TRUE, lty = 1, lwd = 1)
curve(expr = calc_substitutable_ifx(dirN1 = 1, dirN2 = 2, rho = x, Beta = 1), add = TRUE, lty = 2, lwd = 1)
curve(expr = calc_substitutable_ifx(dirN1 = 0, dirN2 = 4, rho = x, Beta = 1), add = TRUE, lty = 1, lwd = 1)
text(0.75, 9, label = expression(+4), cex = 0.9)
text(0.75, -9, label = expression(-4), cex = 0.9)
text(0.99, 6.6, label = expression(+1), cex = 0.9)
text(0.99, 1, label = expression(0), cex = 0.9)
text(0.99, -6.6, label = expression(-1), cex = 0.9)
par(xpd = TRUE)
mtext("Figure 5", side = 1, line = 0, adj = 0, outer = TRUE)

text( 0.48, -14 + -1.5, xlpdf, cex = 1.05 ) 
text( 0.73, -14 + -1.5, "(  )", cex = 1.05)
text(0.73, -14 + -1.42, "\\*r", vfont = c('serif symbol', 'italic'), cex= 1)

dev.off()
