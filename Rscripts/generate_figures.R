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
par(mfrow = c(1, 1), oma = c(2, 2, 2, 2), mar = c(5, 8, 2, 2), las = 1, cex.axis = 0.8)

rho_UTF <- '\u03C1' # UTF code for Greek letter rho for use in plots 
yl1 <- expression("Indirect effects" ~ ~~italic(frac("dN"["F indirect"]^"*", "dS"[1])))
yl2 <- expression("Indirect effects" ~ ~~italic(frac("dN"["F indirect"]^"*", "dS"[2])))
xl <- substitute( "Niche overlap" ~ "(" * italic(rho_UTF) * ")", list(rho_UTF = rho_UTF))

xldirect <- expression("Direct effects" ~ ~"(" * italic(frac(partialdiff * "N"["F"]^"*", partialdiff * "S"[1])) * ")")
xldirect2 <- expression("Focal species sensitivity" ~ "(" * italic(frac("a"[1], "m"[1] * "q"["F2"])) * ")")

title1 <- expression(italic(frac(partialdiff * "N"["F"]^"*", partialdiff * "S"[1]) ~ "="))
title2 <- expression(italic(frac("a"[1], "m"[1] * "q"["F2"])) ~ "=")
title3 <- as.expression( substitute("" * italic(rho_UTF) ~ "=", list( rho_UTF = rho_UTF)))
# ------------------------------------------------------------------------------

# Figure 3A Indirect effect strength as a function of niche overlap ------------
curve(expr = calc_essential_ifx1(x, dfxN1 = 0.5), from = 0, to = 1, n = 100, ylim = c(0, 5), 
      xlab = "", ylab = yl1, lty = 3, lwd = 2, cex.lab = 1)
mtext(xl, side = 1, line = 3, cex.lab = 1)
curve(expr = calc_essential_ifx1(x, dfxN1 = 2), add = TRUE, lty = 2, lwd = 2)
curve(expr = calc_essential_ifx1(x, dfxN1 = 5), add = TRUE, lty = 1)
legend(x = 0.1, y = 4.2, legend = c("0.5", "2", "5"), lty = c(3, 2, 1), lwd = c(2, 2, 1), 
       cex = 0.9, bty = "n", title = title1)
mtext("Figure 3A", side =3, line = 1 )


# Figure 3B Indirect effect strength as a function of direct effects ------------ 
curve(expr = calc_essential_ifx1(rho = 0.75, dfxN1 = x), from = 0, to = 5, n = 10, 
      xlab = "", ylab = yl1, lty = 3, lwd = 2, cex.lab = 1, ylim = c(0, 5))
mtext(xldirect, side = 1, line = 4, cex.lab = 1)
curve(expr = calc_essential_ifx1(rho = 0.5, dfxN1 = x), add = TRUE, lty = 2, lwd = 2)
curve(expr = calc_essential_ifx1(rho = 0.25, dfxN1 = x), add = TRUE, lty = 1)
legend(x = 0.4, y = 4.8, legend = c("0.75", "0.5", "0.25"), lty = c(3, 2, 1), lwd = c(2, 2, 1), cex = 0.9, bty = "n", title = title3)
mtext("Figure 3B", side =3, line = 1 )

# Figure 4A Indirect effects of change in non-limiting resource as a function of rho. 
fSense <- c(0.5, 2, 5)
curve(expr = calc_essential_ifx2(rho = x, focalSensitivity = fSense[1]), from = 0, to = 1, n = 100, 
      xlab = '', ylab = yl2, ylim = c(-5, 0), lty = 3, lwd = 2, cex.lab = 1)
mtext(xl, side = 1, line = 3, cex.lab = 1)
curve(expr = calc_essential_ifx2(rho = x, focalSensitivity = fSense[2]), add = TRUE, lty = 2, lwd = 2)
curve(expr = calc_essential_ifx2(rho = x, focalSensitivity = fSense[3]), add = TRUE, lty = 1, lwd = 1)
legend(x = 0, y = -3, legend = as.character(fSense), lty = c(3, 2, 1), lwd = c(2, 2, 1), cex = 0.9, 
       bty = "n", title = title2)
mtext("Figure 4A", side =3, line = 1 )

# Figure 4B Indirect effects of a change in non-limiting resource as a function of 
# direct effects on the focal species. 
rhos <- c(0.75, 0.5, 0.25)
curve(expr = calc_essential_ifx2(rho = rhos[1], focalSensitivity = x), from = 0, to = 5, n = 10, xlab = "", ylab = yl2, lty = 3, lwd = 2, cex.lab = 1, ylim = c(-5, 0))
curve(expr = calc_essential_ifx2(rho = rhos[2], focalSensitivity = x), add = TRUE, lty = 2, lwd = 2)
curve(expr = calc_essential_ifx2(rho = rhos[3], focalSensitivity = x), add = TRUE, lty = 1)
legend(x = 0, y = -3, legend = as.character(rhos), lty = c(3, 2, 1), lwd = c(2, 2, 1), cex = 0.9, bty = "n", title = title3)
mtext(xldirect2, side = 1, line = 4, cex.lab = 1) 
mtext("Figure 4B", side =3, line = 1 )

## Figure 5 Indirect effects as a function of rho in a substitutable resource model 
curve(expr = calc_substitutable_ifx(dirN1 = 4, dirN2 = 0, rho = x, Beta = 1), 
      from = 0, to = 1, n = 100, xlab = '', ylab = yl1, ylim = c(-10, 10), lty = 1, lwd = 1, cex.lab = 1)
mtext(xl, side = 1, line = 3, cex.lab = 1)
curve(expr = calc_substitutable_ifx(dirN1 = 2, dirN2 = 1, rho = x, Beta = 1), add = TRUE, lty = 2, lwd = 1)
curve(expr = calc_substitutable_ifx(dirN1 = 1, dirN2 = 1, rho = x, Beta = 1), add = TRUE, lty = 1, lwd = 1)
curve(expr = calc_substitutable_ifx(dirN1 = 1, dirN2 = 2, rho = x, Beta = 1), add = TRUE, lty = 2, lwd = 1)
curve(expr = calc_substitutable_ifx(dirN1 = 0, dirN2 = 4, rho = x, Beta = 1), add = TRUE, lty = 1, lwd = 1)
text(0.75, 9, label = expression(+4), cex = 0.9)
text(0.75, -9, label = expression(-4), cex = 0.9)
text(0.99, 6.6, label = expression(+1), cex = 0.9)
text(0.99, 1, label = expression(0), cex = 0.9)
text(0.99, -6.6, label = expression(-1), cex = 0.9)
mtext("Figure 5", side =3, line = 1 )
