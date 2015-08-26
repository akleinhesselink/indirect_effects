install.packages('scatterplot3d')
require(scatterplot3d)
require(ggplot2)


curve(expr = calc_substitutable_ifx(dirN1 = 2, dirN2 = 1, rho = x, Beta = 1), add = TRUE, lty = 2, lwd = 1)

predict_df <-  expand.grid (  dirN1 = 0, dirN2 = seq(0, -10, -0.1), Beta = 1, rho = seq(0, 0.9, 0.05), ifx = NA) 

for ( i in 1:nrow (predict_df)  ) {   
  predict_df$ifx[i] <- calc_substitutable_ifx( predict_df$dirN1[i], predict_df$dirN2[i], predict_df$Beta[i], predict_df$rho[i]  )
}


ifx3d <- ggplot(predict_df, aes(x = rho, y = dirN2, z= log(ifx)))+
  geom_tile(aes(fill= log(ifx)))+
  stat_contour(bins=10,aes(x = rho,y = dirN2,z= log(ifx)), color="black", size=0.6) + 
  scale_fill_continuous( low= 'blue', high = 'red', guide_legend( title= "log indirect effects")) + 
  geom_vline( xintercept = 0.25, linetype = 'dashed') + 
  xlab( 'nich overlap') + ylab( 'direct effects on competitor')

ifx3d

pdf( file = "~/Dropbox/peer_reviews/ifx3d.pdf", width= 7, height = 5) 
ifx3d 
dev.off()

