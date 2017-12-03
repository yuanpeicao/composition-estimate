require(ggplot2)
require(gridExtra)

hist_top <- ggplot() + geom_histogram(aes(rnorm(100)))
empty <- ggplot() + geom_point(aes(1,1), colour = "white")+
  theme(axis.ticks = element_blank(), 
        panel.background = element_blank(), 
        axis.text.x = element_blank(), axis.text.y = element_blank(),           
        axis.title.x = element_blank(), axis.title.y = element_blank())

scatter <- ggplot() + geom_point(aes(rnorm(100), rnorm(100)))
hist_right <- ggplot() + geom_histogram(aes(rnorm(100))) + coord_flip()

grid.arrange(hist_top, empty, scatter, hist_right, ncol = 2, 
             nrow = 2, widths = c(4, 1), heights = c(1, 4))