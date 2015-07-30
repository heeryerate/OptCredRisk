res <- list()
for(i in 1:2){
  lossprob <- lossdata[[i]]
  lossprob$valueprob[is.infinite(lossprob$valueprob)] <- 0
  lossprob$nameprob <- c(rep("Isminus",ind),rep("Isplus",ind),rep("Mcminus",ind),rep("Mcplus",ind))
  lossprob$Simulation <- factor(lossprob$colorprob, levels = c(1,2), labels = c("IS", "MC"))
  a <- ggplot2::ggplot(lossprob, ggplot2::aes(x = loss, y = valueprob, color = Simulation, shape = nameprob))
  a <- a + ggplot2::scale_x_continuous(name = "Loss Level" , limits = c(10,140)) + ggplot2::scale_y_continuous(name = "log(Tail Probability)", limits = c(-7.5,-2))
  # a <- a + ggplot2::geom_point(size = 3)
  a <- a + ggplot2::geom_line(size = 1.2)
  res[[i]] <- a
}
gridExtra::grid.arrange(res[[1]],res[[2]],res[[1]],res[[2]])
# gridExtra::grid.arrange(res[[5]],res[[6]],res[[7]],res[[8]])