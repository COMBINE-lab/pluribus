rel_diff_plot <- function(d, true, x, y) {
  rdx <- rel_diff(d, true, x)
  rdy <- rel_diff(d, true, y)
  dt <- data.table(true=d[[true]], abs_diff=abs(rdx) - abs(rdy))
  l1 <- unlist(strsplit(x, "\\."))[2]
  l2 <- unlist(strsplit(y, "\\."))[2]
  ggplot(dt, aes(true, abs_diff)) + stat_binhex(bins=50) + 
    scale_x_log10(na.value=0) + 
    ylab(bquote("ARD"[.(l1)] ~ " - " ~ "ARD"[.(l2)])) + 
    xlab("True read count") + theme_minimal()
}

scatter_plot <- function(d, x, y) {
  ggplot(m, aes_string(x, y)) + geom_point() + scale_x_log10(limits=c(1e-3,NA)) + 
    scale_y_log10(limits=c(1e-3,NA)) + 
    theme_minimal() + stat_binhex(bins=90) + 
    scale_fill_gradientn(colors=c("blue", "green", "red"), trans="log10")
}