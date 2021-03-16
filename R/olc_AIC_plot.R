#' Optimal Level Collapsing AIC plot
#'
#' @param object object created from olc.eval() function.
#'
#' @return A scatter plot with the AIC by number of collapsed groups from different methods.
#' @export
#'
#' @examples
olc_AIC_plot <- function(object){
  print(plot(object$aic$k,object$aic$Method1,type = "l",
             xlab="number of collapsed groups",ylab="AIC",
             main=paste0("AIC of the model by number of collapsed groups"),
             ylim = c(min(c(object$aic$Method1, object$aic$Method2)), max(c(object$aic$Method1, object$aic$Method2)))))
  print(lines(object$aic$k,object$aic$Method2,col="red"))
  print(points(object$scheffe$k,object$scheffe$aic,col="brown"))
  print(abline(v=object$scheffe$k,col="brown"))
  print(points(object$snk$k,object$snk$aic,col="orange"))
  print(abline(v=object$snk$k,col="orange"))
  print(legend("right",legend=c("Method 1","Method 2","Scheffe","SNK"),col=c("black","red","brown","orange"),lty = 1,cex=0.6))
}