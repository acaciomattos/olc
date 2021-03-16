#' Optimal Level Collapsing LR-test p-value plot
#'
#' @param object object created from olc.eval() function.
#'
#' @return A scatter plot with the LR-test p-value by number of collapsed groups from different methods.
#' @export
#'
#' @examples
olc_LR_plot <- function(object){
  print(plot(object$lr_p$k,object$lr_p$Method1,type = "l",
             xlab="number of collapsed groups",ylab="LR-test p-value",
             main=paste0("LR-test p-value of the model against the full model by number of collapsed groups"),
             ylim = c(min(c(object$lr_p$Method1, object$lr_p$Method2)), max(c(object$lr_p$Method1, object$lr_p$Method2)))))
  print(lines(object$lr_p$k,object$lr_p$Method2,col="red"))
  print(points(object$scheffe$k,object$scheffe$lr_p,col="brown"))
  print(abline(v=object$scheffe$k,col="brown"))
  print(points(object$snk$k,object$snk$lr_p,col="orange"))
  print(abline(v=object$snk$k,col="orange"))
  print(legend("right",legend=c("Method 1","Method 2","Scheffe","SNK"),col=c("black","red","brown","orange"),lty = 1,cex=0.6))
}