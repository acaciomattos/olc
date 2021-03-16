#' Optimal Level collapsing apply
#'
#' @param YY Vector with response variable.
#' @param Levels Vector with the categorical variable to be collapsed.
#' @param model The model distribution to be used. See ?olc.eval for available model distributions.
#' @param k The number of groups to collapse the categorical variable.
#'
#' @return A list with the corresponded collapsed groups for different methods.
#' @export
#'
#' @examples
olc <- function(YY,Levels,model,k){
  options(scipen=999)
  if(!class(Levels) %in% "factor"){Levels<-as.factor(Levels)
  warning("The nominal variable was coerced to factor vector")}
  
  if(!model %in% c("binomial","gaussian","Gamma","inverse.gaussian",
                   "poisson","quasi","quasibinomial","quasipoisson",
                   "nb","poi_zero","neg_zero","geo_zero")){stop("The choosen family is not available.")}
  
  
  
  Data_aux <- data.frame(YY,Levels)
  if(model=="nb"){modelo1 <- glm.nb(YY ~-1+Levels,data=Data_aux,control = list(maxit=50))}
  if(model=="poi_zero"){modelo1 <-  zeroinfl(YY~-1+Levels | 1, data=Data_aux, dist = "poisson")}
  if(model=="neg_zero"){modelo1 <-  zeroinfl(YY~-1+Levels | 1, data=Data_aux, dist = "negbin")}
  if(model=="geo_zero"){modelo1 <-  zeroinfl(YY~-1+Levels | 1, data=Data_aux, dist = "geometric")}
  if(model=="gaussian"){modelo1 <-  lm(YY~-1+Levels,data=Data_aux)}
  if(model%in%c("binomial","gaussian","Gamma","inverse.gaussian",
                "poisson","quasi","quasibinomial","quasipoisson")){modelo1 <-  glm(YY ~-1+Levels,data=Data_aux,family=model,control = list(maxit=50))}
  modelo2 <-  lm(YY~-1+Levels,data=Data_aux)
  
  
  #### Level collapsing for Levelss with Poisson ####
  if(model %in% c("poisson","binomial","gaussian","quasipoisson","quasi","quasibinomial","inverse.gaussian","Gamma","nb")){
    variable_aux <- data.frame(
      Beta_mean = modelo1$coefficients,
      Beta_low = ifelse(summary(modelo1)[[12]][,4] > 0.1,0,modelo1$coefficients-1.96*summary(modelo1)[[12]][,2]),
      Beta_upp = ifelse(summary(modelo1)[[12]][,4] > 0.1,0,modelo1$coefficients+1.96*summary(modelo1)[[12]][,2]),
      Levels = modelo1$xlevels
    )}
  if(model %in% c("poi_zero","neg_zero","geo_zero")){
    variable_aux <- data.frame(
      Beta_mean = modelo1$coefficients$count,
      Beta_low = ifelse(summary(modelo1)$coefficients$count[,4] > 0.1,0,modelo1$coefficients$count-1.96*summary(modelo1)$coefficients$count[,2]),
      Beta_upp = ifelse(summary(modelo1)$coefficients$count[,4] > 0.1,0,modelo1$coefficients$count+1.96*summary(modelo1)$coefficients$count[,2]),
      Levels = levels(as.factor(Data_aux$Levels)))
  }
  
  options(warn=-1)
  #Clustering
  clust_aux1  <- hclust(dist(variable_aux$Beta_mean),method="ward.D2")
  clust_aux2  <- hclust(dist(variable_aux[,c("Beta_low","Beta_upp")]),method="ward.D2")
  
  
  variable_aux$centroids1 <- cutree((clust_aux1),k=k)
  variable_aux$centroids2 <- cutree((clust_aux2),k=k)
  
  centroids_aux1 <- variable_aux %>% group_by(centroids1) %>% summarise(cent=mean(Beta_mean)) %>% arrange(cent)
  centroids_aux2 <- variable_aux %>% group_by(centroids2) %>% summarise(cent1=mean(Beta_low),cent2=mean(Beta_upp)) %>% arrange(cent1,cent2)
  
  variable_aux$centroids1 <- NULL
  variable_aux$centroids2 <- NULL
  
  variable_aux$cluster.met1 <- kmeans(variable_aux$Beta_mean,centers=unique(centroids_aux1$cent),iter.max = 20)$cluster
  variable_aux$cluster.met2 <- kmeans(variable_aux[,c("Beta_low","Beta_upp")],centers=unique(cbind(centroids_aux2$cent1,centroids_aux2$cent2)),iter.max = 20)$cluster
  
  #Standardizing new variable name and setting as factor
  names(variable_aux)[which(names(variable_aux)=="cluster.met1")] <- paste0("olc.met1.",k,".",names(variable_aux)[4])
  names(variable_aux)[which(names(variable_aux)=="cluster.met2")] <- paste0("olc.met2.",k,".",names(variable_aux)[4])
  Data_aux <- left_join(Data_aux,variable_aux[,-c(1:3)],by=c("Levels"="Levels"))
  
  scheffe_test <- scheffe.test(modelo2,"Levels",group = TRUE,alpha = 0.05)
  
  scheffe_collapsed <- data.frame(
    VAR = as.character(substr(names(sort(lm(YY~-1+Levels,data=Data_aux)$coefficients,decreasing=TRUE)),nchar(names(lm(YY~-1+Levels,data=Data_aux)$xlevels))+1,100)),
    Grupo_scheffe = as.character(scheffe_test$groups$groups)
  )
  
  Data_aux <- left_join(Data_aux,scheffe_collapsed,by=c("Levels"="VAR"))
  
  Data_aux$Grupo_scheffe <- as.character(Data_aux$Grupo_scheffe)
  Data_aux$Grupo_scheffe <- factor(Data_aux$Grupo_scheffe)
  
  snk_test <- SNK.test(modelo2,"Levels",group = TRUE,alpha = 0.1)
  
  snk_collapsed <- data.frame(
    VAR = as.character(substr(names(sort(lm(YY~-1+Levels,data=Data_aux)$coefficients,decreasing=TRUE)),nchar(names(lm(YY~-1+Levels,data=Data_aux)$xlevels))+1,100)),
    Grupo_snk = as.character(snk_test$groups$groups)
  )
  
  
  Data_aux <- left_join(Data_aux,snk_collapsed,by=c("Levels"="VAR"))
  
  OUTPUT <- list(
    data_groups = Data_aux,
    variable_aux = variable_aux[,4:6]
  )
  
  return(OUTPUT)
  options(scipen=0)
}