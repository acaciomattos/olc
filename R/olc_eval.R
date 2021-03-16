#' Optimal Level Collapsing Evaluation
#'
#' @param YY Vector with response variable.
#' @param Levels Vector with the categorical variable to be collapsed.
#' @param model The model distribution to be used. Available families are: "binomial", "gaussian", "Gamma", "inverse.gaussian", "poisson", "quasi", "quasibinomial", "quasipoisson", "nb" (Negative Binomial), "poi_zero" (Zero Inflated Poisson), "neg_zero" (Zero Inflated Negative Binomial), "geo_zero" (Zero Inflated Geometric).
#' @param k.max The maximum number of groups to be evaluated.
#'
#' @return A list with the corresponded collapsed groups for different methods and number of collapsed groups.
#' @export
#'
#' @examples
olc.eval <- function(YY,Levels,model,k.max){
  
  # browser()
  
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
  
  deviance_aux_met1 <- c()
  deviance_aux_met2 <- c()
  aic_aux_met1 <- c()
  aic_aux_met2 <- c()
  lr_p_v_met1 <- c()
  lr_p_v_met2 <- c()
  k.max <- ifelse(is.null(k.max[1]),
                  ifelse(length(unique(variable_aux$Beta_mean))>70,70,length(unique(variable_aux$Beta_mean))),k.max)
  
  
  for(ii in 2:k.max){
    options(warn=-1)
    #Clustering
    clust_aux1  <- hclust(dist(variable_aux$Beta_mean),method="ward.D2")
    clust_aux2  <- hclust(dist(variable_aux[,c("Beta_low","Beta_upp")]),method="ward.D2")
    
    
    variable_aux$centroids1 <- cutree((clust_aux1),k=ii)
    variable_aux$centroids2 <- cutree((clust_aux2),k=ii)
    
    centroids_aux1 <- variable_aux %>% group_by(centroids1) %>% summarise(cent=mean(Beta_mean)) %>% arrange(cent)
    centroids_aux2 <- variable_aux %>% group_by(centroids2) %>% summarise(cent1=mean(Beta_low),cent2=mean(Beta_upp)) %>% arrange(cent1,cent2)
    
    variable_aux$centroids1 <- NULL
    variable_aux$centroids2 <- NULL
    
    variable_aux$cluster.met1 <- kmeans(variable_aux$Beta_mean,centers=unique(centroids_aux1$cent),iter.max = 20)$cluster
    variable_aux$cluster.met2 <- kmeans(variable_aux[,c("Beta_low","Beta_upp")],centers=unique(cbind(centroids_aux2$cent1,centroids_aux2$cent2)),iter.max = 20)$cluster
    #aa <- try(kmeans(variable_aux$Beta_mean,centers=unique(centroids_aux4$cent),iter.max = 20)$cluster,silent = TRUE)
    
    #Standardizing new variable name and setting as factor
    names(variable_aux)[which(names(variable_aux)=="cluster.met1")] <- paste0("olc.met1.",ii,".",names(variable_aux)[4])
    names(variable_aux)[which(names(variable_aux)=="cluster.met2")] <- paste0("olc.met2.",ii,".",names(variable_aux)[4])
    Data_aux <- left_join(Data_aux,variable_aux[,-c(1:3)],by=c("Levels"="Levels"))
    Data_aux[,grep(paste0(".",names(variable_aux)[4]),names(Data_aux))] <- apply(Data_aux[,grep(paste0(".",names(variable_aux)[4]),names(Data_aux))],2,as.factor)
    #Creating auxiliar formulas.-.
    formula_aux1 <- as.formula(paste("YY~","-1+",paste0("olc.met1.",ii,".",names(variable_aux)[4])))
    formula_aux2 <- as.formula(paste("YY~","-1+",paste0("olc.met2.",ii,".",names(variable_aux)[4])))
    
    if(model %in% c("poi_zero","neg_zero","geo_zero")){
      formula_aux1 <- as.formula(paste("YY~","-1+",paste0("olc.met1.",ii,".",names(variable_aux)[4]),paste0(" | 1")))
      formula_aux2 <- as.formula(paste("YY~","-1+",paste0("olc.met2.",ii,".",names(variable_aux)[4]),paste0(" | 1")))
    }
    #Fitting models
    if(model=="nb"){
      modelo1.olc <- glm.nb(formula_aux1,data=Data_aux,control = list(maxit=50))
      modelo2.olc <- glm.nb(formula_aux2,data=Data_aux,control = list(maxit=50))
    }
    if(model=="poi_zero"){
      modelo1.olc <- zeroinfl(formula_aux1,data=Data_aux,dist="poisson")
      modelo2.olc <- zeroinfl(formula_aux2,data=Data_aux,dist="poisson")
    }
    
    if(model=="geo_zero"){
      modelo1.olc <- zeroinfl(formula_aux1,data=Data_aux,dist="geometric")
      modelo2.olc <- zeroinfl(formula_aux2,data=Data_aux,dist="geometric")
    }
    
    if(model=="neg_zero"){
      modelo1.olc <- zeroinfl(formula_aux1,data=Data_aux,dist="negbin")
      modelo2.olc <- zeroinfl(formula_aux2,data=Data_aux,dist="negbin")
    }
    
    if(model=="gaussian"){
      modelo1.olc <- lm(formula_aux1,data=Data_aux)
      modelo2.olc <- lm(formula_aux2,data=Data_aux)
    }
    if(model %in% c(c("binomial","gaussian","Gamma","inverse.gaussian",
                      "poisson","quasi","quasibinomial","quasipoisson","nb"))){
      modelo1.olc <- glm(formula_aux1,data=Data_aux,family=model,control = list(maxit=50))
      modelo2.olc <- glm(formula_aux2,data=Data_aux,family=model,control = list(maxit=50))
    }
    if(model %in% c("binomial","gaussian","Gamma","inverse.gaussian",
                    "poisson","quasi","quasibinomial","quasipoisson","nb")){
      aic_aux_met1 <- c(aic_aux_met1,modelo1$aic/modelo1.olc$aic)
      aic_aux_met2 <- c(aic_aux_met2,modelo1$aic/modelo2.olc$aic)
    }
    if(model %in% c("poi_zero","neg_zero","geo_zero")){
      aic_aux_met1 <- c(aic_aux_met1,AIC(modelo1)/AIC(modelo1.olc))
      aic_aux_met2 <- c(aic_aux_met2,AIC(modelo1)/AIC(modelo2.olc))
    }
    lr_p_v_met1 <- c(lr_p_v_met1,lrtest(modelo1,modelo1.olc)[2,5])
    lr_p_v_met2 <- c(lr_p_v_met2,lrtest(modelo1,modelo2.olc)[2,5])
    if(model=="binomial"){
      auc_aux1 <- c(auc_aux1,roc(Data_aux$YY2,modelo1.olc$fitted.values)$auc[1])
      auc_aux2 <- c(auc_aux1,roc(Data_aux$YY2,modelo2.olc$fitted.values)$auc[1])
    }
    #CLeaning data
    Data_aux[,grep(paste0(".",names(variable_aux)[4]),names(Data_aux))] <- NULL
    variable_aux[,grep(paste0(".",names(variable_aux)[4]),names(variable_aux))] <- NULL
    options(warn=0)
  }
  
  # Scheffe 
  scheffe_test <- scheffe.test(modelo2,"Levels",group = TRUE,alpha = 0.05)
  
  scheffe_collapsed <- data.frame(
    VAR = as.character(substr(names(sort(lm(YY~-1+Levels,data=Data_aux)$coefficients,decreasing=TRUE)),nchar(names(lm(YY~-1+Levels,data=Data_aux)$xlevels))+1,100)),
    Grupo_scheffe = as.character(scheffe_test$groups$groups)
  )
  
  Data_aux <- left_join(Data_aux,scheffe_collapsed,by=c("Levels"="VAR"))
  
  Data_aux$Grupo_scheffe <- as.character(Data_aux$Grupo_scheffe)
  Data_aux$Grupo_scheffe <- factor(Data_aux$Grupo_scheffe)
  
  if(length(Data_aux$Grupo_scheffe) == 1){
    if(model=="nb"){Modelo1.scheffe <- glm.nb(YY~-1+Grupo_scheffe,data=Data_aux)}
    if(model=="gaussian"){Modelo1.scheffe <- lm(YY~-1+Grupo_scheffe,data=Data_aux)}
    if(model=="poi_zero"){Modelo1.scheffe <- zeroinfl(YY~-1+Grupo_scheffe | 1,data=Data_aux,dist="poisson")}
    if(model=="geo_zero"){Modelo1.scheffe <- zeroinfl(YY~-1+Grupo_scheffe | 1,data=Data_aux,dist="geometric")}
    if(model=="neg_zero"){Modelo1.scheffe <- zeroinfl(YY~-1+Grupo_scheffe | 1,data=Data_aux,dist="negbin")}
    if(model %in% c("binomial","gaussian","Gamma","inverse.gaussian",
                    "poisson","quasi","quasibinomial","quasipoisson")){Modelo1.scheffe <- glm(YY~-1+Grupo_scheffe,data=Data_aux,family=model)}
    
  } else{
    Modelo1.scheffe <- NULL
  }
  
  
  # SNK 
  snk_test <- SNK.test(modelo2,"Levels",group = TRUE,alpha = 0.1)
  
  snk_collapsed <- data.frame(
    VAR = as.character(substr(names(sort(lm(YY~-1+Levels,data=Data_aux)$coefficients,decreasing=TRUE)),nchar(names(lm(YY~-1+Levels,data=Data_aux)$xlevels))+1,100)),
    Grupo_snk = as.character(snk_test$groups$groups)
  )
  
  
  Data_aux <- left_join(Data_aux,snk_collapsed,by=c("Levels"="VAR"))
  
  
  Data_aux$Grupo_snk <- as.character(Data_aux$Grupo_snk)
  Data_aux$Grupo_snk <- factor(Data_aux$Grupo_snk)
  
  if(length(unique(Data_aux$Grupo_snk))!=1){
    if(model=="gaussian"){Modelo1.snk <- lm(YY~-1+Grupo_snk,data=Data_aux)}
    if(model=="nb"){Modelo1.snk <- glm.nb(YY~-1+Grupo_snk,data=Data_aux)}
    if(model=="poi_zero"){Modelo1.snk <- zeroinfl(YY~-1+Grupo_snk | 1,data=Data_aux,dist="poisson")}
    if(model=="geo_zero"){Modelo1.snk <- zeroinfl(YY~-1+Grupo_snk | 1,data=Data_aux,dist="geometric")}
    if(model=="neg_zero"){Modelo1.snk <- zeroinfl(YY~-1+Grupo_snk | 1,data=Data_aux,dist="negbin")}
    if(model%in%c("binomial","gaussian","Gamma","inverse.gaussian",
                  "poisson","quasi","quasibinomial","quasipoisson")){
      Modelo1.snk <- glm(YY~-1+Grupo_snk,data=Data_aux,family=model)}}
  
  if(is.null(Modelo1.scheffe)){
    OUTPUT <- list(
      model=model,
      k.max=k.max,
      aic=data.frame(
        k=2:k.max,
        Method1=aic_aux_met1,
        Method2=aic_aux_met2
      ),
      lr_p=data.frame(
        k=2:k.max,
        Method1=lr_p_v_met1,
        Method2=lr_p_v_met2
      ),
      scheffe=NULL,
      snk=list(
        k=length(Modelo1.snk$coefficients),
        aic=AIC(Modelo1.snk)/AIC(modelo1),
        lr_p=lrtest(modelo1,Modelo1.snk)[2,5],
        Groups=snk_collapsed
      )
    )
  } else{
    OUTPUT <- list(
      model=model,
      k.max=k.max,
      aic=data.frame(
        k=2:k.max,
        Method1=aic_aux_met1,
        Method2=aic_aux_met2
      ),
      lr_p=data.frame(
        k=2:k.max,
        Method1=lr_p_v_met1,
        Method2=lr_p_v_met2
      ),
      scheffe=list(
        k=length(Modelo1.scheffe$coefficients),
        aic=AIC(Modelo1.scheffe)/AIC(modelo1),
        lr_p=lrtest(modelo1,Modelo1.scheffe)[2,5],
        Groups=scheffe_collapsed
      ),
      snk=list(
        k=length(Modelo1.snk$coefficients),
        aic=AIC(Modelo1.snk)/AIC(modelo1),
        lr_p=lrtest(modelo1,Modelo1.snk)[2,5],
        Groups=snk_collapsed
      )
    )
  }
  
  
  
  
  return(OUTPUT)
  options(scipen=0)
}