#' ridge cox model using glmnet
#'
#'
#' @param r a numeric value, a seed to run this method
#' @param data a dataframe, the data used to performance this survival model
#' @param cvK a numeric value, cross-validation fold
#' @param numm a numeric value, the number of variables,i.e.for example, number of proteins in the data
#' @param topnumm a numeric value, the number of variables selected to be passed into the model, for example, the number of DE genes
#' @param fitform_ogl a Surv object from package survival, the survival function
#' @param formula1 a Surv object from package survival, to calculate a version of the brier score, details please check package pec
#' @param formula2 a Surv object from package survival, to calculate a version of the brier score, details please check package pec
#' @param formula3 a Surv object from package survival, to calculate a version of the brier score, details please check package pec
#' @param formula4 a Surv object from package survival, to calculate a version of the brier score, details please check package pec
#' @param timess a numeric vector, contains time points to get the time-dependent AUC values
#' @return a data.frame with allevaluation measurements in all columns and rows are each fold results from cross-validation
#'
#' @examples
#' data("exampledt", package = "SurvBenchmark")
#' fitform_ogl=survival::Surv(time,status)~.
#' formula1=fitform_ogl
#' formula2=fitform_ogl
#' formula3=survival::Surv(time,status)~1
#' formula4=survival::Surv(time,status)~1
#' form1=as.formula(~.)
#' timess=seq(as.numeric(summary(cancerdt2_1$time)[2]),as.numeric(summary(cancerdt2_1$time)[5]),(as.numeric(summary(cancerdt2_1$time)[5])-as.numeric(summary(cancerdt2_1$time)[2]))/14)
#' want=glmnet_lassocox_fun(1,cancerdt2_1[,-dim(cancerdt2_1)[2]],5,formula1,formula2,formula3,formula4,timess);
#' @export

glmnet_ridge_fun=function(r,data,cvK,formula1,formula2,formula3, formula4, timess){
  set.seed(r)
  print(r)
  cvSets = cvTools::cvFolds(nrow(data), cvK)  # permute all the data, into 5 folds
  bicfun=purrr::possibly(function(j){
    test_id = cvSets$subsets[cvSets$which == j]
    test = data[test_id, ]
    train = data[-test_id, ]
    tr_matrix=data.matrix(train[,!names(train)%in% c("time","status")])
    te_matrix=data.matrix(test[,!names(test)%in% c("time","status")])
    tr_st=train$time
    tr_y=train$status
    te_st=test$time
    te_y=test$status

    fit0=glmnet::cv.glmnet(tr_matrix, survival::Surv(tr_st,tr_y), family="cox", standardize = F,alpha=0, nfolds = 5,type.measure = "C")
    fit=glmnet::glmnet(tr_matrix, survival::Surv(tr_st,tr_y), family="cox", standardize = F,alpha=0,lambda = fit0$lambda.min,type.measure = "C")
    lpnew=predict(fit,newx = te_matrix,type="link") #this is the linear predictors (risk?)


    lp<- predict(fit,newx = tr_matrix,type="link")

    #cindex
    harrelC1 <- Hmisc::rcorr.cens(-lpnew,Surv(test$time,test$status))
    hc_1<-harrelC1["C Index"]
    Surv.rsp <- survival::Surv(train$time, train$status)
    Surv.rsp.new <- survival::Surv(test$time, test$status)
    bc_1 <- survAUC::BeggC(Surv.rsp, Surv.rsp.new,lp, lpnew)
    unoc_1<-survAUC::UnoC(Surv.rsp, Surv.rsp.new, lpnew)
    ghc_1<-survAUC::GHCI(lpnew)


    #br

    briers1 <- survAUC::predErr(Surv.rsp, Surv.rsp.new, lp, lpnew,times=test$time, type = "brier", int.type = "unweighted")$error
    br1<-sum(na.omit(briers1))
    briers2<-survAUC::predErr(Surv.rsp, Surv.rsp.new, lp, lpnew,times=test$time, type = "brier", int.type = "weighted")$error
    br2<-sum(na.omit(briers2))
    ibsfun1=purrr::possibly(function(modell){
      briers3 <- pec::pec(list("cox1"=modell),data=test,formula=fitform_ogl1,cens.model="cox")
      return(crps(briers3)[2])
    },otherwise = NA)
    br3<-NA
    ibsfun2=purrr::possibly(function(modell){
      briers4 <- pec::pec(list("cox1"=modell),data=test,formula=fitform_ogl2,cens.model="marginal")
      return(crps(briers4)[2])
    },otherwise = NA)
    br4<-NA
    ibsfun3=purrr::possibly(function(modell){
      briers5 <- pec::pec(list("cox1"=modell),data=test,formula=fitform_ogl3,cens.model="cox")
      return(crps(briers5)[2])
    },otherwise = NA)
    br5<-NA
    ibsfun4=purrr::possibly(function(modell){
      briers6 <- pec::pec(list("cox1"=modell),data=test,formula=fitform_ogl4,cens.model="marginal")
      return(crps(briers6)[2])
    },otherwise = NA)
    br6<-NA
    #time-dependent auc
    times <- timess
    AUC_CD <- survAUC::AUC.uno(Surv.rsp, Surv.rsp.new, lpnew, times)
    a1=AUC_CD$auc[1]
    a2=AUC_CD$auc[2]
    a3=AUC_CD$auc[3]
    a4=AUC_CD$auc[4]
    a5=AUC_CD$auc[5]
    a6=AUC_CD$auc[6]
    a7=AUC_CD$auc[7]
    a8=AUC_CD$auc[8]
    a9=AUC_CD$auc[9]
    a10=AUC_CD$auc[10]
    a11=AUC_CD$auc[11]
    a12=AUC_CD$auc[12]
    a13=AUC_CD$auc[13]
    a14=AUC_CD$auc[14]
    a15=AUC_CD$auc[15]
    a=AUC_CD$iauc
    return(c(hc_1,bc_1,unoc_1,ghc_1,br1,br2,br3,br4,br5,br6,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a))},otherwise=NA)

  cv5_result=rbind.data.frame(bicfun(1),bicfun(2),bicfun(3),bicfun(4),bicfun(5))
  #colnames(cv5_result)=c("hc_1","bc_1","unoc_1","ghc_1","br1","br2","br3","br4","br5","br6","a1","a2","a3","a4","a5","a6","a7","a8","a9","a10","a11","a12","a13","a14","a15","a")


  return(cv5_result)}
