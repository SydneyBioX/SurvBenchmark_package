
#' penalised coxph model lasso
#'
#'
#' @param r a numeric value, a seed to run this method
#' @param data a dataframe, the data used to performance this survival model
#' @param cvK a numeric value, cross-validation fold
#' @param fitform_ogl a Surv object from package survival, the survival function
#' @param formula1 a Surv object from package survival, to calculate a version of the brier score, details please check package pec
#' @param formula2 a Surv object from package survival, to calculate a version of the brier score, details please check package pec
#' @param formula3 a Surv object from package survival, to calculate a version of the brier score, details please check package pec
#' @param formula4 a Surv object from package survival, to calculate a version of the brier score, details please check package pec
#' @param timess a numeric vector of length 15, contains time points to get the time-dependent AUC values
#' @param lambda_1 a numeric value, the tuning parameter used in the function,check package "penalized" for details
#' @param lambda_2 a numeric value, the tuning parameter used in the function,check package "penalized" for details
#' @return a data.frame with allevaluation measurements in all columns and rows are each fold results from cross-validation
#'
#' @examples
#' data("exampledt2", package = "SurvBenchmark")
#'
#' xnam <- paste(colnames(veteran)[c(1,2,5,6,7,8)], sep="")
#' form=as.formula(paste("survival::Surv(time, status)~ ", paste(xnam, collapse= "+")))
#' fitform_ogl=form
#' formula1=fitform_ogl
#' formula2=fitform_ogl
#' formula3=survival::Surv(time,status)~1
#' formula4=survival::Surv(time,status)~1
#' timess=seq(as.numeric(summary(veteran$time)[2]),as.numeric(summary(veteran$time)[5]),(as.numeric(summary(veteran$time)[5])-as.numeric(summary(veteran$time)[2]))/14)
#' want=p_cox1_fun(1,veteran,5,fitform_ogl,formula1,formula2,formula3, formula4, timess,1,0);
#' @export


p_cox1_fun=function(r,data,cvK,form1,formula1,formula2,formula3, formula4, timess,lambda_1,lambda_2){
  if (! is.numeric(r)) stop("Input seed is wrong")
  if (! is.numeric(cvK)) stop("Input cross-validation fold number is wrong")
  if (is.null(dim(data))) stop("Input data is wrong")
  if (length(timess)!=15) stop("Wrong time vector length")
  if (class(timess)!= "numeric") stop("Wrong time vector type")
  set.seed(r)
  print(r)
  cvSets = cvTools::cvFolds(nrow(data), cvK)  # permute all the data, into 5 folds
  bicfun=purrr::possibly(function(j){
    test_id = cvSets$subsets[cvSets$which == j]
    test = data[test_id, ]
    train = data[-test_id, ]

    fit=penalized::penalized(Surv(time,status),penalized = form1,lambda1 = lambda_1,lambda2 = lambda_2,data = train)
    predd=predict(fit,penalized=form1,data=test)
    predd1=as.matrix(predd)

    lpp<- predict(fit,penalized=form1,data=train)
    lpp1=as.matrix(lpp)

    #cindex
    predd2=predd1[,round(dim(predd1)[2]/2,digits = 0)]
    harrelC1 <- Himsc::rcorr.cens(predd2,Surv(test$time,test$status))
    hc_1<-harrelC1["C Index"]
    lpp2=-lpp1[,round(dim(lpp1)[2]/2,digits = 0)]
    lpnew <- -predd2
    Surv.rsp <- survival::Surv(train$time, train$status)
    Surv.rsp.new <- survival::Surv(test$time, test$status)
    bc_1 <- survAUC::BeggC(Surv.rsp, Surv.rsp.new,lpp2, lpnew)
    unoc_1<-survAUC::UnoC(Surv.rsp, Surv.rsp.new, lpnew)
    ghc_1<-survAUC::GHCI(lpnew)


    #br
    lp=lpp2
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
