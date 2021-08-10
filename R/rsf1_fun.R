#' random survival forest model
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
#' @return a data.frame with allevaluation measurements in all columns and rows are each fold results from cross-validation
#'
#' @examples
#' data("exampledt2", package = "SurvBenchmark")
#'
#' xnam <- paste(colnames(veteran)[c(1,2,5,6,7,8)], sep="")
#' form=as.formula(paste("Surv(time, status)~ ", paste(xnam, collapse= "+")))
#' fitform_ogl=Surv(time, status) ~ trt + celltype + karno + diagtime + age + prior
#' formula1=Surv(time, status) ~ trt + celltype + karno + diagtime + age + prior
#' formula2=Surv(time, status) ~ trt + celltype + karno + diagtime + age + prior
#' formula3=survival::Surv(time,status)~1
#' formula4=survival::Surv(time,status)~1
#' timess=seq(as.numeric(summary(veteran$time)[2]),as.numeric(summary(veteran$time)[5]),(as.numeric(summary(veteran$time)[5])-as.numeric(summary(veteran$time)[2]))/14)
#' want=rsf1_fun(1,veteran,5,fitform_ogl,formula1,formula2,formula3, formula4, timess);
#' @export


rsf1_fun=function(r,data,cvK,fitform_ogl,formula1,formula2,formula3, formula4, timess){
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
    set.seed(123)
    out.rsf.1 <- randomForestSRC::rfsrc(fitform_ogl,
                       data = train,
                       ntree = 1000,
                       mtry = 10,
                       tree.err=TRUE,
                       importance = TRUE)
    print("success")
    pred.test.fin = randomForestSRC::predict( out.rsf.1,
                             newdata = test,
                             importance = "none" ,type="response",block.size = 1)
    pred.test.fin2 = randomForestSRC::predict( out.rsf.1,
                              newdata = train,
                              importance = "none" ,type="response",block.size = 1)
    pred_te=pred.test.fin$predicted
    pred_tr=pred.test.fin2$predicted

    #predd<-predictSurvProb(original_cox1,newdata=test,times=seq(365,365*15,365))
    #harrel cindex
    harrelC1 <- Hmisc::rcorr.cens(-pred_te,with(test,Surv(time,status)))
    hc<-harrelC1["C Index"]
    print(hc)
    #begg cindex
    lp<- pred_tr
    lpnew <- pred_te
    Surv.rsp <- survival::Surv(train$time, train$status)
    Surv.rsp.new <- survival::Surv(test$time, test$status)
    bc <- survAUC::BeggC(Surv.rsp, Surv.rsp.new,lp, lpnew)
    #uno cindex
    unoc<-survAUC::UnoC(Surv.rsp, Surv.rsp.new, lpnew)
    #gh cindex
    ghc<-survAUC::GHCI(lpnew)
    #br
    briers1 <- survAUC::predErr(Surv.rsp, Surv.rsp.new, lp, lpnew,times=test$time, type = "brier", int.type = "unweighted")$error
    br1<-sum(na.omit(briers1))
    briers2<-survAUC::predErr(Surv.rsp, Surv.rsp.new, lp, lpnew,times=test$time, type = "brier", int.type = "weighted")$error
    br2<-sum(na.omit(briers2))
    ibsfun1=purrr::possibly(function(modell){
      briers3 <- pec::pec(list("cox1"=modell),data=test,formula=formula1,cens.model="cox")
      return(crps(briers3)[2])
    },otherwise = NA)
    #briers3 <- pec(list("cox1"=original_cox1),data=test,formula=Surv(tx_gperiod,tx_gstatus)~recip_sex+recip_eth+recip_age+recip_height+recip_weight+recip_smoker+recip_lung+recip_coronary+recip_pvd+recip_cvd+recip_diabetes+recip_waittime+donor_age+donor_sex+donor_height+donor_weight+donor_causedeath_cva+donor_dcd+donor_diabetes+donor_ht+donor_smoker+donor_creatinine+tx_ischaemia+tx_misa+tx_misb+tx_misdr,cens.model="cox")
    #bs3[j]<-crps(briers3)[2]
    br3<-ibsfun1(out.rsf.1)
    ibsfun2=purrr::possibly(function(modell){
      briers4 <- pec::pec(list("cox1"=modell),data=test,formula=formula2,cens.model="marginal")
      return(crps(briers4)[2])
    },otherwise = NA)
    br4<-ibsfun2(out.rsf.1)
    ibsfun3=purrr::possibly(function(modell){
      briers5 <- pec::pec(list("cox1"=modell),data=test,formula=formula3,cens.model="cox")
      return(crps(briers5)[2])
    },otherwise = NA)
    br5<-ibsfun3(out.rsf.1)
    ibsfun4=purrr::possibly(function(modell){
      briers6 <- pec::pec(list("cox1"=modell),data=test,formula=formula4,cens.model="marginal")
      return(crps(briers6)[2])
    },otherwise = NA)
    br6<-ibsfun4(out.rsf.1)
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
    return(c(hc,bc,unoc,ghc,br1,br2,br3,br4,br5,br6,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a))},otherwise=NA)
  print(bicfun(1))
  cv5_result=rbind.data.frame(bicfun(1),bicfun(2),bicfun(3),bicfun(4),bicfun(5))
  print(class(cv5_result))
  #colnames(cv5_result)=c("hc","bc","unoc","ghc","br1","br2","br3","br4","br5","br6","a1","a2","a3","a4","a5","a6","a7","a8","a9","a10","a11","a12","a13","a14","a15","a")

  # want=cbind.data.frame(hc_acc5,bc_acc5,unoc_acc5,ghc_acc5,bs1,bs2,bs3,bs4,bs5,bs6,auc1,auc2,auc3,auc4,auc5,auc6,auc7,auc8,auc9,auc10,auc11,auc12,auc13,auc14,auc15,auc)
  return(cv5_result)}
