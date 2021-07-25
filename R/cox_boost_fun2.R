
#' CoxBoost model
#'
#'
#' @param r a numeric value, a seed to run this method
#' @param data a dataframe, the data used to performance this survival model
#' @param cvK a numeric value, cross-validation fold
#' @param numm a numeric value, the number of variables,i.e.for example, number of proteins in the data
#' @param topnumm a numeric value, the number of variables selected to be passed into the model, for example, the number of DE genes
#' @param timess a numeric vector, contains time points to get the time-dependent AUC values
#' @param time1 a numeric value, the time point to calculate the risk, see package "CoxBoost"
#' @param stepnumber a numeric value, the number of stpes performed in the model, see package "CoxBoost"
#' @param penaltynumber a numeric value, the penalty number used in the model, see package "CoxBoost"
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

#' time1=timess[3]
#' stepnumber=10
#' penaltynumber=100
#' want=cox_boost_fun2(1,cancerdt2_1,5,16047,1000,time1,timess,stepnumber,penaltynumber);
#' @export

cox_boost_fun2=function(r,data,cvK,numm,topnumm, time1,timess,stepnumber,penaltynumber){
  set.seed(r)
  print(r)
  cvSets = cvTools::cvFolds(nrow(data), cvK)  # permute all the data, into 5 folds
  bicfun=purrr::possibly(function(j){
    test_id = cvSets$subsets[cvSets$which == j]
    test = data[test_id, ]
    train = data[-test_id, ]
    ## Limma analysis:
    protein4_1=train[,c(1:numm)]
    protein4_2=as.matrix(t(protein4_1))
    groupname <-as.factor( as.character(train$os_class))
    design <- model.matrix(~ groupname + 0)
    fit <- limma::lmFit(protein4_2, design)
    cont.matrix <- limma::makeContrasts(contrasts = "groupnamegood-groupnamepoor", levels=design)
    fit2 <- limma::contrasts.fit(fit, cont.matrix)
    fit2 <- limma::eBayes(fit2)
    tT <- limma::topTable(fit2,number = nrow(fit2))
    tT = tT %>% rownames_to_column("GeneName") %>% select(GeneName, logFC, P.Value, adj.P.Val)
    # head(tT)
    # sum(tT$P.Value<0.05)
    selectedname=tT$GeneName[1:topnumm]


    train=train[,colnames(train)%in%c(selectedname,"status","time")]
    test=test[,colnames(test)%in%c(selectedname,"status","time")]


    #fitform_ogl=as.formula(paste("Surv(time, status)~ ", paste(colnames(train)[1:(dim(train)[2]-2)], collapse= "+")))
    fitform_ogl=survival::Surv(time,status)~.
    form1=as.formula(~.)
    formula1=fitform_ogl
    formula2=fitform_ogl
    formula3=survival::Surv(time,status)~1
    formula4=survival::Surv(time,status)~1

    #form1=as.formula(paste("~ ",paste(colnames(train)[1:(dim(train)[2]-2)], collapse= "+")))

    tr_predictormatrix=train[,-which(colnames(train)%in% c("time","status"))]
    tr_predictormatrix=data.matrix(tr_predictormatrix)
    te_predictormatrix=test[,-which(colnames(test)%in% c("time","status"))]
    te_predictormatrix=data.matrix(te_predictormatrix)
    coxboost=CoxBoost::CoxBoost(time=train$time,status=train$status,x=tr_predictormatrix,stepno=stepnumber,penalty=penaltynumber)
    lpnew=predict(coxboost,type="risk",times=time1,newdata=te_predictormatrix)
    #it is the survival probability #the predicted probability of not yet having had the event at the time points given in times
    lpnew=-lpnew #change to risk

    lp<- predict(coxboost,type="risk",times=time1,newdata=tr_predictormatrix)
    lp=-lp
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
