
#' dnnsurv model
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
#' @param timess a numeric vector, contains time points to get the time-dependent AUC values
#' @param pickedtime a numeric vector, the picked pesudo time points, see package "DNNSurv" for details
#' @return a data.frame with all evaluation measurements in all columns and rows are each fold results from cross-validation
#' @examples
#' data("exampledt2", package = "SurvBenchmark")
#' xnam <- paste(colnames(veteran)[c(1,2,5,6,7,8)], sep="")
#' form=as.formula(paste("survival::Surv(time, status)~ ", paste(xnam, collapse= "+")))
#' fitform_ogl=form
#' formula1=fitform_ogl
#' formula2=fitform_ogl
#' formula3=survival::Surv(time,status)~1
#' formula4=survival::Surv(time,status)~1
#' pickedtime=c(2,3,5,7,9)
#' timess=seq(as.numeric(summary(veteran$time)[2]),as.numeric(summary(veteran$time)[5]),(as.numeric(summary(veteran$time)[5])-as.numeric(summary(veteran$time)[2]))/14)
#' want=dnnsurv_fun(1,veteran,5,fitform_ogl,formula1,formula2,formula3,formula4,timess,pickedtime);
#' @export


#-----------------------------------------------------------------------------------------------------------------------------------
# get pseudo conditional survival probabilities
# t is the survival time
# d is the censoring indicator
# qt is a vector of time points that are used to divide the time interval
# output has subject id (id) and time points (s) and pseudo conditional survival probabilities (pseudost) for subject=id and at time s
#------------------------------------------------------------------------------------------------------------------------------------
getPseudoConditional <- function(t, d, qt){
  #browser()
  s <- c(0, qt)
  n=length(t)
  ns=length(s)-1  # the number of intervals
  D <- do.call(cbind, lapply(1:ns, function(j)  (s[j] < t) * (t <= s[j+1]) * (d == 1)))
  R <- do.call(cbind, lapply(1:ns, function(j) ifelse(s[j] < t, 1, 0)))
  Delta<-do.call(cbind, lapply(1:ns, function(j) pmin(t,s[j+1])-s[j]))

  # format into long formate
  dd.tmp=cbind.data.frame(id=rep(1:n,ns),s=rep(c(0,qt[-length(qt)]), each=n), y=c(R*Delta),d=c(D))

  dd=dd.tmp[dd.tmp$y>0,]
  pseudost=rep(NA, nrow(dd))
  for (j in 1:ns){
    index= (dd$s==s[j])
    dds=dd[index,]
    pseudost[index]=pseudo::pseudosurv(time=dds$y, event=dds$d, tmax=s[j+1]-s[j])$pseudo
    print(j)
  }
  dd$pseudost=pseudost

  return(dd[,c(1,2,5)])
}


#--------------------------------------------------------------------------------------
# There are two example neural network models used in the paper.
# You should tunning the hyperparamters to find the best neural network for your own study.
#-------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------
# Example of a neural network with two hidden layers implemented in R keras in the paper
#-----------------------------------------------------------------------------------------
pseudoDNN.train <- function(x_train, y_train){

  model <- keras::keras_model_sequential() %>%
    layer_dense(units = 8, kernel_regularizer = regularizer_l2(0.0001), activation = "tanh",
                input_shape = dim(x_train)[[2]]) %>%
    layer_dense(units = 4, kernel_regularizer = regularizer_l2(0.01),
                activation = "tanh") %>%
    layer_dense(units = 1, activation='sigmoid')

  model %>% compile(
    optimizer = optimizer_adam(lr = 0.0025),
    loss = "mse",
    metrics = c("mae")
  )
  model %>% fit(x_train, y_train,
                epochs = 30, batch_size = 64,
                verbose = 0)

  model
}

#----------------------------------------------------------------------------------------
# Another example of a neural network with one hidden layer implemented in R keras in the paper
#-----------------------------------------------------------------------------------------
pseudoDNN.train <- function(x_train, y_train){
  # use selu instead of relu for some studies
  model <- keras::keras_model_sequential() %>%
    layer_dense(units=16,  activation = "selu",bias_initializer = initializer_constant(0.0),
                input_shape = dim(x_train)[[2]]) %>%
    layer_dropout(rate = 0.2) %>%
    layer_dense(units = 1, activation='sigmoid')

  model %>% compile(
    optimizer = optimizer_rmsprop(lr = 0.001),
    loss = "mse",
    metrics = c("mae")
  )
  model %>% fit(x_train, y_train,
                epochs = 1000, batch_size =256,
                verbose = 0)
  model
}

#----------------------------------------------------------------------------
#prediction based on a keras model
#----------------------------------------------------------------------------
pseudoDNN.predict <- function(model, x_test){
  ypred <- model %>% predict(x_test)
  ypred
}


dnnsurv_fun=function(r,data,cvK,fitform_ogl,formula1,formula2,formula3, formula4, timess,pickedtime){
  set.seed(r)
  print(r)
  cvSets = cvTools::cvFolds(nrow(data), cvK)  # permute all the data, into 5 folds
  bicfun=purrr::possibly(function(j){
    test_id = cvSets$subsets[cvSets$which == j]
    test = data[test_id, ]
    train = data[-test_id, ]
    set.seed(123)
    x_train=as.matrix(subset(train, select=-c(status,time)))
    y_train = cbind(time = train$time, status = train$status)

    x_test=as.matrix(subset(test,select=-c(status,time)))
    y_test = cbind(time = test$time, status = test$status)


    #pickTime=c(0.7,1.7,3.2,5.3,8.3)
    pickTime=pickedtime
    x_train  <- x_train
    surv_train <- y_train[,1]
    cen_train <- y_train[,2]

    # data normalization
    mean <- apply(as.matrix(x_train), 2, mean)
    std <- apply(as.matrix(x_train), 2, sd)
    x_train <- scale(x_train, center = mean, scale = std)

    # data normalization
    x_test<- scale(x_test, center = mean, scale = std)
    surv_test <- y_test[,1]
    cen_test =y_test[,2]


    # get the pseudo conditinal survival probability
    pseudoCond  = getPseudoConditional(surv_train, cen_train, pickTime)

    # covaraite
    x <- x_train[pseudoCond$id,]

    # create dummy variables for the time points
    smatrix=model.matrix(~as.factor(pseudoCond$s)+0)

    #create input predictors
    x_train.all <- cbind(x, smatrix)

    y_train.all <- pseudoCond$pseudost

    model = pseudoDNN.train(x_train.all, y_train.all)


    # format the test data
    x_test.all=do.call(rbind, replicate(length(pickTime), x_test, simplify=FALSE))
    s_test=rep(pickTime,each=nrow(x_test))
    smatrix.test=model.matrix(~as.factor(s_test)+0)
    x_test.all=cbind(x_test.all,smatrix.test)

    # predict test data
    ypred.con <- pseudoDNN.predict(model, x_test.all)

    # obtain the marginal survival probability by multiple series of conditional probabilities
    ypred.con <- matrix(ypred.con, nrow=nrow(x_test))
    #ypred <- lapply(1:length(s), function(i) apply(ypred.con[,1:i, drop=FALSE], 1, prod)) #i changed this original one to the one below
    ypred <- lapply(1:ncol(ypred.con), function(i) apply(ypred.con[,1:i, drop=FALSE], 1, prod))
    surv_prob <- Reduce(cbind, ypred)#for those 5 time points

    pred_te =surv_prob

    # predict train data
    ypred.con <- pseudoDNN.predict(model, x_train.all)

    # obtain the marginal survival probability by multiple series of conditional probabilities
    ypred.con <- matrix(ypred.con, nrow=nrow(x_train))
    #ypred <- lapply(1:length(s), function(i) apply(ypred.con[,1:i, drop=FALSE], 1, prod)) #i changed this original one to the one below
    ypred <- lapply(1:ncol(ypred.con), function(i) apply(ypred.con[,1:i, drop=FALSE], 1, prod))
    surv_prob <- Reduce(cbind, ypred)#for those 5 time points

    pred_tr =surv_prob

    #predd<-predictSurvProb(original_cox1,newdata=test,times=seq(365,365*15,365))
    #harrel cindex
    harrelC1 <- Hmisc::rcorr.cens(pred_te[,1],with(test,Surv(time,status)))
    hc<-harrelC1["C Index"]
    #begg cindex
    lp<- -pred_tr[,1]
    lpnew <- -pred_te[,1]
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

  cv5_result=rbind.data.frame(bicfun(1),bicfun(2),bicfun(3),bicfun(4),bicfun(5))
  #print(dim(bicfun(1)))
  #print(dim(bicfun(2)))
  #print(dim(bicfun(3)))
  #print(dim(bicfun(4)))
  #print(dim(bicfun(5)))
  #colnames(cv5_result)=c("hc","bc","unoc","ghc","br1","br2","br3","br4","br5","br6","a1","a2","a3","a4","a5","a6","a7","a8","a9","a10","a11","a12","a13","a14","a15","a")

  # want=cbind.data.frame(hc_acc5,bc_acc5,unoc_acc5,ghc_acc5,bs1,bs2,bs3,bs4,bs5,bs6,auc1,auc2,auc3,auc4,auc5,auc6,auc7,auc8,auc9,auc10,auc11,auc12,auc13,auc14,auc15,auc)
  return(cv5_result)}
