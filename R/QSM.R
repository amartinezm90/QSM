proyectorQSM = function(x,n){
  x = as.matrix(x)
  P = matrix(0,nrow = ncol(x),ncol = ncol(x))
  for (i in 1:n){
    P = P + outer(x[i,],(x[i,]))
  }
  return(P)
}

PositivePoryector = function(proyector,limit=0){
  for (i in 1:ncol(proyector)){
    for (j in 1:nrow(proyector)){
      if (proyector[i,j] < limit){
        proyector[i,j] = 0
      }
    }
  }
  return(proyector)
}

stateVectorEstimationQSM = function(x,A,B){
  abs((t(x)%*%A%*%x)-(t(x)%*%B%*%x))
}

QSM = function(A,B,st){
  return(((norm(B%*%((A%*%st)/norm(A%*%st))))^2)*(norm(A%*%st))^2)
}

InterferenceEffect = function(A,B,st){
  return(norm(A%*%st)^2-norm(B%*%A%*%st)^2-norm(B%*%(diag(1,ncol(A))-A)%*%st)^2)
}

QSFunction = function(BasisA, BasisB, ndimA, ndimB, positive = TRUE, ndimLSA){

  BA = read.csv(BasisA, sep = "", header = F)
  BB = read.csv(BasisB, sep = "", header = F)

  BA_SE = as.matrix(BA[1:ndimA,])
  BB_SE = as.matrix(BB[1:ndimB,])

  library(reshape2)
  library(ggplot2)
  meaningfulVectors = round(BA_SE%*%t(BB_SE),2)
  meaningfulVectors_melted <- melt(meaningfulVectors)
  plotAB = ggplot(data = meaningfulVectors_melted, aes(x=Var1, y=Var2, fill=value)) +
    geom_tile()

  PA = proyectorQSM(BA, ndimA)
  PB = proyectorQSM(BB, ndimB)

  if (positive == TRUE){
    PA = PositivePoryector(PA)
    PB = PositivePoryector(PB)
  }


  stateVectorMean = rep((mean(as.matrix(BA)) + mean(as.matrix(BB)))/2,ndimLSA)
  stVecOptim = optim(par = stateVectorMean, fn = stateVectorEstimationQSM, A=PA, B=PB ,
                     lower = rep(-1,ncol(PA)), upper = rep(1,ncol(PA)), method = "L-BFGS-B")
  stVecOptim = stVecOptim$par
  stateVectorOptimTest = list(zero_optimization = ((t(stVecOptim) %*% PA %*% stVecOptim) -
                                (t(stVecOptim) %*% PB %*% stVecOptim)),state_vector_sum = sum(stVecOptim))

  simAB = QSM(PA,PB,stVecOptim)
  simBA = QSM(PB,PA,stVecOptim)

  interference = InterferenceEffect(PA,PB,stVecOptim)

  simData = data.frame(rbind(simAB, simBA, interference))
  simData <- cbind(rownames(simData), data.frame(simData, row.names=NULL))
  colnames(simData) = c("name","sim")
  simData <- simData[order(simData$sim),]

  simPlot = ggplot(simData, aes(x = reorder(name, -sim), y = sim)) +
    geom_bar(stat = "identity", fill="steelblue")+
    theme_minimal() +
    coord_flip()

  return(list(correlation_matrix = meaningfulVectors, correlation_plot = plotAB, state_vector_test = stateVectorOptimTest,
              AB_similarity = simAB, BA_similarity = simBA, interference_effect = interference, similarity_df = simData,
              similarity_plot = simPlot))
}
