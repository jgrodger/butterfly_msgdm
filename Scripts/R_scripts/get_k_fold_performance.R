#' k-fold validation of Zeta.msgdm for restricted argument values
#' 
#' `get_k_fold_performance` returns a dataframe of performance metrics for
#' msgdm split over k folds
#'
#' This validates msgdms by comparing predicted to observed zeta values.
#' 
#' The code is based on the example for Predict.msgdm
#' 
#'  The procedure for each fold is
#' (1) run an msgdm on training data
#' (2) Get i-spline-transformed predictor data for the training set that is 
#' needed to predict zeta diversity values for the test data, using the model 
#' from the training data, and predict zeta diversity (Predict.msgdm())
#' (3) Run an msgdm on the test data
#' (4) compare predicted zeta diversity values from (2) with observed 
#' values from (3) to calculate model performance metrics, and do same for a 
#' null model with the mean observed value of zeta as the predicted value.
#' 
#' POTENTIAL IMPROVEMENTS/QUERIES
#' Could I get away with a much smaller sam for the predicted data?
#' 
#' RMSE (Root mean square error) is the SD of residuals
#' MAE (Mean absolute error) is the mean of absolute value of residuals. 
#' MAE is less affected by outliers than RMSE because of lack of squaring
#' R2 = r-squared 
#' 
#' The first seven arguments as for Zeta.msgdm.
#' 
#' The following, and defaults not changed, are hardwired below for the time being
#' 
#' family=binomial(link="log"),
#' cons.inter = -1, 
#' glm.init = TRUE)
#' 
#' seeds can be changed so that performance can be assessed across multiple runs of the function with different random samples
#' 
#' @param k number of folds 
#' @param seed_shuffle: seed to shuffle data
#' 
#' @returns A dataframe with RMSE, R2(r-squared), and MSE as columns and k rows 


get_k_fold_performance <- function (
  data.spec,
  data.env,
  xy = NULL,
  order = 1,
  sam = 1000,
  reg.type = "ispline",
  normalize = "Simpson",
  k = 5, 
  seed_shuffle = 1) {

# Divide the data into folds and create a dataframe to store model and null model performance values
# of all the folds.

folds <- cut(seq(1, nrow(data.spec)), breaks = k, labels = FALSE)
set.seed(seed_shuffle)
folds <-folds[sample(nrow(data.spec))]
  
performance <- data.frame((matrix(NA, k, 6)))
names(performance) <- c("train_R2", "model_RMSE", "model_R2", "model_MAE",
                              "null_RMSE", "null_MAE")

#Perform k-fold cross-validation
for(i in 1:k){
  testIndexes <- which(folds == i, arr.ind=TRUE)
  
  train.site.by.species <- data.spec[-testIndexes, ]
  test.site.by.species <- data.spec[testIndexes, ]
  
  train.site.by.env <- data.env[-testIndexes, ]
  test.site.by.env <- data.env[testIndexes, ]   
  
  train.site.by.xy <- xy[-testIndexes, ]
  test.site.by.xy <- xy[testIndexes, ]
  
  # (1) Run msgdm on the training data
  set.seed(i)
  train.msgdm <- Zeta.msgdm(train.site.by.species, train.site.by.env, 
    xy = train.site.by.xy, sam = sam, order = order, reg.type = reg.type, 
    normalize = normalize, family=binomial(link="log"), 
    cons.inter = -1, glm.init = TRUE)
  
  # (2) Get i-spline-transformed predictor data for the training set 
  # and predict zeta diversity (Predict.msgdm())
  # Here we apply ispline transformation to the environmental data, 
  # create a dataframe to store mean distances between ispline transformed variables
  # create a vector to store mean pairwise geographical distance between sites
  # For geographical distance, ispline transformation is calculated after...
  
  data.splines <- Ispline(test.site.by.env) 
  test.data <- data.frame(matrix(NA, sam, ncol(data.splines$splines)))
  names(test.data) <- names(data.splines$splines)
  distance <- rep(NA, sam)
  
  
  # For test.data, the inner apply() first gets mean distance between sites in 
  # values of each spline from environmental data, the outer apply() gets 
  # mean means over all sites in samp 
  # For distance, apply gets the average pairwise euclidean geographical 
  # distance between sites
  
  set.seed(i)
  for(z in 1:sam){
    samp <- sample(1:nrow(test.site.by.species), order, replace = FALSE)
    test.data[z,] <- apply(apply(data.splines$splines[samp,], 2, stats::dist), 2, mean) 
    distance[z] <- apply(as.matrix(c(stats::dist(test.site.by.xy[samp, ]))), 2, mean)
   }
  
  #Ispline transform distance values, add to test data, and rescale test data using 
  # the same rescaling factors used by the training msgdm
  
  distance.splines <- Ispline(data.frame(distance))
  test.data <- cbind(test.data, distance.splines$splines)
  test.data <- test.data/matrix(rep(train.msgdm$rescale.factor,sam),  
    sam,length(train.msgdm$rescale.factor), byrow=TRUE)
  
  # Predict zeta diversity from test data using the model from the training data
  predict.msgdm <- Predict.msgdm(model.msgdm = train.msgdm$model, reg.type = reg.type, newdata = test.data)

  # (3) Run an msgdm on the test data
  
  # estimate zeta diversity for samples in the test data
  set.seed(i)
  zeta_mc <- Zeta.order.mc.sam(data.spec = test.site.by.species, normalize = normalize, order=order, sam = sam)
  observed_zeta <- zeta_mc$zeta.val.vec
  
  # (4) compare predicted zeta diversity values from (2) 
  # with observed values from (3) to calculate model performance metrics
  
  
 # plot(observed_zeta ~ predict.msgdm, main = paste("k =", i), ylab = "Observed",
 #  xlab = "Predicted")
  
  
  train_R2 <- with(summary(train.msgdm$model), 1 - deviance/null.deviance)
  
  mean_obs_zeta <- mean(observed_zeta)

  model_RMSE <- sqrt(mean((observed_zeta - predict.msgdm)^2))
  model_R2 <- (cor(predict.msgdm, observed_zeta))^2
  model_MAE <- mean(abs((observed_zeta - predict.msgdm)))
  
  null_RMSE <- sqrt(mean((observed_zeta - mean_obs_zeta)^2))
  null_MAE <- mean(abs((observed_zeta - mean_obs_zeta)))
  
  performance[i,] <- c(train_R2, model_RMSE, model_R2, model_MAE, 
                      null_RMSE, null_MAE)
}

return(performance)

}

