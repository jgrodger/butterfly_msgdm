#' Simpler validation based on my function get_k_fold_validation
#' 
#' `get_validation_different_seeds_and_samps` returns a dataframe of performance metrics for
#' msgdm split over k folds, which is based on the example for Predict.msgdm
#'
#' This validates msgdms in two ways. First it checks for predictive skill
#' by comparing predicted to observed zeta values and second it checks for 
#' consistency of statistial inference. This method doesn't
#' bother to create folds as the number of possible combinations is vast, so it is 
#' unlikely that any samp will be repeated in training and prediction.Training 
#' and test data here are different samps (randomly selected sets of sites, as
#' per zeta order)
#' 
#' 
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
#' @param seeds: number of different seeds to test
#' @param sam.train: number of samples (samps) for training the model
#' @param sam.test: number of samples for testing the model
#' @returns A list, containing (1) a dataframe with RMSE, R2(r-squared), 
#' and MSE as columns (2) model outputs from MSGDMS 


get_validation_different_seeds_and_samps <- function (
    data.spec,
    data.env,
    xy = NULL,
    order = 1,
    sam.train = 1000,
    sam.test = 1000,    
    reg.type = "ispline",
    normalize = "Simpson",
    seeds = 1) {
  
  # Create a dataframe to store model and null model performance values, 
  # and a list for all output

  
  performance <- data.frame((matrix(NA, seeds, 6)))
  names(performance) <- c("train_R2", "model_RMSE", "model_R2", "model_MAE",
                          "null_RMSE", "null_MAE")
  msgdm_compare <- list()
  msgdm_compare$models <- list()
  
  #Perform validation
  for(i in 1:seeds){

    # (1) Run msgdm on the training data
    set.seed(i)
    train.msgdm <- Zeta.msgdm(site.by.species, site.by.env, 
                    xy = site.by.xy, sam = sam.train, order = order, reg.type = reg.type, 
                    normalize = normalize, family=binomial(link="log"), 
                    cons.inter = -1, glm.init = TRUE)
    
    msgdm_compare$models[[i]] <- summary(train.msgdm$model)
    
    # (2) Get i-spline-transformed predictor data for the training set 
    # and predict zeta diversity (Predict.msgdm())
    # Here we apply ispline transformation to the environmental data, 
    # create a dataframe to store mean distances between ispline transformed variables
    # create a vector to store mean pairwise geographical distance between sites
    # For geographical distance, ispline transformation is calculated after...
    
    data.splines <- Ispline(site.by.env) 
    test.data <- data.frame(matrix(NA, sam.test, ncol(data.splines$splines)))
    names(test.data) <- names(data.splines$splines)
    distance <- rep(NA, sam.test)
    
    
    # For test.data, the inner apply() first gets mean distance between sites in 
    # values of each spline from environmental data, the outer apply() gets 
    # mean means over all sites in samp 
    # For distance, apply gets the average pairwise euclidean geographical 
    # distance between sites
    
    # Set a different seed, which must be the same below for observed zeta sampling
    set.seed(i + 1)
      for(z in 1:sam.test){
      samp <- sample(1:nrow(site.by.species), order, replace = FALSE)
      test.data[z,] <- apply(apply(data.splines$splines[samp,], 2, stats::dist), 2, mean) 
      distance[z] <- apply(as.matrix(c(stats::dist(site.by.xy[samp, ]))), 2, mean)
    }
    
    #Ispline transform distance values, add to test data, and rescale test data using 
    # the same rescaling factors used by the training msgdm
    
    distance.splines <- Ispline(data.frame(distance))
    test.data <- cbind(test.data, distance.splines$splines)
    test.data <- test.data/matrix(rep(train.msgdm$rescale.factor,sam.test),  
                                  sam.test,length(train.msgdm$rescale.factor), byrow=TRUE)
    
    # Predict zeta diversity from test data using the model from the training data
    predict.msgdm <- Predict.msgdm(model.msgdm = train.msgdm$model, reg.type = reg.type, newdata = test.data)
    
    # (3) Run an msgdm on the test data
    
    # estimate zeta diversity for samples in the test data
    set.seed(i + 1)
    zeta_mc <- Zeta.order.mc.sam(data.spec = site.by.species, normalize = normalize, order=order, sam = sam.test)
    observed_zeta <- zeta_mc$zeta.val.vec
    
    # (4) compare predicted zeta diversity values from (2) 
    # with observed values from (3) to calculate model performance metrics
    
    
    # plot(observed_zeta ~ predict.msgdm, main = paste("seed =", i), ylab = "Observed",
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
  
  msgdm_compare$performance <- performance
  
  return(msgdm_compare)
  
}

