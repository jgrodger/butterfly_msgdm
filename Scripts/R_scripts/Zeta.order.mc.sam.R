#' `Zeta.order.mc.sam` is an exact copy of Zeta.order.mc from Zetadiv 1.3.0 
#' except that the vector of sample zeta values is added to the output

Zeta.order.mc.sam <- function (data.spec, xy = NULL, order = 1, sam = 1000, sd.correct = TRUE, 
          sd.correct.adapt = FALSE, rescale = FALSE, normalize = FALSE, 
          NON = FALSE, FPO = NULL, DIR = FALSE, empty.row = "empty", 
          silent = TRUE) 
{
  if (!inherits(data.spec, "data.frame")) {
    stop("Error: ", paste(deparse(substitute(data.spec)), 
                          " is a ", class(data.spec), ". It must be a data frame.", 
                          sep = ""))
  }
  if (order > dim(data.spec)[1]) {
    stop("Error: wrong value for \"order\": it must be equal or lower than the number of sites.")
  }
  if ((NON == TRUE || !is.null(FPO)) & is.null(xy)) {
    stop("Error: if NON = TRUE or !is.null(FPO), xy must be non null.")
  }
  if ((NON == TRUE || !is.null(FPO)) && nrow(data.spec) != 
      nrow(xy)) {
    stop("Error: data.spec and xy must have the same number of rows.")
  }
  if (empty.row == "remove") {
    if (length(which(rowSums(data.spec))) > 0) {
      data.spec <- data.spec[-which(rowSums(data.spec) == 
                                      0), ]
    }
  }
  x <- dim(data.spec)[1]
  if (is.null(FPO)) {
    if (order == 1) {
      zeta.val <- mean(rowSums(data.spec))
      if (sd.correct == TRUE & sd.correct.adapt == FALSE) {
        zeta.val.sd <- stats::sd(rowSums(data.spec))
      }
      else {
        zeta.val.sd <- stats::sd(rowSums(data.spec)) * 
          nrow(data.spec)
      }
      if (rescale == TRUE || normalize != FALSE) {
        zeta.val <- 1
        zeta.val.sd <- zeta.val.sd/mean(rowSums(data.spec))
      }
    }
    else {
      if (NON == FALSE) {
        if (choose(x, order) > sam) {
          if (silent == FALSE) {
            print(paste("Monte Carlo sampling for order", 
                        order))
          }
          u <- rep(NA, sam)
          for (z in 1:sam) {
            samp <- sample(1:x, order, replace = FALSE)
            u[z] <- sum(apply(data.spec[samp, ], 2, prod))
            if (normalize == "Jaccard") {
              toto <- (ncol(data.spec) - sum(apply((1 - 
                                                      data.spec[samp, ]), 2, prod)))
              if (toto == 0) {
                if (empty.row == 0) {
                  u[z] <- 0
                }
                else if (empty.row == 1) {
                  u[z] <- 1
                }
              }
              else u[z] <- u[z]/toto
            }
            else if (normalize == "Sorensen") {
              toto <- (mean(apply(data.spec[samp, ], 
                                  1, sum)))
              if (toto == 0) {
                if (empty.row == 0) {
                  u[z] <- 0
                }
                else if (empty.row == 1) {
                  u[z] <- 1
                }
              }
              else u[z] <- u[z]/toto
            }
            else if (normalize == "Simpson") {
              toto <- (min(apply(data.spec[samp, ], 1, 
                                 sum)))
              if (toto == 0) {
                if (empty.row == 0) {
                  u[z] <- 0
                }
                else if (empty.row == 1) {
                  u[z] <- 1
                }
              }
              else u[z] <- u[z]/toto
            }
          }
        }
        else {
          if (silent == FALSE) {
            print(paste("Exact solution for order", order))
          }
          u <- rep(NA, choose(x, order))
          samp <- utils::combn(1:x, order)
          for (z in 1:dim(samp)[2]) {
            u[z] <- sum(apply(data.spec[samp[, z], ], 
                              2, prod))
            if (normalize == "Jaccard") {
              toto <- (ncol(data.spec) - sum(apply((1 - 
                                                      data.spec[samp[, z], ]), 2, prod)))
              if (toto == 0) {
                if (empty.row == 0) {
                  u[z] <- 0
                }
                else if (empty.row == 1) {
                  u[z] <- 1
                }
              }
              else u[z] <- u[z]/toto
            }
            else if (normalize == "Sorensen") {
              toto <- (mean(apply(data.spec[samp[, z], 
              ], 1, sum)))
              if (toto == 0) {
                if (empty.row == 0) {
                  u[z] <- 0
                }
                else if (empty.row == 1) {
                  u[z] <- 1
                }
              }
              else u[z] <- u[z]/toto
            }
            else if (normalize == "Simpson") {
              toto <- (min(apply(data.spec[samp[, z], 
              ], 1, sum)))
              if (toto == 0) {
                if (empty.row == 0) {
                  u[z] <- 0
                }
                else if (empty.row == 1) {
                  u[z] <- 1
                }
              }
              else u[z] <- u[z]/toto
            }
          }
        }
      }
      else {
        if (x > sam) {
          u <- rep(NA, sam)
          samps <- sample(1:x, sam, replace = FALSE)
          for (z in 1:sam) {
            samp <- samps[z]
            xy.dist <- (xy[, 1] - xy[samp, 1])^2 + (xy[, 
                                                       2] - xy[samp, 2])^2
            samp <- c(samp, order(xy.dist)[2:order])
            u[z] <- sum(apply(data.spec[samp, ], 2, prod))
            if (normalize == "Jaccard") {
              toto <- (ncol(data.spec) - sum(apply((1 - 
                                                      data.spec[samp, ]), 2, prod)))
              if (toto == 0) {
                if (empty.row == 0) {
                  u[z] <- 0
                }
                else if (empty.row == 1) {
                  u[z] <- 1
                }
              }
              else u[z] <- u[z]/toto
            }
            else if (normalize == "Sorensen") {
              toto <- (mean(apply(data.spec[samp, ], 
                                  1, sum)))
              if (toto == 0) {
                if (empty.row == 0) {
                  u[z] <- 0
                }
                else if (empty.row == 1) {
                  u[z] <- 1
                }
              }
              else u[z] <- u[z]/toto
            }
            else if (normalize == "Simpson") {
              toto <- (min(apply(data.spec[samp, ], 1, 
                                 sum)))
              if (toto == 0) {
                if (empty.row == 0) {
                  u[z] <- 0
                }
                else if (empty.row == 1) {
                  u[z] <- 1
                }
              }
              else u[z] <- u[z]/toto
            }
          }
        }
        else {
          u <- rep(NA, x)
          samps <- 1:x
          for (z in 1:x) {
            samp <- samps[z]
            xy.dist <- (xy[, 1] - xy[samp, 1])^2 + (xy[, 
                                                       2] - xy[samp, 2])^2
            samp <- c(samp, order(xy.dist)[2:order])
            u[z] <- sum(apply(data.spec[samp, ], 2, prod))
            if (normalize == "Jaccard") {
              toto <- (ncol(data.spec) - sum(apply((1 - 
                                                      data.spec[samp, ]), 2, prod)))
              if (toto == 0) {
                if (empty.row == 0) {
                  u[z] <- 0
                }
                else if (empty.row == 1) {
                  u[z] <- 1
                }
              }
              else u[z] <- u[z]/toto
            }
            else if (normalize == "Sorensen") {
              toto <- (mean(apply(data.spec[samp, ], 
                                  1, sum)))
              if (toto == 0) {
                if (empty.row == 0) {
                  u[z] <- 0
                }
                else if (empty.row == 1) {
                  u[z] <- 1
                }
              }
              else u[z] <- u[z]/toto
            }
            else if (normalize == "Simpson") {
              toto <- (min(apply(data.spec[samp, ], 1, 
                                 sum)))
              if (toto == 0) {
                if (empty.row == 0) {
                  u[z] <- 0
                }
                else if (empty.row == 1) {
                  u[z] <- 1
                }
              }
              else u[z] <- u[z]/toto
            }
          }
        }
      }
      if (rescale == TRUE & normalize == FALSE) {
        u <- u/ncol(data.spec)
      }
      zeta.val <- mean(u)
      if (sd.correct.adapt == FALSE) {
        if (sd.correct == TRUE) {
          zeta.val.sd <- stats::sd(u)
        }
        else {
          zeta.val.sd <- stats::sd(u) * sqrt((length(u) - 
                                                1)/length(u))
        }
      }
      else {
        if (x > sam) {
          zeta.val.sd <- stats::sd(u) * sqrt((length(u) - 
                                                1)/length(u))
        }
        else {
          zeta.val.sd <- stats::sd(u)
        }
      }
    }
  }
  else {
    if (DIR == FALSE) {
      xy.dist <- (FPO[1] - xy[, 1])^2 + (FPO[2] - xy[, 
                                                     2])^2
      samp <- order(xy.dist)[1:order]
      u <- sum(apply(data.spec[samp, ], 2, prod))
      if (normalize == "Jaccard") {
        toto <- (ncol(data.spec) - sum(apply((1 - data.spec[samp, 
        ]), 2, prod)))
        if (toto == 0) {
          if (empty.row == 0) {
            u <- 0
          }
          else if (empty.row == 1) {
            u <- 1
          }
        }
        else u <- u/toto
      }
      else if (normalize == "Sorensen") {
        toto <- (mean(apply(data.spec[samp, ], 1, sum)))
        if (toto == 0) {
          if (empty.row == 0) {
            u <- 0
          }
          else if (empty.row == 1) {
            u <- 1
          }
        }
        else u <- u/toto
      }
      else if (normalize == "Simpson") {
        toto <- (min(apply(data.spec[samp, ], 1, sum)))
        if (toto == 0) {
          if (empty.row == 0) {
            u <- 0
          }
          else if (empty.row == 1) {
            u <- 1
          }
        }
        else u <- u/toto
      }
    }
    else {
      if (order == 1) {
        zeta.val <- mean(rowSums(data.spec))
        if (sd.correct == TRUE & sd.correct.adapt == 
            FALSE) {
          zeta.val.sd <- stats::sd(rowSums(data.spec))
        }
        else {
          zeta.val.sd <- stats::sd(rowSums(data.spec)) * 
            nrow(data.spec)
        }
        if (rescale == TRUE || normalize != FALSE) {
          zeta.val <- 1
          zeta.val.sd <- zeta.val.sd/mean(rowSums(data.spec))
        }
      }
      else {
        xy.FPO <- as.matrix(xy - FPO)
        if (x > sam) {
          u <- rep(NA, sam)
          samps <- sample(1:x, sam, replace = FALSE)
          for (z in 1:sam) {
            samp <- samps[z]
            xy0 <- xy.FPO[samp, ]
            no <- sqrt(sum(xy0^2))
            R <- matrix(c(xy0[2]/no, xy0[1]/no, -xy0[1]/no, 
                          xy0[2]/no), 2, 2)
            xy.FPO.tr <- apply(xy.FPO, 1, function(xy.FPO, 
                                                   R, xy0, no) {
              R %*% matrix(xy.FPO, 2, 1) - c(0, no)
            }, R, xy0, no)
            xy.FPO.tr <- as.data.frame(t(xy.FPO.tr))
            xy.FPO.tr[samp, ] <- c(0, 0)
            xy.FPO.tr[which(xy.FPO.tr[, 2] <= 0), ] <- NA
            xy.dist <- xy.FPO.tr[, 1]^2 + xy.FPO.tr[, 
                                                    2]^2
            if (length(which(!is.na(xy.dist))) >= (order - 
                                                   1)) {
              samp <- c(samp, order(xy.dist)[1:(order - 
                                                  1)])
              u[z] <- sum(apply(data.spec[samp, ], 2, 
                                prod))
              if (normalize == "Jaccard") {
                toto <- (ncol(data.spec) - sum(apply((1 - 
                                                        data.spec[samp, ]), 2, prod)))
                if (toto == 0) {
                  if (empty.row == 0) {
                    u[z] <- 0
                  }
                  else if (empty.row == 1) {
                    u[z] <- 1
                  }
                }
                else u[z] <- u[z]/toto
              }
              else if (normalize == "Sorensen") {
                toto <- (mean(apply(data.spec[samp, ], 
                                    1, sum)))
                if (toto == 0) {
                  if (empty.row == 0) {
                    u[z] <- 0
                  }
                  else if (empty.row == 1) {
                    u[z] <- 1
                  }
                }
                else u[z] <- u[z]/toto
              }
              else if (normalize == "Simpson") {
                toto <- (min(apply(data.spec[samp, ], 
                                   1, sum)))
                if (toto == 0) {
                  if (empty.row == 0) {
                    u[z] <- 0
                  }
                  else if (empty.row == 1) {
                    u[z] <- 1
                  }
                }
                else u[z] <- u[z]/toto
              }
            }
            else {
              if (silent == FALSE) {
                print(paste("warning: number of sites away from the FPO too low to compute zeta for order", 
                            order))
              }
              u[z] <- NA
            }
          }
        }
        else {
          u <- rep(NA, x)
          samps <- 1:x
          for (z in 1:(x - order + 1)) {
            samp <- samps[z]
            xy0 <- xy.FPO[samp, ]
            no <- sqrt(sum(xy0^2))
            R <- matrix(c(xy0[2]/no, xy0[1]/no, -xy0[1]/no, 
                          xy0[2]/no), 2, 2)
            xy.FPO.tr <- apply(xy.FPO, 1, function(xy.FPO, 
                                                   R, xy0, no) {
              R %*% matrix(xy.FPO, 2, 1) - c(0, no)
            }, R, xy0, no)
            xy.FPO.tr <- as.data.frame(t(xy.FPO.tr))
            xy.FPO.tr[samp, ] <- c(0, 0)
            xy.FPO.tr[which(xy.FPO.tr[, 2] <= 0), ] <- NA
            xy.dist <- xy.FPO.tr[, 1]^2 + xy.FPO.tr[, 
                                                    2]^2
            if (length(which(!is.na(xy.dist))) >= (order - 
                                                   1)) {
              samp <- c(samp, order(xy.dist)[1:(order - 
                                                  1)])
              u[z] <- sum(apply(data.spec[samp, ], 2, 
                                prod))
              if (normalize == "Jaccard") {
                toto <- (ncol(data.spec) - sum(apply((1 - 
                                                        data.spec[samp, ]), 2, prod)))
                if (toto == 0) {
                  if (empty.row == 0) {
                    u[z] <- 0
                  }
                  else if (empty.row == 1) {
                    u[z] <- 1
                  }
                }
                else u[z] <- u[z]/toto
              }
              else if (normalize == "Sorensen") {
                toto <- (mean(apply(data.spec[samp, ], 
                                    1, sum)))
                if (toto == 0) {
                  if (empty.row == 0) {
                    u[z] <- 0
                  }
                  else if (empty.row == 1) {
                    u[z] <- 1
                  }
                }
                else u[z] <- u[z]/toto
              }
              else if (normalize == "Simpson") {
                toto <- (min(apply(data.spec[samp, ], 
                                   1, sum)))
                if (toto == 0) {
                  if (empty.row == 0) {
                    u[z] <- 0
                  }
                  else if (empty.row == 1) {
                    u[z] <- 1
                  }
                }
                else u[z] <- u[z]/toto
              }
            }
            else {
              if (silent == FALSE) {
                print(paste("warning: number of sites away from the FPO too low to compute zeta for order", 
                            order))
              }
              u[z] <- NA
            }
          }
        }
      }
    }
    zeta.val <- mean(u, na.rm = TRUE)
    if (sd.correct.adapt == FALSE) {
      if (sd.correct == TRUE) {
        zeta.val.sd <- stats::sd(u, na.rm = TRUE)
      }
      else {
        zeta.val.sd <- stats::sd(u, na.rm = TRUE) * sqrt((length(which(!is.na(u))) - 
                                                            1)/length(which(!is.na(u))))
      }
    }
    else {
      if (x > sam) {
        zeta.val.sd <- stats::sd(u, na.rm = TRUE) * sqrt((length(which(!is.na(u))) - 
                                                            1)/length(which(!is.na(u))))
      }
      else {
        zeta.val.sd <- stats::sd(u, na.rm = TRUE)
      }
    }
  }
  zeta.order <- list()
  zeta.order$zeta.order <- order
  zeta.order$combinations <- choose(x, order)
  zeta.order$zeta.val <- zeta.val
  zeta.order$zeta.val.sd <- zeta.val.sd
  zeta.order$zeta.val.vec <- u
  return(zeta.order)
}