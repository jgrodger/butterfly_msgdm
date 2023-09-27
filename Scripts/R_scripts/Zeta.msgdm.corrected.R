Zeta.msgdm.corrected <- function (data.spec, data.env, xy = NULL, data.spec.pred = NULL, 
          order = 1, sam = 1000, reg.type = "glm", family = stats::gaussian(), 
          method.glm = "glm.fit.cons", cons = -1, cons.inter = 1, 
          confint.level = 0.95, bs = "mpd", kn = -1, order.ispline = 2, 
          kn.ispline = 1, distance.type = "Euclidean", dist.custom = NULL, 
          rescale = FALSE, rescale.pred = TRUE, method = "mean", normalize = FALSE, 
          silent = FALSE, empty.row = 0, control = list(), glm.init = FALSE) 
{
  if (nrow(data.spec) != nrow(data.env)) {
    stop("Error: data.spec and data.env must have the same number of rows.")
  }
  if (!is.null(xy)) {
    if (nrow(data.spec) != nrow(xy) || nrow(data.env) != 
        nrow(xy)) {
      stop("Error: data.spec, data.env and xy must have the same number of rows.")
    }
  }
  if (!inherits(data.spec, "data.frame")) {
    stop(paste("Error: ", deparse(substitute(data.spec)), 
               " is a ", class(data.spec), ". It must be a data frame.", 
               sep = ""))
  }
  if (!inherits(data.env, "data.frame")) {
    stop(paste("Error: ", deparse(substitute(data.env)), 
               " is a ", class(data.env), ". It must be a data frame.", 
               sep = ""))
  }
  if (order > dim(data.spec)[1]) {
    stop("Error: wrong value for \"order\": it must be equal or lower than the number of sites.")
  }
  if (length(which(sapply(data.env, inherits, c("factor", 
                                                "numeric")) == 0)) > 0) {
    stop("Error: variables must be numeric or factor")
  }
  if (order == 1 & (!is.null(xy) | distance.type == "custom")) {
    stop("Error: cannot include distance for order = 1")
  }
  if (silent == FALSE & order == 1 & sum(sapply(data.env, 
                                                inherits, "factor")) > 0) {
    warning("factor variables should be dummy for order = 1")
  }
  if (silent == FALSE & !is.null(dist.custom)) {
    if (!isSymmetric(dist.custom)) {
      warning("Distance matrix is not symmetrical")
    }
  }
  if (empty.row == "remove") {
    if (length(which(rowSums(data.spec) == 0)) > 0) {
      data.env <- data.env[-which(rowSums(data.spec) == 
                                    0), ]
      if (!is.null(xy)) 
        xy <- xy[-which(rowSums(data.spec) == 0), ]
      if (!is.null(data.spec.pred)) 
        if (inherits(data.spec.pred, "data.frame")) 
          data.spec.pred <- data.spec.pred[-which(rowSums(data.spec) == 
                                                    0), ]
      if (inherits(data.spec.pred, "list")) {
        for (p in 1:length(data.spec.pred)) data.spec.pred[[p]] <- data.spec.pred[[p]][-which(rowSums(data.spec) == 
                                                                                                0), ]
      }
      data.spec <- data.spec[-which(rowSums(data.spec) == 
                                      0), ]
    }
  }
  if (reg.type == "ispline") {
    num <- which(sapply(data.env, inherits, "numeric"))
    if (length(num) > 1) {
      range.min <- apply(data.env[, num], 2, min)
      range.max <- apply(data.env[, num], 2, max)
    }
    else {
      range.min <- min(data.env[, num])
      range.max <- max(data.env[, num])
    }
  }
  if (reg.type == "ispline") {
    data.env.num <- as.data.frame(data.env[, which(sapply(data.env, 
                                                          inherits, "numeric"))])
    names(data.env.num) <- names(data.env)[which(sapply(data.env, 
                                                        inherits, "numeric"))]
    ts <- matrix(NA, ncol(data.env.num), 2 * order.ispline + 
                   kn.ispline)
    for (i in 1:ncol(data.env.num)) {
      data.env.num[, i] <- (data.env.num[, i] - min(data.env.num[, 
                                                                 i]))/(max(data.env.num[, i]) - min(data.env.num[, 
                                                                                                                 i]))
      ts[i, ] <- c(rep(0, order.ispline), stats::quantile(data.env.num[, 
                                                                       i], probs = seq(1/(kn.ispline + 1), 1 - 1/(kn.ispline + 
                                                                                                                    1), 1/(kn.ispline + 1))), rep(1, order.ispline))
    }
    IE <- matrix(NA, nrow(data.env.num), (ncol(data.env.num) * 
                                            (order.ispline + kn.ispline)))
    k = order.ispline
    for (j in 1:ncol(data.env.num)) {
      for (i in 1:(order.ispline + kn.ispline)) {
        xx <- 0
        for (x in data.env.num[, j]) {
          xx <- xx + 1
          if (x == 1) {
            IE[xx, (j - 1) * (order.ispline + kn.ispline) + 
                 i] <- 1
          }
          else {
            IE[xx, (j - 1) * (order.ispline + kn.ispline) + 
                 i] <- zetadiv:::.Ii(i, k, x, ts[j, ])
          }
        }
      }
    }
    IE <- data.frame(IE)
    for (i in 1:(ncol(IE)/(order.ispline + kn.ispline))) {
      for (j in 1:(order.ispline + kn.ispline)) {
        names(IE)[(i - 1) * (order.ispline + kn.ispline) + 
                    j] <- paste(names(data.env.num)[i], j, sep = "")
      }
    }
    Fa <- as.data.frame(data.env[, which(sapply(data.env, 
                                                inherits, "factor"))])
    names(Fa) <- names(data.env)[which(sapply(data.env, 
                                              inherits, "factor"))]
    data.env <- cbind(data.frame(IE), Fa)
  }
  x.dim <- dim(data.spec)[1]
  zeta.val <- numeric()
  zeta.val.sd <- numeric()
  if (order == 1) {
    if (nrow(data.spec) < sam) {
      zeta.val <- rowSums(data.spec)
      data.var <- data.env
    }
    else {
      samp <- sample(nrow(data.spec), sam)
      zeta.val <- rowSums(data.spec[samp, ])
      data.var <- data.env[samp, ]
    }
    if (rescale == TRUE || normalize != FALSE) {
      zeta.val <- zeta.val/max(zeta.val)
    }
  }
  else {
    if (choose(x.dim, order) > sam) {
      u <- rep(NA, sam)
      if (!is.null(data.spec.pred)) {
        if (inherits(data.spec.pred, "data.frame")) 
          u2 <- rep(NA, sam)
        if (inherits(data.spec.pred, "list")) {
          u2 <- list()
          for (p in 1:length(data.spec.pred)) u2[[p]] <- rep(NA, 
                                                             sam)
        }
      }
      data.var <- as.data.frame(matrix(NA, sam, dim(data.env)[2]))
      distance <- rep(NA, sam)
      for (z in 1:sam) {
        samp <- sample(1:x.dim, order, replace = FALSE)
        u[z] <- sum(apply(data.spec[samp, ], 2, prod))
        if (!is.null(data.spec.pred)) {
          if (inherits(data.spec.pred, "data.frame")) 
            u2[z] <- sum(apply(data.spec.pred[samp, 
            ], 2, prod))
          if (inherits(data.spec.pred, "list")) {
            for (p in 1:length(data.spec.pred)) u2[[p]][z] <- sum(apply(data.spec.pred[[p]][samp, 
            ], 2, prod))
          }
        }
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
          else {
            u[z] <- u[z]/toto
          }
          if (!is.null(data.spec.pred)) {
            if (inherits(data.spec.pred, "data.frame")) {
              tata <- (ncol(data.spec.pred) - sum(apply((1 - 
                                                           data.spec.pred[samp, ]), 2, prod)))
              if (tata == 0) {
                if (empty.row == 0) {
                  u2[z] <- 0
                }
                else if (empty.row == 1) {
                  u2[z] <- 1
                }
              }
              else {
                u2[z] <- u2[z]/tata
              }
            }
            if (inherits(data.spec.pred, "list")) {
              for (p in 1:length(data.spec.pred)) {
                tata <- (ncol(data.spec.pred[[p]]) - 
                           sum(apply((1 - data.spec.pred[[p]][samp, 
                           ]), 2, prod)))
                if (tata == 0) {
                  if (empty.row == 0) {
                    u2[[p]][z] <- 0
                  }
                  else if (empty.row == 1) {
                    u2[[p]][z] <- 1
                  }
                }
                else {
                  u2[[p]][z] <- u2[[p]][z]/tata
                }
              }
            }
          }
        }
        else if (normalize == "Sorensen") {
          toto <- (mean(apply(data.spec[samp, ], 1, 
                              sum)))
          if (toto == 0) {
            if (empty.row == 0) {
              u[z] <- 0
            }
            else if (empty.row == 1) {
              u[z] <- 1
            }
          }
          else {
            u[z] <- u[z]/toto
          }
          if (!is.null(data.spec.pred)) {
            if (inherits(data.spec.pred, "data.frame")) {
              tata <- (mean(apply(data.spec.pred[samp, 
              ], 1, sum)))
              if (tata == 0) {
                if (empty.row == 0) {
                  u2[z] <- 0
                }
                else if (empty.row == 1) {
                  u2[z] <- 1
                }
              }
              else {
                u2[z] <- u2[z]/tata
              }
            }
            if (inherits(data.spec.pred, "list")) {
              for (p in 1:length(data.spec.pred)) {
                tata <- (mean(apply(data.spec.pred[[p]][samp, 
                ], 1, sum)))
                if (tata == 0) {
                  if (empty.row == 0) {
                    u2[[p]][z] <- 0
                  }
                  else if (empty.row == 1) {
                    u2[[p]][z] <- 1
                  }
                }
                else {
                  u2[[p]][z] <- u2[[p]][z]/tata
                }
              }
            }
          }
        }
        else if (normalize == "Simpson") {
          toto <- (min(apply(data.spec[samp, ], 1, sum)))
          if (toto == 0) {
            if (empty.row == 0) {
              u[z] <- 0
            }
            else if (empty.row == 1) {
              u[z] <- 1
            }
          }
          else {
            u[z] <- u[z]/toto
          }
          if (!is.null(data.spec.pred)) {
            if (inherits(data.spec.pred, "data.frame")) {
              tata <- (min(apply(data.spec.pred[samp, 
              ], 1, sum)))
              if (tata == 0) {
                if (empty.row == 0) {
                  u2[z] <- 0
                }
                else if (empty.row == 1) {
                  u2[z] <- 1
                }
              }
              else {
                u2[z] <- u2[z]/tata
              }
            }
            if (inherits(data.spec.pred, "list")) {
              for (p in 1:length(data.spec.pred)) {
                tata <- (min(apply(data.spec.pred[[p]][samp, 
                ], 1, sum)))
                if (tata == 0) {
                  if (empty.row == 0) {
                    u2[[p]][z] <- 0
                  }
                  else if (empty.row == 1) {
                    u2[[p]][z] <- 1
                  }
                }
                else {
                  u2[[p]][z] <- u2[[p]][z]/tata
                }
              }
            }
          }
        }
        fac <- which(sapply(data.env, inherits, "factor"))
        num <- which(sapply(data.env, inherits, "numeric"))
        if (order > 2) {
          if (length(num) > 1) {
            toto <- data.env[samp, num]
            rownames(toto) <- c()
            data.var[z, num] <- apply(apply(toto, 2, 
                                            stats::dist), 2, get(method))
          }
          else if (length(num) > 0) {
            toto <- data.env[samp, num]
            rownames(toto) <- c()
            data.var[z, num] <- apply(as.matrix(c(stats::dist(toto))), 
                                      2, get(method))
          }
        }
        else {
          if (length(num) > 1) {
            toto <- data.env[samp, num]
            rownames(toto) <- c()
            data.var[z, num] <- apply(toto, 2, stats::dist)
          }
          else if (length(num) > 0) {
            toto <- data.env[samp, num]
            rownames(toto) <- c()
            data.var[z, num] <- stats::dist(toto)
          }
        }
        if (length(fac) > 1) {
          toto <- data.env[samp, fac]
          rownames(toto) <- c()
          data.var[z, fac] <- apply(toto, 2, function(x) {
            length(unique(x))
          }) - 1
        }
        else if (length(fac) > 0) {
          toto <- data.env[samp, fac]
          rownames(toto) <- c()
          data.var[z, fac] <- length(unique(toto)) - 
            1
        }
        if (!is.null(xy)) {
          if (distance.type == "Euclidean") {
            distance[z] <- apply(as.matrix(c(stats::dist(xy[samp, 
            ]))), 2, get(method))
          }
          else if (distance.type == "ortho") {
            distance[z] <- apply(as.matrix(c(.gdist_matrix(xy[samp, 
            ]))), 2, get(method))
          }
          else {
            stop("Error: invalid distance type")
          }
        }
        else if (distance.type == "custom") {
          if (is.null(dist.custom)) {
            stop("Error: a distance matrix must be provided if distance.type = 'custom'")
          }
          distance[z] <- apply(as.matrix(dist.custom[t(utils::combn(sort(samp), 
                                                                    2))]), 2, get(method))
        }
      }
      if (rescale == TRUE & normalize == FALSE) {
        u <- u/ncol(data.spec)
        if (!is.null(data.spec.pred)) {
          if (inherits(data.spec.pred, "data.frame")) 
            u2 <- u2/ncol(data.spec.pred)
          if (inherits(data.spec.pred, "list")) {
            for (p in 1:length(data.spec.pred)) u2[[p]] <- u2[[p]]/ncol(data.spec.pred[[p]])
          }
        }
      }
    }
    else {
      u <- rep(NA, choose(x.dim, order))
      if (!is.null(data.spec.pred)) {
        if (inherits(data.spec.pred, "data.frame")) 
          u2 <- rep(NA, choose(x.dim, order))
        if (inherits(data.spec.pred, "list")) {
          u2 <- list()
          for (p in 1:length(data.spec.pred)) u2[[p]] <- rep(NA, 
                                                             choose(x.dim, order))
        }
      }
      data.var <- as.data.frame(matrix(NA, choose(x.dim, 
                                                  order), dim(data.env)[2]))
      distance <- rep(NA, choose(x.dim, order))
      samp <- utils::combn(1:x.dim, order)
      for (z in 1:dim(samp)[2]) {
        u[z] <- sum(apply(data.spec[samp[, z], ], 2, 
                          prod))
        if (!is.null(data.spec.pred)) {
          if (inherits(data.spec.pred, "data.frame")) 
            u2[z] <- sum(apply(data.spec.pred[samp[, 
                                                   z], ], 2, prod))
          if (inherits(data.spec.pred, "list")) {
            for (p in 1:length(data.spec.pred)) u2[[p]][z] <- sum(apply(data.spec.pred[[p]][samp[, 
                                                                                                 z], ], 2, prod))
          }
        }
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
          else {
            u[z] <- u[z]/toto
          }
          if (!is.null(data.spec.pred)) {
            if (inherits(data.spec.pred, "data.frame")) {
              tata <- (ncol(data.spec.pred) - sum(apply((1 - 
                                                           data.spec.pred[samp[, z], ]), 2, prod)))
              if (tata == 0) {
                if (empty.row == 0) {
                  u2[z] <- 0
                }
                else if (empty.row == 1) {
                  u2[z] <- 1
                }
              }
              else {
                u2[z] <- u2[z]/tata
              }
            }
            if (inherits(data.spec.pred, "list")) {
              for (p in 1:length(data.spec.pred)) {
                tata <- (ncol(data.spec.pred[[p]]) - 
                           sum(apply((1 - data.spec.pred[[p]][samp[, 
                                                                   z], ]), 2, prod)))
                if (tata == 0) {
                  if (empty.row == 0) {
                    u2[[p]][z] <- 0
                  }
                  else if (empty.row == 1) {
                    u2[[p]][z] <- 1
                  }
                }
                else {
                  u2[[p]][z] <- u2[[p]][z]/tata
                }
              }
            }
          }
        }
        else if (normalize == "Sorensen") {
          toto <- (mean(apply(data.spec[samp[, z], ], 
                              1, sum)))
          if (toto == 0) {
            if (empty.row == 0) {
              u[z] <- 0
            }
            else if (empty.row == 1) {
              u[z] <- 1
            }
          }
          else {
            u[z] <- u[z]/toto
          }
          if (!is.null(data.spec.pred)) {
            if (inherits(data.spec.pred, "data.frame")) {
              tata <- (mean(apply(data.spec.pred[samp[, 
                                                      z], ], 1, sum)))
              if (tata == 0) {
                if (empty.row == 0) {
                  u2[z] <- 0
                }
                else if (empty.row == 1) {
                  u2[z] <- 1
                }
              }
              else {
                u2[z] <- u2[z]/tata
              }
            }
            if (inherits(data.spec.pred, "list")) {
              for (p in 1:length(data.spec.pred)) {
                tata <- (mean(apply(data.spec.pred[[p]][samp[, 
                                                             z], ], 1, sum)))
                if (tata == 0) {
                  if (empty.row == 0) {
                    u2[[p]][z] <- 0
                  }
                  else if (empty.row == 1) {
                    u2[[p]][z] <- 1
                  }
                }
                else {
                  u2[[p]][z] <- u2[[p]][z]/tata
                }
              }
            }
          }
        }
        else if (normalize == "Simpson") {
          toto <- (min(apply(data.spec[samp[, z], ], 
                             1, sum)))
          if (toto == 0) {
            if (empty.row == 0) {
              u[z] <- 0
            }
            else if (empty.row == 1) {
              u[z] <- 1
            }
          }
          else {
            u[z] <- u[z]/toto
          }
          if (!is.null(data.spec.pred)) {
            if (inherits(data.spec.pred, "data.frame")) {
              tata <- (min(apply(data.spec.pred[samp[, 
                                                     z], ], 1, sum)))
              if (tata == 0) {
                if (empty.row == 0) {
                  u2[z] <- 0
                }
                else if (empty.row == 1) {
                  u2[z] <- 1
                }
              }
              else {
                u2[z] <- u2[z]/tata
              }
            }
            if (inherits(data.spec.pred, "list")) {
              for (p in 1:length(data.spec.pred)) {
                tata <- (min(apply(data.spec.pred[[p]][samp[, 
                                                            z], ], 1, sum)))
                if (tata == 0) {
                  if (empty.row == 0) {
                    u2[[p]][z] <- 0
                  }
                  else if (empty.row == 1) {
                    u2[[p]][z] <- 1
                  }
                }
                else {
                  u2[[p]][z] <- u2[[p]][z]/tata
                }
              }
            }
          }
        }
        fac <- which(sapply(data.env, inherits, "factor"))
        num <- which(sapply(data.env, inherits, "numeric"))
        if (order > 2) {
          if (length(num) > 1) {
            toto <- data.env[samp[, z], num]
            rownames(toto) <- c()
            data.var[z, num] <- apply(apply(toto, 2, 
                                            stats::dist), 2, get(method))
          }
          else if (length(num) > 0) {
            toto <- data.env[samp[, z], num]
            rownames(toto) <- c()
            data.var[z, num] <- apply(as.matrix(c(stats::dist(toto))), 
                                      2, get(method))
          }
        }
        else {
          if (length(num) > 1) {
            toto <- data.env[samp[, z], num]
            rownames(toto) <- c()
            data.var[z, num] <- apply(toto, 2, stats::dist)
          }
          else if (length(num) > 0) {
            toto <- data.env[samp[, z], num]
            rownames(toto) <- c()
            data.var[z, num] <- stats::dist(toto)
          }
        }
        if (length(fac) > 1) {
          toto <- data.env[samp[, z], fac]
          rownames(toto) <- c()
          data.var[z, fac] <- apply(toto, 2, function(x) {
            length(unique(x))
          }) - 1
        }
        else if (length(fac) > 0) {
          toto <- data.env[samp[, z], fac]
          rownames(toto) <- c()
          data.var[z, fac] <- length(unique(toto)) - 
            1
        }
        if (!is.null(xy)) {
          if (distance.type == "Euclidean") {
            distance[z] <- apply(as.matrix(c(stats::dist(xy[samp[, 
                                                                 z], ]))), 2, get(method))
          }
          else if (distance.type == "ortho") {
            distance[z] <- apply(as.matrix(c(.gdist_matrix(xy[samp[, 
                                                                   z], ]))), 2, get(method))
          }
          else {
            stop("Error: invalid distance type")
          }
        }
        else if (distance.type == "custom") {
          if (is.null(dist.custom)) {
            stop("Error: a distance matrix must be provided if distance.type = 'custom'")
          }
          distance[z] <- apply(as.matrix(dist.custom[t(utils::combn(sort(samp[, 
                                                                              z]), 2))]), 2, get(method))
        }
      }
      if (rescale == TRUE & normalize == FALSE) {
        u <- u/ncol(data.spec)
        if (!is.null(data.spec.pred)) {
          if (inherits(data.spec.pred, "data.frame")) 
            u2 <- u2/ncol(data.spec.pred)
          if (inherits(data.spec.pred, "list")) {
            for (p in 1:length(data.spec.pred)) u2[[p]] <- u2[[p]]/ncol(data.spec.pred[[p]])
          }
        }
      }
    }
    zeta.val <- u
  }
  if (!is.null(xy) | distance.type == "custom") {
    distance.raw <- distance
    d <- max(distance)
  }
  if (rescale.pred == TRUE) {
    if (order > 1) {
      fac <- apply(data.var, 2, max)
      data.var <- data.var/matrix(rep(apply(data.var, 
                                            2, max), min(choose(x.dim, order), sam)), min(choose(x.dim, 
                                                                                                 order), sam), dim(data.env)[2], byrow = T)
    }
    else {
      if (reg.type != "ispline") {
        num <- which(sapply(data.env, inherits, "numeric"))
        if (length(num) > 1) {
          range.min <- apply(data.var[, num], 2, min)
          range.max <- apply(data.var[, num], 2, max)
          data.var[, num] <- (data.var[, num] - matrix(rep(apply(data.var[, 
                                                                          num], 2, min), min(choose(x.dim, order), 
                                                                                             sam)), min(choose(x.dim, order), sam), dim(data.var[, 
                                                                                                                                                 num])[2], , byrow = T))/matrix(rep((apply(data.var[, 
                                                                                                                                                                                                    num], 2, max) - apply(data.var[, num], 2, 
                                                                                                                                                                                                                          min)), min(choose(x.dim, order), sam)), 
                                                                                                                                                                                min(choose(x.dim, order), sam), dim(data.var[, 
                                                                                                                                                                                                                             num])[2], byrow = T)
        }
        else {
          range.min <- min(data.var[, num])
          range.max <- max(data.var[, num])
          data.var[, num] <- (data.var[, num] - min(data.var[, 
                                                             num]))/(max(data.var[, num]) - min(data.var[, 
                                                                                                         num]))
        }
      }
    }
    if (!is.null(xy) | distance.type == "custom") {
      distance <- distance/max(distance)
    }
  }
  names(data.var) <- names(data.env)
  if ((!is.null(xy) | distance.type == "custom") & reg.type == 
      "ispline") {
    dist2 <- distance/max(distance)
    ts <- c(rep(0, order.ispline), stats::quantile(dist2, 
                                                   probs = seq(1/(kn.ispline + 1), 1 - 1/(kn.ispline + 
                                                                                            1), 1/(kn.ispline + 1))), rep(1, order.ispline))
    IE <- matrix(NA, length(distance), (order.ispline + 
                                          kn.ispline))
    k = order.ispline
    for (i in 1:(order.ispline + kn.ispline)) {
      xx <- 0
      for (x in dist2) {
        xx <- xx + 1
        if (x == 1) {
          IE[xx, i] <- 1
        }
        else {
          IE[xx, i] <- zetadiv:::.Ii(i, k, x, ts)
        }
      }
    }
    distance <- data.frame(IE)
    for (i in 1:(order.ispline + kn.ispline)) {
      names(distance)[i] <- paste("distance", i, sep = "")
    }
  }
  if (reg.type == "ispline") {
    if (!is.null(data.spec.pred)) {
      if (inherits(data.spec.pred, "data.frame")) {
        sp <- 1 - u2
        spp <- matrix(sp, length(sp), 1)
        ts <- c(rep(0, order.ispline), stats::quantile(sp, 
                                                       probs = seq(1/(kn.ispline + 1), 1 - 1/(kn.ispline + 
                                                                                                1), 1/(kn.ispline + 1))), rep(1, order.ispline))
        IE <- matrix(NA, length(u2), (order.ispline + 
                                        kn.ispline))
        k = order.ispline
        for (i in 1:(order.ispline + kn.ispline)) {
          xx <- 0
          for (x in sp) {
            xx <- xx + 1
            if (x == 1) {
              IE[xx, i] <- 1
            }
            else {
              IE[xx, i] <- zetadiv:::.Ii(i, k, x, ts)
            }
          }
        }
        sp.prev <- data.frame(IE)
        for (i in 1:(order.ispline + kn.ispline)) {
          names(sp.prev)[i] <- paste("Biotic", i, sep = "")
        }
      }
      if (inherits(data.spec.pred, "list")) {
        for (p in 1:length(data.spec.pred)) {
          sp <- 1 - u2[[p]]
          if (p == 1) {
            spp <- matrix(sp, length(sp), 1)
          }
          else {
            spp <- cbind(spp, sp)
          }
          ts <- c(rep(0, order.ispline), stats::quantile(sp, 
                                                         probs = seq(1/(kn.ispline + 1), 1 - 1/(kn.ispline + 
                                                                                                  1), 1/(kn.ispline + 1))), rep(1, order.ispline))
          IE <- matrix(NA, length(u2[[p]]), (order.ispline + 
                                               kn.ispline))
          k = order.ispline
          for (i in 1:(order.ispline + kn.ispline)) {
            xx <- 0
            for (x in sp) {
              xx <- xx + 1
              if (x == 1) {
                IE[xx, i] <- 1
              }
              else {
                IE[xx, i] <- zetadiv:::.Ii(i, k, x, ts)
              }
            }
          }
          if (p == 1) {
            sp.prev <- data.frame(IE)
          }
          else {
            sp.prev <- cbind(sp.prev, data.frame(IE))
          }
        }
        pp <- 0
        for (p in 1:length(data.spec.pred)) {
          for (i in 1:(order.ispline + kn.ispline)) {
            pp <- pp + 1
            names(sp.prev)[pp] <- paste("Biotic_", p, 
                                        "_", i, sep = "")
          }
        }
      }
    }
  }
  else {
    if (!is.null(data.spec.pred)) {
      if (inherits(data.spec.pred, "data.frame")) {
        sp.prev <- data.frame(1 - u2)
        names(sp.prev) <- c("Biotic")
      }
      if (inherits(data.spec.pred, "list")) {
        for (p in 1:length(data.spec.pred)) {
          if (p == 1) {
            sp.prev <- data.frame(1 - u2[[p]])
          }
          else {
            sp.prev <- cbind(sp.prev, data.frame(1 - 
                                                   u2[[p]]))
          }
        }
        for (p in 1:length(data.spec.pred)) {
          names(sp.prev)[p] <- paste("Biotic", p, sep = "")
        }
      }
    }
  }
  # if (is.null(xy) & distance.type != "custom" & is.null(data.spec.pred)) {
  #   data.tot <- data.var
  # }
  # else if (is.null(xy) & distance.type != "custom" & !is.null(data.spec.pred)) {
  #   data.tot <- cbind(data.var, sp.prev)
  # }
  # else if (!is.null(xy) & distance.type != "custom" & is.null(data.spec.pred)) {
  #   data.tot <- cbind(data.var, distance)
  # }
  # else {
  #   data.tot <- cbind(data.var, sp.prev, distance)
  # }
  if (is.null(xy) & distance.type != "custom" & is.null(data.spec.pred)) {
    data.tot <- data.var
  }
  else if (is.null(xy) & distance.type != "custom" & !is.null(data.spec.pred)) {
    data.tot <- cbind(data.var, sp.prev)
  }
  else if ((!is.null(xy) | distance.type == "custom") & is.null(data.spec.pred)) {
    data.tot <- cbind(data.var, distance)
  }
  else {
    data.tot <- cbind(data.var, sp.prev, distance)
  }
  zeta.msgdm <- list()
  zeta.msgdm$val <- zeta.val
  zeta.msgdm$predictors <- data.tot
  if (reg.type == "ispline") {
    zeta.msgdm$range.min <- range.min
    zeta.msgdm$range.max <- range.max
    if (!is.null(data.spec.pred)) 
      zeta.msgdm$biotic <- spp
    if (!is.null(xy) | distance.type == "custom") 
      zeta.msgdm$distance <- distance.raw
  }
  if (rescale.pred == TRUE) {
    if (order > 1) {
      if (!is.null(xy) | distance.type == "custom") {
        zeta.msgdm$rescale.factor <- c(fac, d)
      }
      else {
        zeta.msgdm$rescale.factor <- fac
      }
    }
    else {
      if (reg.type != "ispline") {
        num <- which(sapply(data.env, inherits, "numeric"))
        if (length(num) > 1) {
          zeta.msgdm$range.min <- range.min
          zeta.msgdm$range.max <- range.max
        }
      }
    }
  }
  zeta.msgdm$my.order <- order
  zeta.msgdm$order.ispline <- order.ispline
  zeta.msgdm$kn.ispline <- kn.ispline
  if (reg.type == "glm") {
    if (method.glm == "glm.fit.cons") {
      zeta.msgdm.model <- glm.cons(zeta.val ~ ., data = data.tot, 
                                   family = family, method = method.glm, cons = cons, 
                                   cons.inter = cons.inter, control = control)
    }
    else {
      zeta.msgdm.model <- glm2::glm2(zeta.val ~ ., data = data.tot, 
                                     family = family, method = method.glm, control = control)
      zeta.msgdm.confint <- suppressMessages(stats::confint(zeta.msgdm.model, 
                                                            level = confint.level))
    }
    if (dim(data.env)[2] > 1) {
      zeta.msgdm.vif <- car::vif(zeta.msgdm.model)
    }
    else {
      zeta.msgdm.vif <- NA
    }
    zeta.msgdm$model <- zeta.msgdm.model
    if (method.glm == "glm.fit2") 
      zeta.msgdm$confint <- zeta.msgdm.confint
    zeta.msgdm$vif <- zeta.msgdm.vif
  }
  else if (reg.type == "ngls") {
    data.tot2 <- cbind(rep(1, nrow(data.tot)), data.tot)
    start <- c(1, rep(-1, ncol(data.tot)))
    zeta.msgdm$model <- nnls::nnnpls(as.matrix(data.tot2), 
                                     zeta.val, con = start)
  }
  else if (reg.type == "gam") {
    xnam <- names(data.tot)
    fm <- stats::as.formula(paste("zeta.val ~ s(", paste(xnam, 
                                                         collapse = paste(", k = ", kn, ") + s(", sep = ""), 
                                                         sep = ""), ", k = ", kn, ")", sep = ""))
    zeta.msgdm$model <- mgcv::gam(fm, data = data.tot, family = family)
  }
  else if (reg.type == "scam") {
    xnam <- names(data.tot)
    fm <- stats::as.formula(paste("zeta.val ~ s(", paste(xnam, 
                                                         collapse = paste(", k = ", kn, ", bs = '", bs, "') + s(", 
                                                                          sep = ""), sep = ""), ", k = ", kn, ",bs='", 
                                  bs, "')", sep = ""))
    zeta.msgdm$model <- scam::scam(fm, data = data.tot, 
                                   family = family)
  }
  else if (reg.type == "ispline") {
    if (glm.init == TRUE) {
      tutu <- stats::cor(data.tot, zeta.val)
      tutu[which(tutu > 0)] <- 0
      zeta.msgdm$model <- glm.cons(zeta.val ~ ., data = data.tot, 
                                   family = family, method = "glm.fit.cons", cons = cons, 
                                   cons.inter = cons.inter, control = control, 
                                   start = c(-1, tutu))
    }
    else {
      zeta.msgdm$model <- glm.cons(zeta.val ~ ., data = data.tot, 
                                   family = family, method = "glm.fit.cons", cons = cons, 
                                   cons.inter = cons.inter, control = control)
    }
  }
  else {
    stop("Error: unknown regression type.")
  }
  return(zeta.msgdm)
}
