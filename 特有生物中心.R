
load("C:/Users/User/Downloads/特有生物中心/Incidence_data.RData")
library(glmnet)
library(compiler)

## 計算物種數的
r.score = function(x) {
  return( ChaoRichness(x,datatype='incidence_raw')$Estimator)
}
r.score = cmpfun(r.score)

#計算覆蓋度
cv.score = function(x) {
  return(DataInfo(x, datatype = "incidence_raw")$SC)
}
cv.score = cmpfun(cv.score)

## training set optimization 
opt.train = function(X,n,iter_time) {
  
  ## initual condition
  N = ncol(X)
  Nc = N 
  cand = seq(N)
  lambda = 1
  
  ## set solution size 
  if ( n >= 10) { 
    ss = n
  } else {
    ss = 10
  }
  
  ## solution matrix 
  sol = c()
  for ( i in 1:ss) {
    sol = cbind(sol, sample(c(rep(1, n), rep(0, Nc - n))))
  }
  # random.set = cand[which(sol[ ,1] == 1)]
  
  ## start GA
  iter = 0
  stop = 0
  while(stop == 0) {
    iter = iter + 1
    cat(iter, "...", sep = "")
    
    ## r_score
    score = c()
    if (iter == 1) {
      max.score = c()
    }
    for (i in seq(ss)) {
      #calculate r-score 
      score = c(score, r.score(X[,cand[which(sol[ ,i] == 1)]]))
    }
    max.score = c(max.score, max(score)[1])
    
    ## stop criteria
    if (iter >= iter_time) {
      #if (max(score) - max.score[iter - 100] < threshold) {
      stop = stop + 1
      #}
    }
    
    ## elite solution(挑出菁英)
    elite = which(rank(-score, ties.method = "random") <= floor(n/sqrt(ss)))
    cat(which(sol[, elite ] == 1), "\n")
    
    ## Delete solution(刪除兩個)
    del = sample(seq(ss)[-elite], 2, prob = (1/score[-elite]) / sum(1/score[-elite]))
    
    ## crossover(選兩個交叉)
    for (i in del) {
      chr = sample(seq(ss)[-del], 2)
      pos = sample(seq(Nc - 1), 1)
      sol[ ,i] = c(sol[1:pos, chr[1]], sol[(pos + 1):Nc, chr[2]])
    }
    
    ## mutation(突變，菁英集不變，其他隨機突變，2個1變0，2個0變1，
    #如果突變完比舊的還高就取代)
    for ( i in seq(ss)) {
      sol.new = sol[, i]
      n.sol = length(which(sol[, i] == 1))
      
      if (n.sol == n) {
        pos = c(sample(which(sol[, i] == 0), 2), sample(which(sol[, i] == 1), 2))
      } else if (n.sol > n) {
        pos = c(sample(which(sol[, i] == 0), 2), sample(which(sol[, i] == 1), (2 + n.sol - n)))
      } else if (n.sol < n) {
        pos = c(sample(which(sol[, i] == 0), (2 - n.sol + n)), sample(which(sol[, i] == 1), 2))
      }
      
      sol.new[pos] = abs(sol.new[pos] - 1)
      
      if (!(i %in% elite)) {
        sol[, i] = sol.new
      } else {
        old = r.score(X[,cand[which(sol[, i] == 1)]])
        new = r.score(X[,cand[which(sol.new == 1)]])
        if (new > old) {
          sol[, i] = sol.new
        }
      }
    }
    
    if (stop != 0) {
      cat("\n Genatic Algorithm ended!\n")
    }
  } ## END GA
  
  ## OPT Training set 
  sol = sol[, order(score, decreasing = TRUE)[1:5]] # Get top 5 solutions
  opt.set = lapply(1:5, function(i) cand[which(sol[, i] == 1)])
  
  return(opt = opt.set)
}
opt.train = cmpfun(opt.train)

library(iNEXT)

result=list()
result_nos=list()
result_cv=list()
add=matrix(0,26)

for (i in 1:26) {
  dat <- incidence_dataset[[i]]$`1 魚類`
  riv <- substring(colnames(incidence_dataset[[i]][["1 魚類"]]), 1, 6)
  num_1 <- length(unique(riv))
  out <- estimateD(
    dat,
    q = 0,
    datatype = 'incidence_raw',
    base = "coverage",
    level = 0.9,
    nboot = 0,
    conf = 0.95
  )
  
  n <- ncol(dat)
  best <- list()
  nos <- list()
  cv <- list()
  
  if (out$t > n) {
    add[i] <- floor(out$t - n)
    out$t <- ncol(dat)
    best <- colnames(dat)
    nos <- r.score(dat)
    cv <- cv.score(dat)
  } else if (choose(n, ceiling(out$t)) > 10^5) {
    ga <- opt.train(X = dat, n = ceiling(out$t), iter_time = 1000)
    
    for (k in 1:5) {
      if (k == 1) { t <- 1 }
      num_2 <- length(unique(substring(colnames(dat[, ga[[k]]]), 1, 6)))
      if (num_1 == num_2) {
        best[[t]] <- colnames(dat[, ga[[k]]])
        nos[[t]] <- r.score(dat[, ga[[k]]])
        cv[[t]] <- cv.score(dat[, ga[[k]]])
        t <- t + 1
      }
    }
  } else {
    all <- combn(n, ceiling(out$t))
    r_all <- c(0, 0, 0, 0, 0)
    for (j in 1:ncol(all)) {
      if (r.score(dat[, all[, j]]) > min(r_all)) {
        loc <- which.min(r_all)
        r_all[loc] <- r.score(dat[, all[, j]])
        names(r_all)[loc] <- j
        r_all <- sort(r_all, decreasing = TRUE)
      }
    }
    for (k in 1:5) {
      if (k == 1) { t <- 1 }
      num_2 <- length(unique(substring(colnames(dat)[all[, as.numeric(names(r_all)[k])]], 1, 6)))
      if (num_1 == num_2) {
        best[[t]] <- colnames(dat)[all[, as.numeric(names(r_all)[k])]]
        nos[[t]] <- r.score(dat[, all[, as.numeric(names(r_all)[k])]])
        cv[[t]] <- cv.score(dat[, all[, as.numeric(names(r_all)[k])]])
        t <- t + 1
      }
    }
  }
  
  result_cv[[i]] <- cv
  result[[i]] <- best
  result_nos[[i]] <- nos
}

miss <- c()
for (i in 1:26) {
  if ((length(result[[i]]) == 0)) {
    miss <- c(miss, i)
  }
}

for (i in miss) {
  dat <- incidence_dataset[[i]]$`1 魚類`
  riv <- substring(colnames(incidence_dataset[[i]][["1 魚類"]]), 1, 6)
  num_1 <- length(unique(riv))
  out <- estimateD(
    dat,
    q = 0,
    datatype = 'incidence_raw',
    base = "coverage",
    level = 0.9,
    nboot = 0,
    conf = 0.95
  )
  
  n <- ncol(dat)
  best <- list()
  nos <- list()
  cv <- list()
  
  if (choose(n, ceiling(out$t)) > 10^5) {
    ga <- opt.train(X = dat, n = ceiling(out$t), iter_time = 1000)
    
    for (k in 1:5) {
      if (k == 1) { t <- 1 }
      a <- unique(riv)[which(!unique(riv) %in% unique(substring(colnames(dat[, ga[[k]]]), 1, 6)))]
      b <- list()
      for (l in 1:length(a)) {
        b[[l]] <- grep(a[l], colnames(dat))
      }
      num_3 <- expand.grid(b)
      r_max <- 0
      for (m in 1:nrow(num_3)) {
        dat2 <- cbind(dat[, ga[[k]]], dat[, as.numeric(num_3[m, ])])
        if (r.score(dat2) > r_max) {
          best[[t]] <- colnames(dat2)
          cv[[t]] <- cv.score(dat2)
          nos[[t]] <- r_max <- r.score(dat2)
        }
      }
      t <- t + 1
    }
  } else {
    all <- combn(n, ceiling(out$t))
    r_all <- c(0, 0, 0, 0, 0)
    for (j in 1:ncol(all)) {
      if (r.score(dat[, all[, j]]) > min(r_all)) {
        loc <- which.min(r_all)
        r_all[loc] <- r.score(dat[, all[, j]])
        names(r_all)[loc] <- j
        r_all <- sort(r_all, decreasing = TRUE)
      }
    }
    for (k in 1:5) {
      if (k == 1) { t <- 1 }
      a <- unique(riv)[which(!unique(riv) %in% unique(substring(colnames(dat)[all[, as.numeric(names(r_all)[k])]], 1, 6)))]
      b <- list()
      for (l in 1:length(a)) {
        b[[l]] <- grep(a[l], colnames(dat))
      }
      num_3 <- expand.grid(b)
      r_max <- 0
      for (m in 1:nrow(num_3)) {
        dat2 <- cbind(dat[, ga[[k]]], dat[, as.numeric(num_3[m, ])])
        if (r.score(dat2) > r_max) {
          best[[t]] <- colnames(dat2)
          cv[[t]] <- cv.score(dat2)
          nos[[t]] <- r_max <- r.score(dat2)
        }
      }
      t <- t + 1
    }
  }
  
  result_cv[[i]] <- cv
  result[[i]] <- best
  result_nos[[i]] <- nos
}

save(result,result_cv,add,result_nos file = 'result_fish')









