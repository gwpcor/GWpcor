gwpcor <-
  function(sdata,
           res_dp,
           vars,
           method = c("pearson", "spearman"),
           kernel = "bisquare",
           adaptive = FALSE,
           bw,
           dMat,
           foreach = FALSE) {

    `%+%` <- function (x, y) {
      paste0(x, y)
    }
    `%!in%` <- function (x, y) {
      !(x %in% y)
    }
    
    ## only sf acceptable
    dist_mat <- function(dp_sf, res_sf){  
      if (missing(dp_sf) || !any(class(dp_sf) == "sf")) 
        stop("Input data points should be sf format")
      if (is.na(st_is_longlat(dp_sf)))
        stop("Input data points should be projected")
      if (missing(res_sf)) 
        res_sf <- dp_sf
      if (!any(class(res_sf) == "sf")) 
        stop("Input resultant points should be sf format")
      if (is.na(st_is_longlat(res_sf)))
        stop("Input resultant points should be projected")
      
      if (isTRUE(st_is_longlat(sdata))){
        dp_locat <- st_coordinates(st_centroid(sdata))
        res_locat <-  st_coordinates(st_centroid(res_sf))
        res_dist <- geodist(dp_locat, res_locat, measure = "cheap") #the mapbox 'cheap' ruler
      }else{
        res_dist <- st_distance(sdata, res_sf) #slow when latlon
      }

      return(res_dist)  
    }
    
    weight_func <- function(type, adapt, dist_vec, bw ){
      
      dist_vec <- as.double(dist_vec)
      
      if(adapt){
        bw_size <- as.integer(length(dist_vec) * bw)
        bw_dist <- as.double(sort(dist_vec)[bw_size])
        
        if(type=="gaussian"){
          weight <- exp((-0.5) * ((dist_vec^2)/(bw_dist^2)))
        }else if(type=="exponential"){
          weight <- exp((-1) * dist_vec/bw_dist)
        }else if(type=="bisquare"){
          weight <-  ifelse((dist_vec > bw_dist), 0 , (1 - (dist_vec^2/bw_dist^2))^2)
        }else if(type=="tricube"){
          weight <-  ifelse((dist_vec > bw_dist), 0 , (1 - (dist_vec^3/bw_dist^3))^3)
        }else if(type=="boxcar"){
          weight <-  ifelse((dist_vec > bw_dist), 0 , 1)
        }
        
      } else{ ##fixed kernel
        if(type=="gaussian"){
          weight <- exp((-0.5) * ((dist_vec^2)/(bw^2)))
        }else if(type=="exponential"){
          weight <- exp((-1) * dist_vec/bw)
        }else if(type=="bisquare"){
          weight <-  ifelse((dist_vec > bw), 0 , (1 - (dist_vec^2/bw^2))^2)
        }else if(type=="tricube"){
          weight <-  ifelse((dist_vec > bw), 0 , (1 - (dist_vec^3/bw^3))^3)
        }else if(type=="boxcar"){
          weight <-  ifelse((dist_vec > bw), 0 , 1)
        }
        
      }  
      return (weight) 
    }
    
    S_Rho <- function(x, y, w) {
      n <- length(w)
      xr <- rank(x)
      yr <- rank(y)
      scorr <- 1 - (6 * w %*% ((xr - yr) ^ 2)) / (n * (n ^ 2 - 1))
    }
    
    #function for p-values
    corr.pval <- function (m, n) {
      r <- cov2cor(m)
      t <- r * sqrt((n - 2) / (1 - r ^ 2))
      p.value <- 2 * pt(-abs(t), (n - 2))
      diag(p.value) <- 0
      return(p.value)
    }
    
    pcorr.pval <- function (m, n) {
      gp <- ncol(m) - 2
      pcor <- corpcor::cor2pcor(m)
      s <- pcor * sqrt((n - 2 - gp) / (1 - pcor ^ 2))
      p.value <- 2 * pt(-abs(s), (n - 2 - gp))
      diag(p.value) <- 0
      return(p.value)
    }
    
    if (missing(vars))
      stop("Variables input error")
    if (length(vars) < 3)
      stop("At least three variables are needed for partial correlation!")
    if (missing(bw) || bw <= 0)
      stop("Bandwidth is not specified correctly")
    len_var <- length(vars)
    
    if (length(vars) == 0)
      stop("Variables input doesn't match with data")
    
    if (is(sdata, "Spatial")) {
      sdata <- st_as_sf(sdata)
      
    } else if (is(sdata, "sf")) {
      
    # nothing
      
    } else{
      stop("Given data must be a Spatial*DataFrame or sf object")}
    
    if (missing(res_dp)) {
        res_dp_given <- FALSE
        res_dp <- sdata
       } else {
        res_dp_given <- T
        if (is(res_dp, "Spatial")) {
          res_dp <- st_as_sf(res_dp)
  
        } else if (is(res_dp, "sf")) {
          ##nothing
         } else {
            stop("res_dp data should be sp or sf format")
         }
        
        if (st_crs(sdata)$proj4string != st_crs(res_dp)$proj4string){
          stop("coordination is not the same.")
        }
      }
    
    dp_n <- dim(sdata)[1]
    res_dp_n <- dim(res_dp)[1]
        
    if (missing(dMat)) {
      
      dMat <- dist_mat(sdata, res_dp)
      
    } else {

      dim.dMat <- dim(dMat)
      if (dim.dMat[1] != dp_n || dim.dMat[2] != res_dp_n)
        stop("Dimensions of dMat are not correct")
    }
    
    ##data extraction
    if (is(sdata, "Spatial")) {
      x <- as(sdata, "data.frame") %>% 
        dplyr::select(vars) %>% 
        as.matrix()
    } else if (is(sdata, "sf")) {
      geom_name <- attr(sdata, "sf_column")
      x <- data.frame(sdata) %>% 
        dplyr::select(-geom_name) %>% 
        dplyr::select(vars) %>% 
        as.matrix()
    }
      
    var_n <- ncol(x)
    
    if (x %>% anyNA())
      stop(" NA values are not allowed")
    if (len_var > var_n)
      warning("Invalid variables have been specified, please check them again!")
    if (method %!in% c("pearson", "spearman")) {
      stop("The method option should be 'pearson' or 'spearman'.")
    }
    
    name.comb.df <- expand.grid(vars[1:(len_var - 1)], vars[2:len_var])
    name.comb_which <-
      apply(name.comb.df, 1, function(x) {
        (x[1] != x[2])
      })
    name.comb <-
      apply(name.comb.df[name.comb_which, ], 1, function(x) {
        x[1] %+% "." %+% x[2]
      }) %>%
      as.vector()
 
    if (foreach == TRUE) {
      
      cl <- makeCluster(detectCores() - 1)
      registerDoParallel(cl)
      
      foreach_out <- foreach(i = 1:res_dp_n) %dopar% {
        
        dist.vi <- dMat[, i] #distance from i
        
        #W.i <- GWmodel::gw.weight(dist.vi, bw, kernel, adaptive)
        W.i <- weight_func(type = kernel, adapt = adaptive, dist_vec =dist.vi, bw =bw)
        sum.w <- sum(W.i)
        Wi <- c(W.i / sum.w)
        
        #calculation
        if (method == "pearson") {
          ans_covwt <-  cov.wt(x, wt = Wi, cor = TRUE)
          cov.matrix <- ans_covwt$cov
          corr.matrix <- ans_covwt$cor
          corr.pval.matrix <-
            corr.pval(ans_covwt$cov, length(Wi[Wi != 0]))
          pcorr.matrix <- corpcor::cor2pcor(ans_covwt$cov)
          pcorr.pval.matrix <-
            pcorr.pval(ans_covwt$cov, length(Wi[Wi != 0]))
          
          cov.mat_i__ <-  cov.matrix[lower.tri(cov.matrix)]
          corr.mat_i__ <- corr.matrix[lower.tri(corr.matrix)]
          corr.pval.mat_i__ <-
            corr.pval.matrix[lower.tri(corr.matrix)]
          pcorr.mat_i__ <- pcorr.matrix[lower.tri(pcorr.matrix)]
          pcorr.pval.mat_i__ <-
            pcorr.pval.matrix[lower.tri(pcorr.pval.matrix)]
          
          out <-   list(
            cov.mat_i__ = cov.mat_i__,
            corr.mat_i__ = corr.mat_i__,
            corr.pval.mat_i__ = corr.pval.mat_i__,
            pcorr.mat_i__ = pcorr.mat_i__,
            pcorr.pval.mat_i__ = pcorr.pval.mat_i__
          )
          
        } else if (method == "spearman") {
          ans_rank_covwt <- cov.wt(apply(x, 2, rank), wt = Wi, cor = TRUE)
          s.cov.matrix <- ans_rank_covwt$cov
          s.corr.matrix <- ans_rank_covwt$cor
          s.corr.pval.matrix <-
            corr.pval(ans_rank_covwt$cor, length(Wi[Wi != 0]))
          s.pcorr.matrix <- corpcor::cor2pcor(ans_rank_covwt$cov)
          s.pcorr.pval.matrix <-
            pcorr.pval(ans_rank_covwt$cor, length(Wi[Wi != 0]))
          
          s.cov.mat_i__ <- s.cov.matrix[lower.tri(s.cov.matrix)]
          s.corr.mat_i__ <- s.corr.matrix[lower.tri(s.corr.matrix)]
          s.corr.pval.mat_i__ <-
            s.corr.pval.matrix[lower.tri(s.corr.pval.matrix)]
          s.pcorr.mat_i__ <- s.pcorr.matrix[lower.tri(s.pcorr.matrix)]
          s.pcorr.pval.mat_i__ <-
            s.pcorr.pval.matrix[lower.tri(s.pcorr.matrix)]
          
          out <-   list(
            s.cov.mat_i__ = s.cov.mat_i__,
            s.corr.mat_i__ = s.corr.mat_i__,
            s.corr.pval.mat_i__ = s.corr.pval.mat_i__,
            s.pcorr.mat_i__ = s.pcorr.mat_i__,
            s.pcorr.pval.mat_i__ = s.pcorr.pval.mat_i__
          )
          
        }
        
        out
        
      }
      
      stopCluster(cl)
      
      
      if (method == "pearson") {
        cov.mat = matrix(t(sapply(foreach_out, "[[", 'cov.mat_i__')), ncol = length(name.comb))
        corr.mat = matrix(t(sapply(foreach_out, "[[", 'corr.mat_i__')), ncol =
                            length(name.comb))
        corr.pval.mat = matrix(t(sapply(
          foreach_out, "[[", 'corr.pval.mat_i__'
        )), ncol = length(name.comb))
        pcorr.mat = matrix(t(sapply(foreach_out, "[[", 'pcorr.mat_i__')), ncol =
                             length(name.comb))
        pcorr.pval.mat =  matrix(t(sapply(
          foreach_out, "[[", 'pcorr.pval.mat_i__'
        )), ncol = length(name.comb))
        
        colnames(cov.mat) <- "cov_" %+% name.comb
        colnames(corr.mat) <- "corr_" %+% name.comb
        colnames(corr.pval.mat) <- "corr_pval_" %+% name.comb
        colnames(pcorr.mat) <- "pcorr_" %+% name.comb
        colnames(pcorr.pval.mat) <- "pcorr_pval_" %+% name.comb
        
        res.df <- data.frame(cov.mat,
                             corr.mat,
                             corr.pval.mat,
                             pcorr.mat,
                             pcorr.pval.mat)
        
      } else if (method == "spearman") {
        s.cov.mat = matrix(t(sapply(foreach_out, "[[", 's.cov.mat_i__')), ncol =
                             length(name.comb))
        s.corr.mat = matrix(t(sapply(
          foreach_out, "[[", 's.corr.mat_i__'
        )), ncol = length(name.comb))
        s.corr.pval.mat = matrix(t(sapply(
          foreach_out, "[[", 's.corr.pval.mat_i__'
        )), ncol = length(name.comb))
        s.pcorr.mat = matrix(t(sapply(
          foreach_out, "[[", 's.pcorr.mat_i__'
        )), ncol = length(name.comb))
        s.pcorr.pval.mat =  matrix(t(sapply(
          foreach_out, "[[", 's.pcorr.pval.mat_i__'
        )), ncol = length(name.comb))
        
        colnames(s.cov.mat) <- "scov_" %+% name.comb
        colnames(s.corr.mat) <- "scorr_" %+% name.comb
        colnames(s.corr.pval.mat) <- "scorr_pval_" %+% name.comb
        colnames(s.pcorr.mat) <- "spcorr_" %+% name.comb
        colnames(s.pcorr.pval.mat) <- "spcorr_pval_" %+% name.comb
        
        res.df <- data.frame(s.cov.mat,
                             s.corr.mat,
                             s.corr.pval.mat,
                             s.pcorr.mat,
                             s.pcorr.pval.mat)
        
      }
      
    } else{
      #non-parallel
      
      if (method == "pearson") {
        cov.mat <- matrix(NA, nrow = res_dp_n, ncol = sum(1:(var_n - 1)))
        corr.mat <- matrix(NA, nrow = res_dp_n, ncol = sum(1:(var_n - 1)))
        corr.pval.mat <-
          matrix(NA, nrow = res_dp_n, ncol = sum(1:(var_n - 1)))
        pcorr.mat <- matrix(NA, nrow = res_dp_n, ncol = sum(1:(var_n - 1)))
        pcorr.pval.mat <-
          matrix(NA, nrow = res_dp_n, ncol = sum(1:(var_n - 1)))
        
      } else if (method == "spearman") {
        s.cov.mat <- matrix(NA, nrow = res_dp_n, ncol = sum(1:(var_n - 1)))
        s.corr.mat <- matrix(NA, nrow = res_dp_n, ncol = sum(1:(var_n - 1)))
        s.corr.pval.mat <-
          matrix(NA, nrow = res_dp_n, ncol = sum(1:(var_n - 1)))
        s.pcorr.mat <- matrix(NA, nrow = res_dp_n, ncol = sum(1:(var_n - 1)))
        s.pcorr.pval.mat <-
          matrix(NA, nrow = res_dp_n, ncol = sum(1:(var_n - 1)))
      }
      
      for (i in 1:res_dp_n) {
        dist.vi <- dMat[, i] #distance from i
        
        #W.i <- gw.weight(dist.vi, bw, kernel, adaptive)
        W.i <- weight_func(type = kernel, adapt = adaptive, dist_vec =dist.vi, bw =bw)
        sum.w <- sum(W.i)
        Wi <- c(W.i / sum.w)
        
        #calculation
        if (method == "pearson") {
          ans_covwt <-  cov.wt(x, wt = Wi, cor = TRUE)
          cov.matrix <- ans_covwt$cov
          corr.matrix <- ans_covwt$cor
          corr.pval.matrix <-
            corr.pval(ans_covwt$cov, length(Wi[Wi != 0]))
          pcorr.matrix <- corpcor::cor2pcor(ans_covwt$cov)
          pcorr.pval.matrix <-
            pcorr.pval(ans_covwt$cov, length(Wi[Wi != 0]))
          
          cov.mat[i,] <-  cov.matrix[lower.tri(cov.matrix)]
          corr.mat[i,] <- corr.matrix[lower.tri(corr.matrix)]
          corr.pval.mat[i,] <-
            corr.pval.matrix[lower.tri(corr.matrix)]
          pcorr.mat[i,] <- pcorr.matrix[lower.tri(pcorr.matrix)]
          pcorr.pval.mat[i,] <-
            pcorr.pval.matrix[lower.tri(pcorr.pval.matrix)]
          
        } else if (method == "spearman") {
          ans_rank_covwt <- cov.wt(apply(x, 2, rank), wt = Wi, cor = TRUE)
          s.cov.matrix <- ans_rank_covwt$cov
          s.corr.matrix <- ans_rank_covwt$cor
          s.corr.pval.matrix <-
            corr.pval(ans_rank_covwt$cor, length(Wi[Wi != 0]))
          s.pcorr.matrix <- corpcor::cor2pcor(ans_rank_covwt$cov)
          s.pcorr.pval.matrix <-
            pcorr.pval(ans_rank_covwt$cor, length(Wi[Wi != 0]))
          
          s.cov.mat[i,] <- s.cov.matrix[lower.tri(s.cov.matrix)]
          s.corr.mat[i,] <- s.corr.matrix[lower.tri(s.corr.matrix)]
          s.corr.pval.mat[i,] <-
            s.corr.pval.matrix[lower.tri(s.corr.pval.matrix)]
          s.pcorr.mat[i,] <- s.pcorr.matrix[lower.tri(s.pcorr.matrix)]
          s.pcorr.pval.mat[i,] <-
            s.pcorr.pval.matrix[lower.tri(s.pcorr.matrix)]
          
        } else {
          stop("The method option should be 'pearson' or 'spearman'.")
        }
        
      }
      
      if (method == "pearson") {
        colnames(cov.mat) <- "cov_" %+% name.comb
        colnames(corr.mat) <- "corr_" %+% name.comb
        colnames(corr.pval.mat) <- "corr_pval_" %+% name.comb
        colnames(pcorr.mat) <- "pcorr_" %+% name.comb
        colnames(pcorr.pval.mat) <- "pcorr_pval_" %+% name.comb
        
        res.df <- data.frame(cov.mat,
                             corr.mat,
                             corr.pval.mat,
                             pcorr.mat,
                             pcorr.pval.mat)
        
      } else if (method == "spearman") {
        colnames(s.cov.mat) <- "scov_" %+% name.comb
        colnames(s.corr.mat) <- "scorr_" %+% name.comb
        colnames(s.corr.pval.mat) <- "scorr_pval_" %+% name.comb
        colnames(s.pcorr.mat) <- "spcorr_" %+% name.comb
        colnames(s.pcorr.pval.mat) <- "spcorr_pval_" %+% name.comb
        
        res.df <- data.frame(s.cov.mat,
                             s.corr.mat,
                             s.corr.pval.mat,
                             s.pcorr.mat,
                             s.pcorr.pval.mat)
      }
      
    }
    
    
    rownames(res.df) <- rownames(res_dp)
    griddedObj <- FALSE
    
    empty_sf <-
      st_sf(id = 1:nrow(res_dp),
            geometry = st_sfc(st_geometry(res_dp)))
    SDF <-
      merge(empty_sf, res.df %>% dplyr::mutate(id = rownames(res.df)))

    res <-
      list(
        SDF = SDF,
        vars = vars,
        kernel = kernel,
        adaptive = adaptive,
        bw = bw
      )
    class(res) <- "gwpcor"
    invisible(res)
  }
