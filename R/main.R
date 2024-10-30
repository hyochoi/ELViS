#' Sample filtering threshold examination plot.
#'
#' @param mtrx Matrix or data.frame. Rows are positions and columns are samples.
#' @param th Numeric. Sample filtering threshold
#' @param title_txt figure title.
#' @param smry_fun function to calculate summary metric to apply sample filter threshold to
#' @param ... additional argument for `smry_fun` argument.
#'
#' @return ggplot2 object
#' @export
#' @import ggplot2
#' @rawNamespace import(glue, except=c(trim))
#' @import stringr
#' @import scales
#' @import knitr
#'
#' @examples
#' data(mtrx_samtools_reticulate)
#' th <- 50
#' depth_hist(mtrx_samtools_reticulate,th=th,smry_fun=max)
#' depth_hist(mtrx_samtools_reticulate,th=th,smry_fun = quantile,prob=0.95)
#' depth_hist(mtrx_samtools_reticulate,th=th,smry_fun = mean)
depth_hist <- function(mtrx,th=50,title_txt=NULL,smry_fun=max,...){
  if(is.null(title_txt)){
    dot_info <- capture_params_glue(...)
    fun_name <- str_to_title(deparse(substitute(smry_fun)))
    if(length(dot_info)==0){
      title_txt <- (glue("{fun_name} Depth Distribution"))
    }else{
      title_txt <- (glue("{fun_name} {dot_info} Depth Distribution"))
    }

  }

  ggplot(
    data.frame(Max_Depth = apply(mtrx,2,smry_fun,...))
  ) +
    aes(x = .data$Max_Depth) +
    geom_histogram(col="gray30",fill="gray80",) +
    geom_rect(aes(xmin = 0, xmax = th, ymin = 0, ymax = Inf),
              fill = alpha("gray",0.01)) +
    geom_vline(xintercept = th,col="red",linetype="dashed") +
    scale_x_continuous(trans="log10") +
    ggtitle(title_txt)

}


capture_params_glue <- function(...) {
  # Capture the names and values of the arguments
  args <- as.list(match.call())[-1]  # Exclude the function name from the list

  # Create a formatted string using glue and paste
  formatted_params <- paste0(names(args), " = ", vapply(args, deparse,""), collapse = " ")

  # Use glue to return the formatted string
  if(length(args)==0){
    return(NULL)
  }else{
    return(glue("{formatted_params}"))
  }

}



#' @import scales
#' @noRd
log10p1_trans <- scales::trans_new(
  name = "log10p1",
  transform = function(x) log10(x + 1),
  inverse = function(x) 10^x - 1,
  breaks = function(x) log_breaks(base = 10)(x + 1)
)


#' Filtering samples based on summary statistic
#'
#' @param mtrx matrix or data.frame. Rows are positions and columns are samples.
#' @param th Sample filtering threshold (Default : 50)
#' @param smry_fun function to generate summary value of samples, which is used for filtering. (Default : `max`)
#'
#' @return matrix or data.frame. data matrix with low depth samples filtered out.
#' @export
#'
#' @examples
#' data(mtrx_samtools_reticulate)
#' th<-50
#' filtered_mtrx <- filt_samples(mtrx_samtools_reticulate,th=th,smry_fun=max)
filt_samples <- function(mtrx,th=50,smry_fun=max){
  mtrx[,apply(mtrx,2,smry_fun) > th]
}














#' Normalization - scaling by certain quantile
#'
#' @param x numeric vector to normalize.
#' @param probs a single numeric value of probabilities in `[0,1]` used for normalization.(Default = `0.75`)
#'
#' @return numeric vector of normalized values
#'
#' @export
#' @examples
#'
#' norm.fun(seq_len(5))
#' # [1] 0.25 0.50 0.75 1.00 1.25
#'
norm.fun <- function(x,probs=0.75) {
  m <- quantile(x,probs=probs)
  if (m>0) {
    return(x/m)
  } else {
    return(rep(0,length(x)))
  }
}


#' get windows cut at points where depth profile crosses upper quartile
#'
#' @param Y numeric matrix containing normalized data
#' @param sam sample indexindex
#'
#' @return window list
#' @noRd
get_window_v1 <- function(Y,sam) {
  ips <- which((Y[,sam]-1)[seq_len(nrow(Y)-1)]*(Y[,sam]-1)[2:nrow(Y)]<=0)
  win <- matrix(c(1,rep(ips,each=2),dim(Y)[1]),ncol=2,byrow=TRUE)
  win <- win[which((win[,2] - win[,1])>10),]
  return(win)
}


#' get average per window and sample
#'
#' @param Y numeric matrix containing normalized data
#' @param sam sample indexindex
#'
#' @return list containing matrices of averages per window and sample
#' @noRd
get_outmat <- function(Y,sam) {
  win <- get_window(Y,sam)
  outmat <- t(apply(win,1,function(x) pd.rate.hy(apply(Y[c(x[1]:x[2]),],2,mean),qrsc = TRUE)))
  return(outmat)
}




#' get average per window and sample - another version
#'
#' @param Y numeric matrix containing normalized data
#' @param win window matrix obtained by running get_window_v1
#'
#' @return list containing matrices of averages per window and sample
#' @noRd
get_outmat_v1 <- function(Y,win) {
  outmat <- t(apply(win,1,function(x) pd.rate.hy(apply(Y[c(x[1]:x[2]),],2,mean),qrsc = TRUE)))
  return(outmat)
}







#' Function to normalize raw coverages
#'
#' @param X numeric matrix with raw coverages
#' @param q.cutoff a numeric value. This is used as a cut.off to exclude extremely outlying windows when medians are calculated.
#'
#' @return list of two matrices Y and Z. Y and Z contain normalized and standardized data, respectively.
#' @noRd
get_normalized_data <- function(X,q.cutoff=4){
  n <- NCOL(X)
  d <- NROW(X)
  # upper quartile normalizaiton
  Y <- apply(X,2,norm.fun)

  # get windows(intervals) cut at UQ-crossing points
  window.list <- lapply(c(seq_len(ncol(Y))),function(j) get_window_v1(Y=Y,sam=j)) # get windows
  # get average per window and sample, standardize across samples to get relative position in the cohort
  outmat.list <- lapply(window.list,function(j) get_outmat_v1(Y=Y,win=j)) # get per-window mean for each sample and get relative score of it, upper and lower part SDs are different


  Y2 <- matrix(0,ncol=ncol(Y),nrow=nrow(Y))
  colnames(Y2) <- colnames(Y)
  for (sam in seq_len(n) ) {
    # segments cut based on crossing over UQ and deviation from the center of entire cohort
    sammat <- cbind(window.list[[sam]],outmat.list[[sam]][,sam])
    outreg <- c()

    # when calculating median, exclude extreme outlier segments (compared to other samples not compared to other segments) for low coverage
    j <- which.min(sammat[,3])
    if (sammat[j,3]<(-q.cutoff)) {
      outreg <- c(sammat[j,1]:sammat[j,2])
    }
    # renoramlize based on median
    Y2[,sam] <- Y[,sam]/median(Y[which(! c(seq_len(d)) %in% outreg),sam])
  }

  Z2 <- t(apply(Y2,1,function(x) pd.rate.hy(x,qrsc=TRUE))) # stabilize deviations seperately for higher than median and lower than median

  # final normalized data
  output <- list(sample_Ids=colnames(X),Y=Y2,Z=Z2)

  return(output)
}




#' Updating reference segments
#'
#' @param normalized_data list containing 1 character vector and 2 numeric matrices. The eirst element is sample IDs, the second element contain normalized depth and the third element contain standardized depth.
#' @param cores integer; number of threads to use for the analysis. (Default : 10L)
#' @param verbose logical whether to print out information for debugging
#'
#' @return a list containing updated data
#' @import parallel
#' @import segclust2d
#' @rawNamespace import(circlize, except=c(degree))
#' @rawNamespace import(zoo, except=c(index,yearmon,yearqtr,"index<-"))
#' @noRd
update_reference_segments <- function(normalized_data,cores=10L,verbose=FALSE){

  n <- NCOL(normalized_data$Y)
  d <- NROW(normalized_data$Y)

  sample_id_vec <- normalized_data$sample_Ids
  Y <- normalized_data$Y
  Z <- normalized_data$Z
  n <- length(sample_id_vec)

  # 1 get initial segmentation
  # time taken 10 cores : 6 sec
  #   user  system elapsed
  # 37.262   6.004   6.163
  if(verbose) message("    segment.K_initial")
  segment.K_initial <-
    seq_len(n) %>%
    mclapply(
      mc.cores=cores,

      function(sam) {
        message(sam)
        sampleID <- sample_id_vec[sam]
        testdata <- data.frame(z=Z[,sam],y=Y[,sam])
        output <- tryCatch({

          shift_seg <- segmentation(testdata,
                                    lmin=300,
                                    Kmax = 10,
                                    seg.var = c("z","y"), subsample_by = 60, scale.variable = FALSE)
          K <- shift_seg$Kopt.lavielle

          msg <- paste0(sam,"| done"); message(msg)

          K

        },error=function(err) {
          msg <- paste0(sam,"|",err); message(msg)
          message("Running 1D segmentation with robust-scaled data only")

          # proceed only with Z
          shift_seg <- segmentation(testdata,
                                    lmin=300,
                                    Kmax = 10,
                                    seg.var = c("z"), subsample_by = 60, scale.variable = FALSE)
          K <- shift_seg$Kopt.lavielle

          msg <- paste0(sam,"| done"); message(msg)

          #output
          K

        }) # end of tryCatch

        return(output)
      }
    ) %>%
    unlist



  # 2 recentering data, setting least deviating segment as reference level
  # time taken 10 cores : 21.765 sec
  #    user  system elapsed
  # 120.368   6.716  21.765
  if(verbose) message("    Y_recentered")
  Y_recentered <-
    seq_len(n) %>%
    mclapply(mc.cores=cores,

             function(sam) {
               sampleID <- colnames(Y)[sam]
               testdata <- data.frame(z=Z[,sam],y=Y[,sam])
               output <- tryCatch({


                 K <- segment.K_initial[sam]
                 if (K>1) {
                   clust_seg <- segclust(testdata, lmin=300, Kmax=10, ncluster = (2:K),
                                         seg.var = c("z","y"), scale.variable = FALSE, subsample_by = 60)
                   result <- segment(clust_seg)

                   # Update mean
                   # Setting the segment that shows the least deviation from median(=1, since the data is median-normalized)
                   # as new reference segment
                   mean.update <- result$mu.y[which.min(abs(result$mu.y-1))]
                   new_y <- testdata$y/mean.update
                 } else {
                   new_y <- testdata$y
                 }
                 msg <- paste0(sam,"| done"); message(msg)

                 #output
                 new_y
               },error=function(err) {
                 msg <- paste0(sam,"|",err); message(msg)
                 message("Running 1D segmentation with robust-scaled data only")


                 K <- segment.K_initial[sam]
                 if (K>1) {
                   clust_seg <- segclust(testdata, lmin=300, Kmax=10, ncluster = (2:K),
                                         seg.var = c("z"), scale.variable = FALSE, subsample_by = 60)
                   result <- segment(clust_seg)

                   result$mu.y <-
                     apply(
                       result[,c("begin","end")],
                       1,
                       \(crd){
                         mean(testdata$y[crd[1]:crd[2]])
                       }
                     )

                   mu.y_nonzero <- result$mu.y[result$mu.y!=0]

                   # Update mean
                   # Setting the segment that shows the least deviation from median(=1, since the data is median-normalized)
                   # as new reference segment
                   mean.update <- mu.y_nonzero[which.min(abs(mu.y_nonzero-1))]
                   new_y <- testdata$y/mean.update
                 } else {
                   new_y <- testdata$y
                 }
                 msg <- paste0(sam,"| done"); message(msg)

                 #output
                 new_y

               }) # end of tryCatch

               return(output)
             }
    ) %>%
    simplify2array

  Z_recentered <- t(apply(Y_recentered,1,function(x) pd.rate.hy(x,qrsc=TRUE)))

  return(
    list(
      sample_Ids = normalized_data$sample_Ids,
      Y_recentered = Y_recentered,
      Z_recentered = Z_recentered,
      segment.K_initial = segment.K_initial
    )
  )


}











#' get segments and clusters from ref-updated data
#'
#' @param refupate_data list. Output of function [update_reference_segments]. list containing 2 numeric matrices. First element contain normalized depth and the second element contain standardized depth.
#' @param cores integer; number of threads to use for the analysis. (Default : 10L)
#'
#' @return a list containing segment and cluster information
#' @noRd
get_segments_and_clusters <- function(refupate_data,cores=10L){
  sample_id_vec <- refupate_data$sample_Ids
  Y_recentered <- refupate_data$Y_recentered
  Z_recentered <- refupate_data$Z_recentered
  segment.K_initial <- refupate_data$segment.K_initial
  n <- length(sample_id_vec)
  d <- NROW(Y_recentered)

  seg_result <-
    seq_len(n) %>%
    mclapply(
      mc.cores=cores,

      function(sam) {
        sampleID <- sample_id_vec[sam]
        testdata <- data.frame(z=Z_recentered[,sam],y=Y_recentered[,sam])
        output <- tryCatch({

          K <- segment.K_initial[sam]
          if (K>1) {

            shift_seg <- segmentation(testdata, lmin=300, Kmax = 10, seg.var = c("z","y"), subsample_by = 60, scale.variable = FALSE)
            K <- shift_seg$Kopt.lavielle

            # K can be 1 here
            if(K>1){
              clust_seg <- segclust(testdata, lmin=300, Kmax=10, ncluster = (2:K),
                                    seg.var = c("z","y"), scale.variable = FALSE, subsample_by = 60)
              out <- list(K=K,clust = clust_seg)
            }else{
              out <- list(K=K,clust = NULL)
            }
          } else {
            out <- list(K=K,clust = NULL)
          }
          msg <- paste0(sam,"| done"); message(msg)

          #output
          out
        },error=function(err) {
          msg <- paste0(sam,"|",err); message(msg)

          K <- segment.K_initial[sam]
          if (K>1) {

            shift_seg <- segmentation(testdata, lmin=300, Kmax = 10, seg.var = c("z"), subsample_by = 60, scale.variable = FALSE)
            K <- shift_seg$Kopt.lavielle
            # K can be 1 here
            if(K>1){
              clust_seg <- segclust(testdata, lmin=300, Kmax=10, ncluster = (2:K),
                                    seg.var = c("z"), scale.variable = FALSE, subsample_by = 60)
              out <- list(K=K,clust = clust_seg)
            }else{
              out <- list(K=K,clust = NULL)
            }
          } else {
            out <- list(K=K,clust = NULL)
          }
          msg <- paste0(sam,"| done"); message(msg)

          #output
          out


        }) # end of tryCatch

        return(output)
      }
    )

  clust.list <- seg_result %>% lapply(function(x) x$clust)
  segment.K <- seg_result %>% vapply(function(x) x$K,0)

  return(
    list(
      sample_Ids =refupate_data$sample_Ids,
      clust.list = clust.list,
      segment.K = segment.K
    )
  )
}





#' Run Z-score-only segmentation for samples for which 2d segmentation failed
#'
#' @param segupdated_data list containing 5 elements : sample_ids,Y,Z,segment.K, and clust.list. Y, Z corresponds to Y_recentered and Z_recentered in the return value of function update_reference_segments. segment.K and clust.list corresponds to those in the return value of function get_segments_and_clusters.
#'
#' @return list containing 5 elements. Same format as input argument segupdated_data but the segmentation results(segment.K and clust.list) are updated
#' @noRd
finalize_segments_and_clusters <- function(segupdated_data){


  rescued_data <- segupdated_data
  for (sam in which(segupdated_data$segment.K==0)) {
    message(sam)
    sampleID <- segupdated_data$sample_Ids[sam]
    testdata <- data.frame(z=segupdated_data$Z[,sam],y=segupdated_data$Y[,sam])

    shift_seg <- segmentation(testdata, lmin=300, Kmax = 10, seg.var = c("z"), subsample_by = 60, scale.variable = FALSE)
    K <- shift_seg$Kopt.lavielle
    rescued_data$segment.K[sam] <- K
    if (K>1) {
      clust_seg <- segclust(testdata, lmin=300, Kmax=10, ncluster = (2:K),
                            seg.var = c("z"), scale.variable = FALSE, subsample_by = 60)
      rescued_data$clust.list[[sam]] <- clust_seg
    } else {
      rescued_data$segment.K[sam] <- K
    }
  }

  return(rescued_data)


}



#' Title
#'
#' @param rescued_data get_segments_and_clusters output
#' @param X raw depth matrix
#' @param cores number of threads to use
#' @param prob.cutoff breakpoint detection threshold
#' @param BP.ydepth minimum change of normalized depth
#' @param X.cutoff minimum change of raw depth
#' @param q.Y.BP1 breakpoint detection criteria 1
#' @param q.Y.BP2 breakpoint detection criteria 2
#' @param half.width half of breakpoint search window size
#' @param verbose logical whether to print out information for debugging
#'
#' @return updated Y matrix
#' @noRd
detect_bp__update_ref <- function(
    rescued_data
    ,X
    ,cores=10L
    ,prob.cutoff = 0.95
    ,BP.ydepth = 0.08
    ,X.cutoff = 10
    ,q.Y.BP1 = 5
    ,q.Y.BP2 = 10
    ,half.width = 250
    ,verbose=FALSE
){
  n <- NCOL(X)
  d <- NROW(X)

  Xdiff <- apply(X, 2, FUN=function(x) diff(x))
  Xdiff_list <-
    Xdiff |>
    data.table() |>
    as.list()
  Ydiff <- apply(rescued_data$Y, 2, FUN=function(x) diff(x))
  Ydiff_list <-
    Ydiff |>
    data.table() |>
    as.list()
  Y <- rescued_data$Y
  Z <- rescued_data$Z



  if(verbose) message("    BPs2.test05")
  BPs2.test05 <-
    parallel::mclapply(mc.cores=cores,
                       X=Ydiff_list,FUN = function(y) detect_BPs(y=y, q.Y=q.Y.BP1))

  # detect_BPs
  if(verbose) message("    BPs2.test10")
  BPs2.test10 <- parallel::mclapply(mc.cores=cores,X=Ydiff_list,FUN = function(y) detect_BPs(y=y, q.Y=q.Y.BP2))
  # sequential
  #    user  system elapsed
  # 243.039   4.520 247.916

  # parallel with 10 cores
  # user  system elapsed
  # 210.077   6.695  28.099


  if(verbose) message("    BPcut1")
  BPcut1 <- lapply(seq_len(dim(X)[2]), FUN=function(sam) which(abs(Xdiff[,sam]) > X.cutoff))

  if(verbose) message("    BPcut2")
  BPcut2 <- lapply(seq_len(dim(X)[2]), FUN=function(sam) which(abs(Ydiff[,sam]) > BP.ydepth))

  if(verbose) message("    BPs2.level1")
  BPs2.level1 <- lapply(seq_len(dim(X)[2]), FUN=function(sam) base::intersect(BPs2.test05[[sam]], BPcut1[[sam]]))

  if(verbose) message("    BPs2.level2")
  BPs2.level2 <- lapply(seq_len(dim(X)[2]), FUN=function(sam) base::intersect(BPs2.test10[[sam]], BPcut1[[sam]]))

  if(verbose) message("    BPs2.level1.plus")
  BPs2.level1.plus <-
    lapply(
      seq_len(dim(X)[2]),
      FUN=function(sam) base::intersect(BPs2.level1[[sam]], BPcut2[[sam]]))
  #     user   system  elapsed
  # 1911.389    5.509 1921.437


  # initial run
  if(verbose) message("    Y_new")
  Y_new <- seq_len(ncol(X)) %>%
    mclapply(mc.cores = cores,
             \(i) {message(i);detect_bp__update_ref__sub_updateY(i,X_input=X,Y=Y,Z=Z,BPs2.level1.plus=BPs2.level1.plus,BPs2.level2=BPs2.level2)}
    ) %>%
    as.data.frame() %>% as.matrix() %>% magrittr::set_names(NULL)
  Z_new <-  t(apply(Y_new,1,function(x) pd.rate.hy(x,qrsc=TRUE)))

  # second-round run
  if(verbose) message("    Y_new2")
  Y_new2 <- seq_len(ncol(X)) %>%
    mclapply(mc.cores = cores,
             detect_bp__update_ref__sub_updateY,X_input=X,Y=Y_new,Z=Z_new,BPs2.level1.plus=BPs2.level1.plus,BPs2.level2=BPs2.level2
    ) %>%
    as.data.frame() %>% as.matrix() %>% magrittr::set_names(NULL)
  Z_new2 <-  t(apply(Y_new2,1,function(x) pd.rate.hy(x,qrsc=TRUE)))
  # user   system  elapsed
  # 1403.340   30.310  207.518

  return(
    Y_new2
  )
}



# segmenting
detect_bp__update_ref__sub_updateY <- function(sam,X_input,Y,Z,BPs2.level1.plus,BPs2.level2){

  sampleID <- colnames(X_input)[sam]
  n <- NCOL(X_input)
  d <- NROW(X_input)
  testdata <- data.frame(x=X_input[,sam], y=Y[,sam],z=Z[,sam])
  Ydiff_sam <- diff(testdata$y)

  BPs1 <- BPs2.level1.plus[[sam]]
  BPs2 <- BPs2.level2[[sam]]

  BPs <- unique(c(BPs1, BPs2))
  BPs.sam <- length(BPs)
  if (BPs.sam > 0) {
    if (BPs.sam %% 2 > 0) {
      BPs.sam <- BPs.sam + 1
      newBP <- c(seq_along(Ydiff_sam)[-BPs])[which.max(abs(Ydiff_sam[-BPs]))]
      BPs <- c(BPs,newBP)
    }


    nseg <- length(BPs)
    nseg2 <- nseg + 1
    srt.BPs <- sort(BPs)
    segment.table <- data.frame(seg = seq_len(nseg2),
                                begin = c(1,srt.BPs+1),
                                end = c(srt.BPs, d))

    segbps <- lapply(seq_len(nseg), FUN=function(i) segment.table$begin[i]:segment.table$end[i])
    #first and last segments are connected
    segbps[[1]] <- c(segbps[[1]], segment.table$begin[nseg2]:segment.table$end[nseg2])

    ## Find a baseline
    Y.tmp <- Y
    Q1.zscore <- max.zscore <- mean.zscore.outwinreg <- rep(1000,nseg)
    med.outwinreg <- len.outwinreg <- rep(0,nseg)
    names(mean.zscore.outwinreg) <- seq_len(nseg)
    for (js in seq_len(nseg)) {
      outwinreg <- segbps[[js]]
      len.outwinreg[js] <- length(outwinreg)
      newmed <- median(Y[outwinreg,sam])
      med.outwinreg[js] <- newmed
      if (newmed > 0.05) {
        Y.tmp[,sam] <- Y[,sam]/newmed
        Z.tmp <- t(apply(Y.tmp,1,function(x) pd.rate.hy(x,qrsc=TRUE)))
        mean.zscore.outwinreg[js] <- mean(abs(Z.tmp[outwinreg,sam]))
        Q1.zscore[js] <- quantile(abs(Z.tmp[,sam]),probs=0.25)
        max.zscore[js] <- max(abs(Z.tmp[,sam]))
      }
    }

    ireg <- which((max.zscore < 30) & (len.outwinreg > 500))
    if (length(ireg) > 0) {
      baseseg <- ireg[which.min(Q1.zscore[ireg])]
    } else {
      ireg <- which(max.zscore < 30)
      baseseg <- ireg[which.min(Q1.zscore[ireg])]
    }
    baseseg
    js <- baseseg
    outwinreg <- segbps[[js]]
    newmed <- median(Y[outwinreg,sam])
    return(Y[,sam]/newmed)
  }else{
    return(Y[,sam])
  }
  message(sam, "\n")

}




#' Final results
#'
#' @param X raw depth matrix
#' @param Y detect_bp__update_ref output
#' @param verbose logical whether to print out information for debugging
#'
#' @return list containing segmentation result
#' @noRd
segmentation_2nd_phase <- function(X,Y,cores=10,verbose=FALSE){
  n <- NCOL(X)
  d <- NROW(X)

  Z <- t(apply(Y,1,function(x) pd.rate.hy(x,qrsc=TRUE)))
  segtable.list  <- vector("list", length=ncol(X))
  segtable2.list <- vector("list", length=ncol(X))


  if(verbose) message("    segtable.lists")
  segtable.lists <-
    seq_len(ncol(X)) %>%
    mclapply(mc.cores=cores,function(sam){
      sampleID <- colnames(X)[sam]

      if(verbose) message("    get_BPsegment_v2")
      segment.table <- get_BPsegment_v2(X=X, Y=Y,Z=Z, sampleID=sampleID)
      if(verbose) message("    compute_CN1")
      slst <-
        data.frame(segment.table, cn=compute_CN(segment.table))

      if(verbose) message("    smooth_segment")
      segment.table2 <- smooth_segment(segment.table)

      if(verbose) message("    compute_CN2")
      slst2 <-
        data.frame(segment.table2, cn=compute_CN(segment.table2, base.state=1))


      return(list(slst=slst,slst2=slst2))
    })


  # 10 cores
  #     user   system  elapsed
  # 1770.612   31.677  268.775
  for (sam in seq_len(ncol(X))) {
    sampleID <- colnames(X)[sam]
    segtable.list[[sam]] <- data.frame(id=sam, sampleID=colnames(X)[sam], segtable.lists[[sam]]$slst)
    segtable2.list[[sam]] <- data.frame(id=sam,sampleID=colnames(X)[sam], segtable.lists[[sam]]$slst2)
  }

  return(
    list(
      segtable.list  = segtable.list,
      segtable2.list = segtable2.list

    )
  )
}








#' Run ELViS using input raw depth matrix
#'
#' @param X Raw depth matrix of position x samples
#' @param N_cores THe number of cores to use (Default : `min(10L,parallel::detectCores())`)
#' @param reduced_output logical indicating whether to return only reduced output
#' @param verbose logical whether to print out information for debugging
#'
#' @return list containing ELViS run results
#' @export
#'
#' @examples
#'
#' data(mtrx_samtools_reticulate)
#' th<-50
#' filtered_mtrx <- filt_samples(mtrx_samtools_reticulate,th=th,smry_fun=max)
#'
#' result <- run_ELViS(filtered_mtrx[,seq_len(5)],N_cores=1L)
#'
run_ELViS <-
  function(
    X
    ,N_cores=min(10L,parallel::detectCores())
    ,reduced_output=TRUE
    ,verbose=FALSE
  ){
    run_ELViS_core(
      X
      ,N_cores=N_cores
      ,save_intermediate_data = FALSE
      ,save_idx=NULL
      ,save_dir="save_dir"
      ,overwrite=FALSE
      ,reduced_output=reduced_output
      ,verbose=verbose
    )


  }


#' Run ELViS using input raw depth matrix, parameters for devel included
#'
#' @param X Raw depth matrix of position x samples
#' @param N_cores THe number of cores to use (Default : `min(10L,parallel::detectCores())`)
#' @param save_intermediate_data Logical indicating whether to save intermediate data.
#' @param save_idx save file index
#' @param save_dir directory to save intermediate files
#' @param overwrite logical whether to overwrite intermediate files rather than reading from them
#' @param reduced_output logical indicating whether to return only reduced output
#' @param verbose logical whether to print out information for debugging
#'
#' @return list containing ELViS run results
#' @noRd
run_ELViS_core <- function(
    X
    ,N_cores=min(10,parallel::detectCores())
    ,save_intermediate_data = FALSE
    ,save_idx=NULL
    ,save_dir="save_dir"
    ,overwrite=FALSE
    ,reduced_output=TRUE
    ,verbose=TRUE
){

  is_save <- save_intermediate_data
  if(is_save){
    dir.create(save_dir)
    if(is.null(save_idx)){save_idx <- 1}
  }

  message("ELViS run starts.")

  Name<-"normalized_data"
  q.cutoff<-4
  if(verbose) message(Name)
  if(is_save){
    rds_file <- glue("{save_dir}/{Name}_{save_idx}.rds")
    if(verbose) message(rds_file)
    if(overwrite|(!file.exists(rds_file))){
      normalized_data <- get_normalized_data(X,q.cutoff=q.cutoff)
      saveRDS(normalized_data,rds_file)
    }else{
      normalized_data <- readRDS(rds_file)
    }
  }else{
    normalized_data <- get_normalized_data(X,q.cutoff=q.cutoff)
  }



  Name<-"refupate_data"
  if(verbose) message(Name)
  if(is_save){
    rds_file <- glue("{save_dir}/{Name}_{save_idx}.rds")
    if(verbose) message(rds_file)
    if(overwrite|(!file.exists(rds_file))){
      refupate_data <-
        update_reference_segments(normalized_data,cores = N_cores,verbose=verbose)
      saveRDS(refupate_data,rds_file)
    }else{
      refupate_data <- readRDS(rds_file)
    }
  }else{
    refupate_data <-
      update_reference_segments(normalized_data,cores = N_cores,verbose=verbose)
  }

  Name<-"refupate_data_secondrun"
  if(verbose) message(Name)
  if(is_save){
    rds_file <- glue("{save_dir}/{Name}_{save_idx}.rds")
    if(verbose) message(rds_file)
    if(overwrite|(!file.exists(rds_file))){
      refupate_data_secondrun <-
        update_reference_segments(
          list(
            sample_Ids=refupate_data$sample_Ids,
            Y=refupate_data$Y_recentered,
            Z=refupate_data$Z_recentered
          )
          ,cores = N_cores
          ,verbose=verbose
        )
      saveRDS(refupate_data_secondrun,rds_file)
    }else{
      refupate_data_secondrun <- readRDS(rds_file)
    }
  }else{
    refupate_data_secondrun <-
      update_reference_segments(
        list(
          sample_Ids=refupate_data$sample_Ids,
          Y=refupate_data$Y_recentered,
          Z=refupate_data$Z_recentered
        )
        ,cores = N_cores
        ,verbose=verbose
      )
  }

  message("Normalization done.")


  Name<-"refupate_data_secondrun_seg_cls"
  if(verbose) message(Name)
  if(is_save){
    rds_file <- glue("{save_dir}/{Name}_{save_idx}.rds")
    if(verbose) message(rds_file)
    if(overwrite|(!file.exists(rds_file))){
      refupate_data_secondrun_seg_cls <-
        get_segments_and_clusters(refupate_data_secondrun,cores=N_cores)

      saveRDS(refupate_data_secondrun_seg_cls,rds_file)
    }else{
      refupate_data_secondrun_seg_cls <- readRDS(rds_file)
    }
  }else{
    refupate_data_secondrun_seg_cls <-
      get_segments_and_clusters(refupate_data_secondrun,cores=N_cores)
  }



  Name<-"rescued_data_p2"
  if(verbose) message(Name)
  if(is_save){
    rds_file <- glue("{save_dir}/{Name}_{save_idx}.rds")
    if(verbose) message(rds_file)
    if(overwrite|(!file.exists(rds_file))){
      rescued_data_p2 <- segupdated_data_p2 <-
        list(
          sample_Ids = refupate_data_secondrun$sample_Ids,
          Y =          refupate_data_secondrun$Y_recentered,
          Z =          refupate_data_secondrun$Z_recentered,
          segment.K =  refupate_data_secondrun_seg_cls$segment.K,
          clust.list = refupate_data_secondrun_seg_cls$clust.list
        )

      saveRDS(rescued_data_p2,rds_file)
    }else{
      rescued_data_p2 <- readRDS(rds_file)
    }
  }else{
    rescued_data_p2 <- segupdated_data_p2 <-
      list(
        sample_Ids = refupate_data_secondrun$sample_Ids,
        Y =          refupate_data_secondrun$Y_recentered,
        Z =          refupate_data_secondrun$Z_recentered,
        segment.K =  refupate_data_secondrun_seg_cls$segment.K,
        clust.list = refupate_data_secondrun_seg_cls$clust.list
      )
  }



  Name<-"new_Y_p2"
  if(verbose) message(Name)
  if(is_save){
    rds_file <- glue("{save_dir}/{Name}_{save_idx}.rds")
    if(verbose) message(rds_file)
    if(overwrite|(!file.exists(rds_file))){
      new_Y_p2 <- detect_bp__update_ref(
        rescued_data = rescued_data_p2,
        X,
        cores=N_cores,
        prob.cutoff = 0.95,
        BP.ydepth = 0.08,
        X.cutoff = 10,
        q.Y.BP1 = 5,
        q.Y.BP2 = 10,
        half.width = 250,
        verbose=verbose
      )
      saveRDS(new_Y_p2,rds_file)
    }else{
      new_Y_p2 <- readRDS(rds_file)
    }
  }else{
    new_Y_p2 <- detect_bp__update_ref(
      rescued_data = rescued_data_p2,
      X,
      cores=N_cores,
      prob.cutoff = 0.95,
      BP.ydepth = 0.08,
      X.cutoff = 10,
      q.Y.BP1 = 5,
      q.Y.BP2 = 10,
      half.width = 250,
      verbose=verbose
    )
  }


  Name<-"segtable_p2"
  if(verbose) message(Name)
  if(is_save){
    rds_file <- glue("{save_dir}/{Name}_{save_idx}.rds")
    if(verbose) message(rds_file)
    if(overwrite|(!file.exists(rds_file))){
      segtable_p2 <- segmentation_2nd_phase(X,Y=new_Y_p2,cores=N_cores,verbose=verbose)
      saveRDS(segtable_p2,rds_file)
    }else{
      segtable_p2 <- readRDS(rds_file)
    }
  }else{
    segtable_p2 <- segmentation_2nd_phase(X,Y=new_Y_p2,cores=N_cores,verbose=verbose)
  }

  message("Segmentation done.")


  final_output <- rbindlist(segtable_p2$segtable2.list )
  final_call_tmp <-
    final_output %>%
    dplyr::group_by(
      .data$sampleID
    ) %>%
    dplyr::summarise(N_seg = n(),N_state = length(unique(.data$state)))
  final_call <-
    list(segmented_samples = which(final_call_tmp$N_seg!=1), cnv_samples = which(final_call_tmp$N_state!=1))
  if(reduced_output){
    return(
      list(
        is_reduced_output = TRUE,
        final_output = final_output,
        final_call = final_call,
        new_Y_p2 = new_Y_p2 %>% magrittr::set_colnames(colnames(X))
      )
    )
  }else{
    return(
      list(
        is_reduced_output = FALSE,
        final_output = final_output,
        final_call = final_call,
        normalized_data = normalized_data,
        refupate_data = refupate_data,
        refupate_data_secondrun = refupate_data_secondrun,
        refupate_data_secondrun_seg_cls = refupate_data_secondrun_seg_cls,
        rescued_data_p2 = rescued_data_p2,
        new_Y_p2 = new_Y_p2,
        segtable_p2 = segtable_p2
      )

    )
  }

}











#' Get new baselines according to criteria user designates
#'
#' @param result Run result
#' @param mode Indecate how new baseline should be set ("longest","shortest")
#'
#' @return a integer vector indicating new baseline index for each sample
#' @export
#' @rawNamespace import(dplyr, except=c(as_data_frame,between,collapse,desc,first,groups,intersect,last,setdiff,union))
#'
#' @examples
#'
#' # its usage example is given in vignette in detail
#'
#' data(ELViS_toy_run_result)
#' result <- ELViS_toy_run_result
#'
#' get_new_baseline(result,mode="longest")
#'
get_new_baseline <- function(result,mode="longest"){

  tmp_data <-
    result$final_output %>%
    dplyr::group_by(id,.data$state) %>%
    dplyr::summarise(total_length = sum(.data$end-.data$begin+1)) %>%
    dplyr::ungroup()


  if(mode == "longest"){

    new_baseline <-
      tmp_data %>%
      dplyr::group_by(id) %>%
      dplyr::slice_max(.data$total_length,n = 1) %>%
      {.$state}

  }else if(mode == "shortest"){

    new_baseline <-
      tmp_data %>%
      dplyr::group_by(id) %>%
      dplyr::slice_min(.data$total_length,n = 1) %>%
      {.$state}

  }else{
    stop(glue("{mode} is not an available mode."))
  }

  return(new_baseline)

}




anno_color_var <- function(variable, palette=NULL, palette_id=NULL, palette_order=NULL,
                           paired=FALSE, cont_midval=NULL) {
  variable0 <- variable
  if (sum(class(variable) %in% c("double","numeric")) > 0) {
    if (is.null(palette) & is.null(palette_id)) palette_id <- 1
    if (is.null(palette)) palette <- names(choi_continuous_palettes)[palette_id]
    colors <- choi_continuous_palettes[[eval(palette)]]
    if (! is.null(cont_midval)) {
      var_col <- colorRamp2(c(min(variable), cont_midval, max(variable)),
                            colors)
    } else {
      var_col <- colorRamp2(quantile(variable,
                                     probs=seq(0,1,length=length(colors))),
                            colors)
    }
    return(list(var=variable0, var_col=var_col))
  } else {
    if (sum(class(variable) %in% "factor") == 0) variable <- factor(variable)
    if (length(levels(variable))>2) paired <- FALSE
    if (paired) {
      if (is.null(palette) & is.null(palette_id)) palette_id <- 1
      if (is.null(palette)) palette <- names(choi_paired_palettes)[palette_id]
      colors <- choi_paired_palettes[[eval(palette)]]
      var_col <- colors[seq_len(length(levels(variable)))]
      names(var_col) <- levels(variable)
    } else {
      if (sum(class(variable) %in% c("ordered")) > 0) {
        if (is.null(palette) & is.null(palette_id)) palette_id <- 3
        if (is.null(palette)) palette <- names(choi_ordered_palettes)[palette_id]
        colors <- choi_ordered_palettes[[eval(palette)]]
        if (! is.null(palette_order)) {
          colors <- colors[palette_order]
        }
        var_col <- colors[seq_len(length(levels(variable)))]
        names(var_col) <- levels(variable)
      } else {
        if (is.null(palette) & is.null(palette_id)) palette_id <- 19
        if (is.null(palette)) palette <- names(choi_discrete_palettes)[palette_id]
        colors <- choi_discrete_palettes[[eval(palette)]]
        if (! is.null(palette_order)) {
          colors <- colors[palette_order]
        }
        var_col <- colors[seq_len(length(levels(variable)))]
        names(var_col) <- levels(variable)
      }
    }
    return(list(var=variable0, var_col=var_col,
                col=as.character(factor(variable0, levels=names(var_col), labels=var_col))))
  }
}



#' @noRd
detect_BPs <- function(y, half.width=250, prob.cutoff=0.95, q.Y=5) {
  rollYquant0 <- rollapply(abs(y),width=2*half.width,FUN=function(x) quantile(x,probs=prob.cutoff))
  rollYquant <- c(rep(rollYquant0[1],half.width), rollYquant0, rep(tail(rollYquant0,1),half.width-1))

  ## remove consecutive BPs
  Yexceed <- abs(y)/abs(q.Y*rollYquant)
  BPcandi <- which(Yexceed > 1)
  # breakpoint candidates : absolute difference larger than 5 times of Q95 among +- half.width positions

  if (length(BPcandi) > 2) {
    BPcandi.info <- data.frame(BPs=BPcandi[2:length(BPcandi)],
                               BPs.sign=y[seq_len(length(BPcandi)-1)]*y[2:length(BPcandi)],
                               BPs.diff=diff(BPcandi))
    BPcandi.info$rm.BPs <- rep(0,dim(BPcandi.info)[1])
    BPcandi.info$rm.BPs[which((BPcandi.info[,3] < 5) & (BPcandi.info[,2] > 0))] <- 1
    BPs2 <- BPcandi[! BPcandi %in% BPcandi.info$BPs[which(BPcandi.info$rm.BPs==1)]]
  } else {
    BPs2 <- BPcandi
  }
  return(BPs2)
}

#' calculate pooled standard deviation
#' @noRd
pooled_sd <- function(s1, s2, n1, n2) {
  sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
}


#' Get breakpoints
#'
#' @param X raw depth
#' @param Y normalized depth
#' @param Z standardized depth
#' @param sampleID sample ID
#' @param q.Y.BP1 breakpoint detection criteria 1
#' @param q.Y.BP2 breakpoint detection criteria 2
#' @param BP.ydepth minimum change of normalized depth
#' @param BP.xdepth minimum change of raw depth
#'
#' @return data.frame containing segmentation results
#' @rawNamespace import(igraph, except=c(as_data_frame,degree,groups,union))
#' @import stringr
#' @noRd
get_BPsegment_v2 <- function(X, Y, Z=NULL, sampleID,
                             q.Y.BP1 = 5,
                             q.Y.BP2 = 10,
                             BP.ydepth = 0.08, BP.xdepth = 10) {
  if (is.numeric(sampleID)) {
    sam <- sampleID
  } else {
    sam <- which(colnames(X)==sampleID)
  }
  if (is.null(Z)) {
    Z <- t(apply(Y,1,function(x) pd.rate.hy(x,qrsc=TRUE)))
  }

  d <- NROW(X)

  xdiff <- diff(X[,sam])
  ydiff <- diff(Y[,sam])
  ymed <- apply(Y, 1, median)

  testdata <- data.frame(x=X[,sam], y=Y[,sam], z=Z[,sam])

  BPs.q10 <- detect_BPs(y=ydiff, q.Y=q.Y.BP2, half.width=250, prob.cutoff = 0.95)
  BPs.q05 <- detect_BPs(y=ydiff, q.Y=q.Y.BP1, half.width=250, prob.cutoff = 0.95)
  BPs.cutx <- which(abs(xdiff) > BP.xdepth)
  BPs.cuty <- which(abs(ydiff) > BP.ydepth)


  BPs.level1 <- base::intersect(BPs.q05, BPs.cutx)
  BPs.level2 <- base::intersect(BPs.q10, BPs.cutx)
  BPs.level1.plus <- base::intersect(BPs.level1, BPs.cuty)

  BPs <- unique(c(BPs.level1.plus, BPs.level2))

  BPs.sam <- length(BPs)
  if (BPs.sam > 0) {
    if (BPs.sam %% 2 > 0) {
      BPs.sam <- BPs.sam + 1
      newBP <- c(c(seq_len(length(ydiff)))[-BPs])[which.max(abs(ydiff[-BPs]))]
      BPs <- c(BPs,newBP)

    }
    BPs <- sort(BPs)

    nseg <- length(BPs)
    nseg2 <- nseg + 1
    srt.BPs <- sort(BPs)
    segment.table <- data.frame(seg = seq_len(nseg2),
                                begin = c(1,srt.BPs+1),
                                end = c(srt.BPs, d))
    segbps <- lapply(seq_len(nseg), FUN=function(i) segment.table$begin[i]:segment.table$end[i])
    segbps[[1]] <- c(segbps[[1]], segment.table$begin[nseg2]:segment.table$end[nseg2])

    ## Find a baseline
    Y.tmp <- Y
    Q1.zscore <- max.zscore <- mean.zscore.outwinreg <- rep(1000,nseg)
    med.outwinreg <- len.outwinreg <- rep(0,nseg)
    names(mean.zscore.outwinreg) <- seq_len(nseg)
    for (js in seq_len(nseg)) {
      outwinreg <- segbps[[js]]
      len.outwinreg[js] <- length(outwinreg)
      newmed <- median(Y[outwinreg,sam])
      med.outwinreg[js] <- newmed
      if (newmed > 0.05) {
        Y.tmp[,sam] <- Y[,sam]/newmed
        Z.tmp <- t(apply(Y.tmp,1,function(x) pd.rate.hy(x,qrsc=TRUE)))
        mean.zscore.outwinreg[js] <- mean(abs(Z.tmp[outwinreg,sam]))
        Q1.zscore[js] <- quantile(abs(Z.tmp[,sam]),probs=0.25)
        max.zscore[js] <- max(abs(Z.tmp[,sam]))
      }
    }
    ireg <- which((max.zscore < 30) & (len.outwinreg > 500))
    if (length(ireg) > 0) {
      baseseg <- ireg[which.min(Q1.zscore[ireg])]
    } else {
      ireg <- which(max.zscore < 30)
      baseseg <- ireg[which.min(Q1.zscore[ireg])]
    }

    ## Make a table for each segment
    seg_len_vec <- vapply(segbps,FUN=length,0)
    segtable <- segment.table[seq_len(nseg), ] # for now, drop the last segment (this is the same as the first segment)
    segtable$mu.x <- vapply(segbps, FUN=function(ibps) mean(testdata$x[ibps]),0)
    segtable$sd.x <- vapply(segbps, FUN=function(ibps) sd(testdata$x[ibps]),0)
    segtable$mu.y <- vapply(segbps, FUN=function(ibps) mean(testdata$y[ibps]),0)
    segtable$sd.y <- vapply(segbps, FUN=function(ibps) sd(testdata$y[ibps]),0)
    segtable$mu.z <- vapply(segbps, FUN=function(ibps) mean(testdata$z[ibps]),0)
    segtable$sd.z <- vapply(segbps, FUN=function(ibps) sd(testdata$z[ibps]),0)
    segtable$cohort.y <- vapply(segbps, FUN=function(ibps) mean(ymed[ibps]),0)

    segtable$baseline <- rep(0,nseg)
    segtable$baseline[baseseg] <- 1

    ## Find more base segments
    baseseg2 <- which(abs(segtable$mu.y - segtable$mu.y[baseseg]) < segtable$sd.y[baseseg])
    baseseg3 <- which(abs(segtable$mu.z - segtable$mu.z[baseseg]) < segtable$sd.z[baseseg])
    segtable$baseline[base::intersect(baseseg2,baseseg3)] <- 1
    segtable$state <- segtable$baseline

    i_nonbase <- which(segtable$baseline==0)
    n_nonbase <- length(i_nonbase)
    if(n_nonbase!=0){
      if(n_nonbase==1){
        segtable$state[which(segtable$baseline==0)] <-
          seq_len(length(which(segtable$baseline==0))) + 1
      }else if(n_nonbase>=2){

        is_above <- segtable$mu.y[i_nonbase] > segtable$mu.y[baseseg]
        i_nonbase_sp <-
          split(
            i_nonbase,is_above,drop = TRUE)

        # find close copynumber segment groups
        i_grouped_lst <-
          i_nonbase_sp |>
          lapply(\(i_tmp){

            if(length(i_tmp)==1){
              out <- list(i_tmp)
            }else if(length(i_tmp)==2){

              psd_y <- pooled_sd(segtable$sd.y[i_tmp[1]],segtable$sd.y[i_tmp[2]],seg_len_vec[i_tmp[1]],seg_len_vec[i_tmp[2]])
              is_simil_y <- abs(diff(segtable$mu.y[i_tmp])) < psd_y
              psd_z <- pooled_sd(segtable$sd.z[i_tmp[1]],segtable$sd.z[i_tmp[2]],seg_len_vec[i_tmp[1]],seg_len_vec[i_tmp[2]])
              is_simil_z <- abs(diff(segtable$mu.z[i_tmp])) < psd_z

              if(is_simil_y&is_simil_z){
                out <- list(i_tmp)  # same seg
              }else{
                out <- as.list(i_tmp) # different seg
              }
            }else if(length(i_tmp)>=3){
              i_tmp_cmb <- utils::combn(i_tmp,2)
              max_rel_d <- i_tmp_cmb |>
                apply(2,\(j_tmp){
                  psd_y <- pooled_sd(segtable$sd.y[j_tmp[1]],segtable$sd.y[j_tmp[2]],seg_len_vec[j_tmp[1]],seg_len_vec[j_tmp[2]])
                  simil_y <- abs(diff(segtable$mu.y[j_tmp]))/psd_y
                  psd_z <- pooled_sd(segtable$sd.z[j_tmp[1]],segtable$sd.z[j_tmp[2]],seg_len_vec[j_tmp[1]],seg_len_vec[j_tmp[2]])
                  simil_z <- abs(diff(segtable$mu.z[j_tmp]))/psd_z
                  return(c(max(simil_y,simil_z) ))
                })

              g <- gen_dist_graph(i_tmp,i_tmp_cmb,max_rel_d,th_d=1)


              subn_grps <-
                get_subcl(g) |>
                lapply(str_replace,"^v","") |>
                lapply(as.numeric)

              out <- subn_grps

            }

            return(out)
          }) |> unlist(recursive = FALSE)

        # update group
        for( k_tmp in  i_grouped_lst[order(vapply(i_grouped_lst,min,0))]){
          segtable$state[k_tmp] <- max(segtable$state) + 1
        }

      }
    }


    segtable2 <- rbind(segtable,segtable[1,])
    segtable2[nseg2,seq_len(3)] <- segment.table[nseg2,seq_len(3)]
    return(segtable2)
  } else {
    segtable <- data.frame(seg = 1, begin = 1, end = dim(testdata)[1],
                           mu.x = mean(testdata$x), sd.x = sd(testdata$x),
                           mu.y = mean(testdata$y), sd.y = sd(testdata$y),
                           mu.z = mean(testdata$z), sd.z = sd(testdata$z),
                           cohort.y = mean(ymed),
                           baseline = 1, state = 1)
    return(segtable)
  }
}


gen_dist_graph <- function(i_tmp,i_tmp_cmb,max_rel_d,th_d){
  max_rel_s <- 1/(max_rel_d+1)

  is_close <- max_rel_d < th_d


  g <- make_empty_graph(n = length(i_tmp),directed = FALSE)
  V(g)$name <- paste0("v",i_tmp)

  if(any(is_close)){
    edges <- paste0("v",i_tmp_cmb[,is_close] |> as.vector())
    g <- add_edges(g, edges)  # Connect node 1 to 2 and node 3 to 4
    E(g)$weight <- max_rel_s[is_close]
  }


  return(g)
}


clique_weight <- function(clique, graph) {
  subg <- induced_subgraph(graph, vids = clique)
  return(sum(E(subg)$weight))
}

get_subcl <- function(g){
  cliques_list <- cliques(g)

  clique_weights <- vapply(X = cliques_list, FUN = clique_weight,FUN.VALUE = numeric(0), graph = g)

  ordered_cliques <- cliques_list[order(-clique_weights)]

  subnetworks <- list()
  assigned_nodes <- c()

  for (clique in ordered_cliques) {
    # Check if any nodes in the clique are already assigned
    if (all(!(clique %in% assigned_nodes))) {
      subnetworks <- append(subnetworks, list(clique))
      assigned_nodes <- c(assigned_nodes, clique)
    }
  }

  out_subnetworks <- subnetworks |> lapply(\(x) induced_subgraph(g,vids=x)) |> lapply(V) %>% lapply(names)
  return(out_subnetworks)
}




get_window2 <- function(y,cutoff=1,min.length=10) {
  d <- length(y)
  ips <- which((y-cutoff)[seq_len(d-1)]*(y-cutoff)[2:d]<=0)
  win <- matrix(c(1,rep(ips,each=2),d),ncol=2,byrow=TRUE)
  win <- win[which((win[,2] - win[,1])>min.length),]
  return(win)
}

get_window <- function(y,pos=NULL,cutoff=1,min.length=10) {
  if (is.null(pos)) {
    return(get_window2(y=y,cutoff=cutoff,min.length=min.length))
  } else {
    wincan <- get_window2(y=y,cutoff=cutoff,min.length=min.length)
    windiff <- wincan - pos
    win.tmp <- wincan[which((windiff[,1]*windiff[,2])<=0),]
    return(win.tmp)
  }
}



#' Title
#'
#' @param color  color name
#' @param percent % transparency
#' @param name  an optional name for the color
#'
#' @return transparency-adjusted color
#' @noRd
t_col <- function(color, percent = 50, name = NULL) {

  ## Get RGB values for named color
  rgb.val <- col2rgb(color)

  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               maxColorValue = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)

  ## Save the color
  invisible(t.col)
}

compute_CN <- function(segment.table, base.state=1) {
  if (is.null(base.state)) base.state <- 1
  base.idx <- which(segment.table$state==base.state)
  w <- (segment.table$end - segment.table$begin)[base.idx]
  w <- w / sum(w)
  cn <- segment.table$mu.y/segment.table$cohort.y
  base.cn <- sum(w*segment.table$mu.y[base.idx])/sum(w*segment.table$cohort.y[base.idx])
  return(cn/base.cn)
}


smooth_segment <- function(segment.table) {
  output <- segment.table
  state <- unique(output$state)
  for (js in state) {
    idx <- which(segment.table$state==js)
    w <- (segment.table$end - segment.table$begin)[idx]
    w <- w / sum(w)
    output$mu.x[idx] <- sum(w*segment.table$mu.x[idx])
    output$mu.y[idx] <- sum(w*segment.table$mu.y[idx])
    output$mu.z[idx] <- sum(w*segment.table$mu.z[idx])
    output$sd.x[idx] <- sum(w*segment.table$sd.x[idx]) # this is not correct, just for convenience
    output$sd.y[idx] <- sum(w*segment.table$sd.y[idx]) # this is not correct, just for convenience
    output$sd.z[idx] <- sum(w*segment.table$sd.z[idx]) # this is not correct, just for convenience
    output$cohort.y[idx] <- sum(w*segment.table$cohort.y[idx])
  }
  return(output)
} # end of function
