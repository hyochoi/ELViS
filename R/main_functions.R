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
    # print((dot_info))
    if(length(dot_info)==0){
      title_txt <- (glue("{fun_name} Depth Distribution"))
    }else{
      title_txt <- (glue("{fun_name} {dot_info} Depth Distribution"))
    }

  }
  # print(title_txt)

  #suppressWarnings({
    ggplot(
      data.frame(Max_Depth = apply(mtrx,2,smry_fun,...))
    ) +
      aes(x = .data$Max_Depth) +
      geom_histogram(col="gray30",fill="gray80",) +
      geom_rect(aes(xmin = 0, xmax = th, ymin = 0, ymax = Inf),
                fill = alpha("gray",0.01)) +
      geom_vline(xintercept = th,col="red",linetype="dashed") +
      # scale_x_continuous(trans=log10p1_trans) +
      scale_x_continuous(trans="log10") +
      ggtitle(title_txt)
  #})

}


capture_params_glue <- function(...) {
  # Capture the names and values of the arguments
  args <- as.list(match.call())[-1]  # Exclude the function name from the list
  # print("capture_params_glue")
  # print(args)

  # Create a formatted string using glue and paste
  formatted_params <- paste0(names(args), " = ", sapply(args, deparse), collapse = " ")

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
  # breaks = log_breaks(base = 10)
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
    # return(x/(2*median(x)))
    return(x/m)
  } else {
    return(rep(0,length(x)))
  }
}


# get_window = function(Y,sam) {
#   ips = which((Y[,sam]-1)[seq_len(nrow(Y)-1)]*(Y[,sam]-1)[2:nrow(Y)]<=0)
#   win = matrix(c(1,rep(ips,each=2),dim(Y)[1]),ncol=2,byrow=TRUE)
#   win = win[which((win[,2] - win[,1])>10),]
#   return(win)
# }


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
  # ips = which((Y[,sam]-1)[seq_len(nrow(Y)-1)]*(Y[,sam]-1)[2:nrow(Y)]<=0)
  # win = matrix(c(1,rep(ips,each=2),dim(Y)[1]),ncol=2,byrow=TRUE)
  # win = win[which((win[,2] - win[,1])>10),]
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
  # ips = which((Y[,sam]-1)[seq_len(nrow(Y)-1)]*(Y[,sam]-1)[2:nrow(Y)]<=0)
  # win = matrix(c(1,rep(ips,each=2),dim(Y)[1]),ncol=2,byrow=TRUE)
  # win = win[which((win[,2] - win[,1])>10),]
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
  window.list <- sapply(c(seq_len(ncol(Y))),function(j) get_window_v1(Y=Y,sam=j)) # get windows
  # get average per window and sample, standardize across samples to get relative position in the cohort
  outmat.list <- sapply(window.list,function(j) get_outmat_v1(Y=Y,win=j)) # get per-window mean for each sample and get relative score of it, upper and lower part SDs are different

  # # max absolute deviation by sample
  # outstat is not used so do not need to calculate this
  # outstat = diag(
  #   # max deviation for each sample - no need to calculate other samples... but it's not that slow
  #   sapply(outmat.list,function(x) apply(abs(x),2,max))
  # )


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
  # window.list2 = sapply(c(seq_len(ncol(Y2))),function(j) get_window_v1(Y=Y2,sam=j))
  # # outmat.list2 = sapply(c(seq_len(ncol(Y2))),function(j) get_outmat(Y=Y2,sam=j))
  # outmat.list2 = sapply(window.list2,function(j) get_outmat_v1(Y=Y,win=j))
  # # maximal deviation
  # outstat2 = diag(sapply(outmat.list2,function(x) apply(abs(x),2,max)))
  # # minimum (max deviation in negative direction)
  # outstat.neg2 = diag(sapply(outmat.list2,function(x) apply(x,2,min)))
}




#' Updating reference segments
#'
#' @param normalized_data list containing 1 character vector and 2 numeric matrices. The eirst element is sample IDs, the second element contain normalized depth and the third element contain standardized depth.
#' @param cores integer; number of threads to use for the analysis. (Default : 10L)
#'
#' @return a list containing updated data
#' @import parallel
#' @import segclust2d
#' @rawNamespace import(circlize, except=c(degree))
#' @rawNamespace import(zoo, except=c(index,yearmon,yearqtr,"index<-"))
#' @noRd
update_reference_segments <- function(normalized_data,cores=10L){

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
  #system.time({
  # message("    segment.K_initial")
  segment.K_initial <-
    # (1:4) %>%
    seq_len(n) %>%
    mclapply(
      mc.cores<-cores,
      # lapply(

      function(sam) {
        message(sam)
        # for (sam in 75:79) {
        sampleID <- sample_id_vec[sam]
        # result =
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
          # K + 100000
          K

        }) # end of tryCatch

        return(output)
      }
    ) %>%
    unlist
  #})



  # 2 recentering data, setting least deviating segment as reference level
  # time taken 10 cores : 21.765 sec
  #    user  system elapsed
  # 120.368   6.716  21.765
  # system.time({
  # for (sam in 1:n)
  # seg_result =
  # message("    Y_recentered")
  Y_recentered <-
    # (1:4) %>%
    seq_len(n) %>%
    mclapply(
      mc.cores=cores,
      # lapply(

      function(sam) {
        # for (sam in 75:79) {
        sampleID <- colnames(Y)[sam]
        # result <-
        testdata <- data.frame(z=Z[,sam],y=Y[,sam])
        output <- tryCatch({


          K <- segment.K_initial[sam]
          if (K>1) {
            clust_seg <- segclust(testdata, lmin=300, Kmax=10, ncluster = (2:K),
                                  seg.var = c("z","y"), scale.variable = FALSE, subsample_by = 60)
            # plot(clust_seg)
            result <- segment(clust_seg)

            # Update mean
            # Setting the segment that shows the least deviation from median(=1, since the data is median-normalized)
            # as new reference segment
            mean.update <- result$mu.y[which.min(abs(result$mu.y-1))] ##########
            new_y <- testdata$y/mean.update
          } else {
            new_y <- testdata$y
          }
          # plot(clust_seg)
          # result = segment(clust_seg)
          msg <- paste0(sam,"| done"); message(msg)
          # rm(testdata,shift_seg,clust_seg,K)

          #output
          new_y
        },error=function(err) {
          msg <- paste0(sam,"|",err); message(msg)
          message("Running 1D segmentation with robust-scaled data only")


          K <- segment.K_initial[sam]
          if (K>1) {
            clust_seg <- segclust(testdata, lmin=300, Kmax=10, ncluster = (2:K),
                                  seg.var = c("z"), scale.variable = FALSE, subsample_by = 60)
            # plot(clust_seg)
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
            mean.update <- mu.y_nonzero[which.min(abs(mu.y_nonzero-1))] ##########
            new_y <- testdata$y/mean.update
          } else {
            new_y <- testdata$y
          }
          # plot(clust_seg)
          # result = segment(clust_seg)
          msg <- paste0(sam,"| done"); message(msg)
          # rm(testdata,shift_seg,clust_seg,K)

          #output
          new_y

        }) # end of tryCatch

        return(output)
      }
    ) %>%
    simplify2array

  Z_recentered <- t(apply(Y_recentered,1,function(x) pd.rate.hy(x,qrsc=TRUE)))
  # })

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
        # for (sam in 75:79) {
        sampleID <- sample_id_vec[sam]
        # result <-
        testdata <- data.frame(z=Z_recentered[,sam],y=Y_recentered[,sam])
        output <- tryCatch({

          K <- segment.K_initial[sam]
          if (K>1) {

            shift_seg <- segmentation(testdata, lmin=300, Kmax = 10, seg.var = c("z","y"), subsample_by = 60, scale.variable = FALSE)
            K <- shift_seg$Kopt.lavielle
            # segment.K[sam] = K
            # K can be 1 here
            if(K>1){
              clust_seg <- segclust(testdata, lmin=300, Kmax=10, ncluster = (2:K),
                                    seg.var = c("z","y"), scale.variable = FALSE, subsample_by = 60)
              # clust.list[[sam]] = clust_seg
              out <- list(K=K,clust = clust_seg)
            }else{
              out <- list(K=K,clust = NULL)
            }
          } else {
            # segment.K[sam] = K
            out <- list(K=K,clust = NULL)
          }
          # plot(clust_seg)
          # result = segment(clust_seg)
          msg <- paste0(sam,"| done"); message(msg)
          # rm(testdata,shift_seg,clust_seg,K)

          #output
          out
        },error=function(err) {
          msg <- paste0(sam,"|",err); message(msg)

          K <- segment.K_initial[sam]
          if (K>1) {

            shift_seg <- segmentation(testdata, lmin=300, Kmax = 10, seg.var = c("z"), subsample_by = 60, scale.variable = FALSE)
            K <- shift_seg$Kopt.lavielle
            # segment.K[sam] = K
            # K can be 1 here
            if(K>1){
              clust_seg <- segclust(testdata, lmin=300, Kmax=10, ncluster = (2:K),
                                    seg.var = c("z"), scale.variable = FALSE, subsample_by = 60)
              # clust.list[[sam]] = clust_seg
              out <- list(K=K,clust = clust_seg)
            }else{
              out <- list(K=K,clust = NULL)
            }
          } else {
            # segment.K[sam] = K
            out <- list(K=K,clust = NULL)
          }
          # plot(clust_seg)
          # result = segment(clust_seg)
          msg <- paste0(sam,"| done"); message(msg)
          # rm(testdata,shift_seg,clust_seg,K)

          #output
          out


        }) # end of tryCatch

        return(output)
      }
    )

  clust.list <- seg_result %>% lapply(function(x) x$clust)
  segment.K <- seg_result %>% sapply(function(x) x$K)

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
  # segupdated_data %>%
  # within({
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
    # rm(testdata,sampleID,K,shift_seg,clust_seg)
  }
  # rm(sam)

  # })

  return(rescued_data)
  # all(rescued_data$Z == segupdated_data$Z)
  # rescued_data$segment.K[rescued_data$segment.K != segupdated_data$segment.K]


}

#
#
# Ydiff = Ydiff_s2
# Y = rescued_data_s2$Y
# Z = rescued_data_s2$Z


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
#'
#' @return updated Y matrix
#' @noRd
detect_bp__update_ref <- function(
    rescued_data,
    X,
    cores=10L,
    prob.cutoff = 0.95,
    BP.ydepth = 0.08,
    X.cutoff = 10,
    q.Y.BP1 = 5,
    q.Y.BP2 = 10,
    half.width = 250
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



  # time5 = system.time({
  # BPs2.test05_old <- apply(Ydiff,2,FUN=function(y) detect_BPs(y=y, q.Y=q.Y.BP1))
  # message("BPs2.test05")
  BPs2.test05 <-
    parallel::mclapply(mc.cores=cores,
                       # lapply(
                       X=Ydiff_list,FUN = function(y) detect_BPs(y=y, q.Y=q.Y.BP1))
  # })
  # print(time5)
  # sequential
  #    user  system elapsed
  # 240.711   3.543 244.635

  # dim(Ydiff)
  # detect_BPs
  # time10 = system.time({
  # BPs2.test10_old <- apply(Ydiff,2,FUN=function(y) detect_BPs(y=y, q.Y=q.Y.BP2))
  # message("BPs2.test10")
  BPs2.test10 <- parallel::mclapply(mc.cores=cores,X=Ydiff_list,FUN = function(y) detect_BPs(y=y, q.Y=q.Y.BP2))
  # })
  # print(time10)
  #    user  system elapsed
  # 243.039   4.520 247.916

  # parallel with 10 cores
  # user  system elapsed
  # 210.077   6.695  28.099


  # time1 = system.time({
  # message("BPcut1")
  BPcut1 <- sapply(seq_len(dim(X)[2]), FUN=function(sam) which(abs(Xdiff[,sam]) > X.cutoff))
  # })
  # print(time1)
  # time2 = system.time({
  # message("BPcut2")
  BPcut2 <- sapply(seq_len(dim(X)[2]), FUN=function(sam) which(abs(Ydiff[,sam]) > BP.ydepth))
  # })
  # print(time2)
  # time21 = system.time({
  # message("BPs2.level1")
  BPs2.level1 <- sapply(seq_len(dim(X)[2]), FUN=function(sam) base::intersect(BPs2.test05[[sam]], BPcut1[[sam]]))
  # })
  # print(time21)
  # time22 = system.time({
  # message("BPs2.level2")
  BPs2.level2 <- sapply(seq_len(dim(X)[2]), FUN=function(sam) base::intersect(BPs2.test10[[sam]], BPcut1[[sam]]))
  # })
  # print(time22)
  # time21p = system.time({
  # message("BPs2.level1.plus")
  BPs2.level1.plus <- sapply(seq_len(dim(X)[2]), FUN=function(sam) base::intersect(BPs2.level1[[sam]], BPcut2[[sam]]))
  # })
  # print(time21p)

  #     user   system  elapsed
  # 1911.389    5.509 1921.437


  #start from here
  # save.image("assets/session_20240503_p1doing.Rdata")
  # load("assets/session_20240503_p1doing.Rdata")

  # find breakpoints and update
  # input
  #  - sample_Ids
  #  - X,Y,Z
  #  - BPs2.level1.plus
  #  - BPs2.level2
  # output
  # - Y,Z get updated
  # 1:200 %>% lapply(function(x,y,z) print(x+y+z),y=1,z=2)
  # detect_bp__update_ref__sub_updateY(1,X=X,Y=Y,Z=Z,BPs2.level1.plus=BPs2.level1.plus,BPs2.level2=BPs2.level2)

  # single run
  # time_test =
  # system.time({
  # message("Y_new")
  Y_new <- seq_len(ncol(X)) %>%
    # lapply(function(x) print(x))
    # lapply(
    mclapply(mc.cores = cores,
             # detect_bp__update_ref__sub_updateY,X_input=X,Y=Y,Z=Z,BPs2.level1.plus=BPs2.level1.plus,BPs2.level2=BPs2.level2
             \(i) {message(i);detect_bp__update_ref__sub_updateY(i,X_input=X,Y=Y,Z=Z,BPs2.level1.plus=BPs2.level1.plus,BPs2.level2=BPs2.level2)}
    ) %>%
    as.data.frame() %>% as.matrix() %>% magrittr::set_names(NULL)
  # })
  Z_new <-  t(apply(Y_new,1,function(x) pd.rate.hy(x,qrsc=TRUE))) ###############################################
  # Y_new = test %>% as.data.frame() %>% as.matrix() %>% magrittr::set_names(NULL)
  # print(time_test)

  # re run
  # time_test =
  # system.time({
  # message("Y_new2")
  Y_new2 <- seq_len(ncol(X)) %>%
    # lapply(function(x) print(x))
    # lapply(
    mclapply(mc.cores = cores,
             detect_bp__update_ref__sub_updateY,X_input=X,Y=Y_new,Z=Z_new,BPs2.level1.plus=BPs2.level1.plus,BPs2.level2=BPs2.level2
    ) %>%
    as.data.frame() %>% as.matrix() %>% magrittr::set_names(NULL)
  # })
  Z_new2 <-  t(apply(Y_new2,1,function(x) pd.rate.hy(x,qrsc=TRUE))) ###############################################
  # print(time_test)
  # Y_new = test %>% as.data.frame() %>% as.matrix() %>% magrittr::set_names(NULL)
  # print(time_test)
  # user   system  elapsed
  # 1403.340   30.310  207.518

  return(
    Y_new2
  )
}

# detect_bp__update_ref__sub_updateY <- function(sam,X,Y,Z,BPs2.level1.plus,BPs2.level2){
#   print(sam)
# }

# BPs1 <- BPs2.level1.plus[[sam]]
# BPs2 <- BPs2.level2[[sam]]

# (1:ncol(X)) %>%
# lapply(function(sam) {list(sampleID=colnames(X)[sam])})

# segmenting
detect_bp__update_ref__sub_updateY <- function(sam,X_input,Y,Z,BPs2.level1.plus,BPs2.level2){

  # time3 = system.time({
  #   for (sam in 1:ncol(X)) {
  # sam <- 14
  # par(mfrow=c(3,3))
  # Y_new = copy(Y)
  sampleID <- colnames(X_input)[sam]
  # sampleID = input$SampleID
  n <- NCOL(X_input)
  d <- NROW(X_input)
  testdata <- data.frame(x=X_input[,sam], y=Y[,sam],z=Z[,sam])
  Ydiff_sam <- diff(testdata$y)

  BPs1 <- BPs2.level1.plus[[sam]]
  BPs2 <- BPs2.level2[[sam]]

  BPs <- unique(c(BPs1, BPs2))
  # BPs <- newBPs
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

    segbps <- sapply(seq_len(nseg), FUN=function(i) segment.table$begin[i]:segment.table$end[i])
    #first and last segments are connected
    segbps[[1]] <- c(segbps[[1]], segment.table$begin[nseg2]:segment.table$end[nseg2])
    sapply(segbps, length)

    ## Find a baseline
    Y.tmp <- Y
    Q1.zscore <- max.zscore <- mean.zscore.outwinreg <- rep(1000,nseg)
    med.outwinreg <- len.outwinreg <- rep(0,nseg)
    names(mean.zscore.outwinreg) <- seq_len(nseg)
    # par(mfrow=c(3,3))
    for (js in seq_len(nseg)) {
      # js <- 4
      outwinreg <- segbps[[js]]
      len.outwinreg[js] <- length(outwinreg)
      # if (length(outwinreg)>200) {
      newmed <- median(Y[outwinreg,sam])
      med.outwinreg[js] <- newmed
      if (newmed > 0.05) {
        Y.tmp[,sam] <- Y[,sam]/newmed
        Z.tmp <- t(apply(Y.tmp,1,function(x) pd.rate.hy(x,qrsc=TRUE)))
        mean.zscore.outwinreg[js] <- mean(abs(Z.tmp[outwinreg,sam]))
        Q1.zscore[js] <- quantile(abs(Z.tmp[,sam]),probs=0.25)
        max.zscore[js] <- max(abs(Z.tmp[,sam]))
      }
      # }
    }
    # mean.zscore.outwinreg
    # Q1.zscore
    # max.zscore
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
    # length(outwinreg)
    newmed <- median(Y[outwinreg,sam])
    return(Y[,sam]/newmed)  #######################################################
    # Z <- t(apply(Y,1,function(x) pd.rate.hy(x,qrsc=TRUE))) ###############################################
  }else{
    return(Y[,sam])
  }
  message(sam, "\n")
  # }
  # })

}




#' Final results
#'
#' @param X raw depth matrix
#' @param Y detect_bp__update_ref output
#'
#' @return list containing segmentation result
#' @noRd
segmentation_2nd_phase <- function(X,Y,cores=10){
  n <- NCOL(X)
  d <- NROW(X)

  Z <- t(apply(Y,1,function(x) pd.rate.hy(x,qrsc=TRUE)))
  segtable.list  <- vector("list", length=ncol(X))
  segtable2.list <- vector("list", length=ncol(X))

  # time41p = system.time({

  # for (sam in 1:ncol(X)) {
  # message("segtable.lists")
  segtable.lists <-
    seq_len(ncol(X)) %>%
    # .[-c(1:63)] %>%
    # 1:10 %>%
    mclapply(mc.cores=cores,function(sam){
      # lapply(function(sam){
      # print("sam")
      # print(sam)
      # cat(sam, "\n")
      sampleID <- colnames(X)[sam]

      # message("get_BPsegment_v2")
      segment.table <- get_BPsegment_v2(X=X, Y=Y,Z=Z, sampleID=sampleID)
      # segtable.list[[sam]] <-
      # message("compute_CN1")
      slst <-
        data.frame(segment.table, cn=compute_CN(segment.table))

      # message("smooth_segment")
      segment.table2 <- smooth_segment(segment.table)
      # segtable2.list[[sam]] <-
      #message("compute_CN2")
      slst2 <-
        data.frame(segment.table2, cn=compute_CN(segment.table2, base.state=1))


      return(list(slst=slst,slst2=slst2))
    })


  # 10 cores
  #     user   system  elapsed
  # 1770.612   31.677  268.775
  # time42p = system.time({
  for (sam in seq_len(ncol(X))) {
    sampleID <- colnames(X)[sam]
    segtable.list[[sam]] <- data.frame(id=sam, sampleID=colnames(X)[sam], segtable.lists[[sam]]$slst)
    segtable2.list[[sam]] <- data.frame(id=sam,sampleID=colnames(X)[sam], segtable.lists[[sam]]$slst2)
    # cat(sam, "\n")
  }
  # })
  # print(time41p)
  # print(time42p)

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
    X,
    N_cores=min(10L,parallel::detectCores()),
    reduced_output=TRUE #,
    #verbose=FALSE
  ){
    #if(verbose){
      run_ELViS_core(
        X
        ,N_cores=N_cores
        ,save_intermediate_data = FALSE
        ,save_idx=NULL
        ,save_dir="save_dir"
        ,overwrite=FALSE
        ,reduced_output=reduced_output
      )

    #}else{
      #suppressMessages(
        run_ELViS_core(
          X
          ,N_cores=N_cores
          ,save_intermediate_data = FALSE
          ,save_idx=NULL
          ,save_dir="save_dir"
          ,overwrite=FALSE
          ,reduced_output=reduced_output
        )
      #)


#    }


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
#'
#' @return list containing ELViS run results
#' @noRd
run_ELViS_core <- function(X,N_cores=min(10,parallel::detectCores())
                           # ,q.cutoff=4
                           ,save_intermediate_data = FALSE
                           ,save_idx=NULL,save_dir="save_dir",overwrite=FALSE,reduced_output=TRUE
                           # ,verbose=TRUE
){

  is_save <- save_intermediate_data
  if(is_save){
    dir.create(save_dir)
    if(is.null(save_idx)){save_idx <- 1}
  }

  message("ELViS run starts.")

  Name<-"normalized_data"
  q.cutoff<-4
  # message(Name)
  if(is_save){
    rds_file <- glue("{save_dir}/{Name}_{save_idx}.rds")
    # message(rds_file)
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
  # message(Name)
  if(is_save){
    rds_file <- glue("{save_dir}/{Name}_{save_idx}.rds")
    # message(rds_file)
    if(overwrite|(!file.exists(rds_file))){
      refupate_data <-
        update_reference_segments(normalized_data,cores = N_cores)
      saveRDS(refupate_data,rds_file)
    }else{
      refupate_data <- readRDS(rds_file)
    }
  }else{
    refupate_data <-
      update_reference_segments(normalized_data,cores = N_cores)
  }

  Name<-"refupate_data_secondrun"
  # message(Name)
  if(is_save){
    rds_file <- glue("{save_dir}/{Name}_{save_idx}.rds")
    # message(rds_file)
    if(overwrite|(!file.exists(rds_file))){
      refupate_data_secondrun <-
        update_reference_segments(
          # refupate_data %>%
          # with(list(
          list(
            sample_Ids=refupate_data$sample_Ids,
            Y=refupate_data$Y_recentered,
            Z=refupate_data$Z_recentered
            # )),
          ),
          cores = N_cores
        )
      saveRDS(refupate_data_secondrun,rds_file)
    }else{
      refupate_data_secondrun <- readRDS(rds_file)
    }
  }else{
    refupate_data_secondrun <-
      update_reference_segments(
        # refupate_data %>%
        # with(list(
        list(
          sample_Ids=refupate_data$sample_Ids,
          Y=refupate_data$Y_recentered,
          Z=refupate_data$Z_recentered
          # )),
        ),
        cores = N_cores
      )
  }

  message("Normalization done.")


  Name<-"refupate_data_secondrun_seg_cls"
  # message(Name)
  if(is_save){
    rds_file <- glue("{save_dir}/{Name}_{save_idx}.rds")
    # message(rds_file)
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
  # message(Name)
  if(is_save){
    rds_file <- glue("{save_dir}/{Name}_{save_idx}.rds")
    # message(rds_file)
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
  # message(Name)
  if(is_save){
    rds_file <- glue("{save_dir}/{Name}_{save_idx}.rds")
    # message(rds_file)
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
        half.width = 250
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
      half.width = 250
    )
  }


  Name<-"segtable_p2"
  # message(Name)
  if(is_save){
    rds_file <- glue("{save_dir}/{Name}_{save_idx}.rds")
    # message(rds_file)
    if(overwrite|(!file.exists(rds_file))){
      segtable_p2 <- segmentation_2nd_phase(X,Y=new_Y_p2,cores=N_cores)
      saveRDS(segtable_p2,rds_file)
    }else{
      segtable_p2 <- readRDS(rds_file)
    }
  }else{
    segtable_p2 <- segmentation_2nd_phase(X,Y=new_Y_p2,cores=N_cores)
  }

  message("Segmentation done.")


  final_output <- rbindlist(segtable_p2$segtable2.list )
  final_call_tmp <- dplyr::group_by(
    .data$final_output,
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
