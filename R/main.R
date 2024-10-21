# source("R/anno_color.R")

require(nord)
require(RColorBrewer)
require(wesanderson)
require(jcolors)
# jcolors_contin
choi_discrete_palettes <- c(wes_palettes[c("Rushmore1","Royal1","Royal2","Zissou1",
                                           "Darjeeling1","Darjeeling2",
                                           "Chevalier1","FantasticFox1","Moonrise2",
                                           "Moonrise3","Cavalcanti1",
                                           "GrandBudapest1","GrandBudapest2")],
                            nord_palettes["aurora"])
choi_discrete_palettes$Accent <- brewer.pal(8,"Accent")
choi_discrete_palettes$Dark2 <- brewer.pal(8,"Dark2")
choi_discrete_palettes$Pastel1 <- brewer.pal(9,"Pastel1")
choi_discrete_palettes$Pastel2 <- brewer.pal(8,"Pastel2")
choi_discrete_palettes$Set1 <- brewer.pal(9,"Set1")
choi_discrete_palettes$Set2 <- brewer.pal(8,"Set2")
choi_discrete_palettes$Set3 <- brewer.pal(12,"Set3")

choi_paired_palettes <- as.list(NULL)
choi_paired_palettes$blue <- brewer.pal(12,"Paired")[c(1,2)]
choi_paired_palettes$green <- brewer.pal(12,"Paired")[c(3,4)]
choi_paired_palettes$red <- brewer.pal(12,"Paired")[c(5,6)]
choi_paired_palettes$orange <- brewer.pal(12,"Paired")[c(7,8)]
choi_paired_palettes$purple <- brewer.pal(12,"Paired")[c(9,10)]
choi_paired_palettes$brown <- brewer.pal(12,"Paired")[c(11,12)]

choi_ordered_palettes <- nord_palettes[c(5,7,8,9,10,11,12,14,16)]
# lumina, silver_mine, lake_superior, victory_bonds, halifax_harbor,
# moose_pond, algoma_forest, red_mountain, afternnon_prairie
choi_continuous_palettes <- as.list(NULL)
choi_continuous_palettes$RdYlBu <- brewer.pal(11,"RdYlBu")[c(1,6,11)]
choi_continuous_palettes$RdGy <- brewer.pal(11,"RdGy")[c(1,6,11)]
choi_continuous_palettes$PRGn <- brewer.pal(11,"PRGn")[c(1,6,11)]
choi_continuous_palettes$PiYG <- brewer.pal(11,"PiYG")[c(1,6,11)]
choi_continuous_palettes$BrBG <- brewer.pal(11,"BrBG")[c(1,6,11)]
choi_continuous_palettes$YlGB <- c("yellow","grey","black")
choi_continuous_palettes$BlGn <- c("#1a1334", "#26294a", "#01545a", "#017351", "#03c383",
                                   "#aad962")
choi_continuous_palettes$PuYl <- c("#110141", "#710162", "#a12a5e", "#ed0345", "#ef6a32",
                                   "#fbbf45")
choi_continuous_palettes$BuGn <- c("#3e71a8", "#577f9f", "#698e96", "#779d8d", "#84ad83",
                                   "#8fbd77", "#99cd6b", "#a2dd5c", "#aaee49", "#b2ff2e")
choi_continuous_palettes$rainbow <- c(rosso_corsa = "#D12600", spanish_orange = "#DB6A00",
                                      green_yellow = "#B2FF2E", green = "#00AD00", pale_cerulean = "#9CCADE",
                                      sea_blue = "#005B94", st_patricks_blue = "#1E2085", tyrian_purple = "#610052",
                                      amaranth_deep_purple = "#953272")

anno_color_var <- function(variable, palette=NULL, palette_id=NULL, palette_order=NULL,
                           paired=FALSE, cont_midval=NULL) {
  require(circlize)
  variable0 <- variable
  # variable <- variable[! is.na(variable)]
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
    if (length(levels(variable))>2) paired = FALSE
    if (paired) {
      if (is.null(palette) & is.null(palette_id)) palette_id = 1
      if (is.null(palette)) palette <- names(choi_paired_palettes)[palette_id]
      colors <- choi_paired_palettes[[eval(palette)]]
      var_col <- colors[1:length(levels(variable))]
      names(var_col) <- levels(variable)
    } else {
      if (sum(class(variable) %in% c("ordered")) > 0) {
        if (is.null(palette) & is.null(palette_id)) palette_id = 3
        if (is.null(palette)) palette <- names(choi_ordered_palettes)[palette_id]
        colors <- choi_ordered_palettes[[eval(palette)]]
        if (! is.null(palette_order)) {
          colors <- colors[palette_order]
        }
        var_col <- colors[1:length(levels(variable))]
        names(var_col) <- levels(variable)
      } else {
        if (is.null(palette) & is.null(palette_id)) palette_id = 19
        if (is.null(palette)) palette <- names(choi_discrete_palettes)[palette_id]
        colors <- choi_discrete_palettes[[eval(palette)]]
        if (! is.null(palette_order)) {
          colors <- colors[palette_order]
        }
        var_col <- colors[1:length(levels(variable))]
        names(var_col) <- levels(variable)
      }
    }
    return(list(var=variable0, var_col=var_col,
                col=as.character(factor(variable0, levels=names(var_col), labels=var_col))))
  }
}



#' @import zoo
#' @noRd
detect_BPs <- function(y, half.width=250, prob.cutoff=0.95, q.Y=5) {
  # y <- Ydiff[,sam]
  require(zoo)
  # rollYquant0 <- rollapply(y,width=2*half.width,FUN=function(x) quantile(x,probs=prob.cutoff))
  rollYquant0 <- rollapply(abs(y),width=2*half.width,FUN=function(x) quantile(x,probs=prob.cutoff))
  rollYquant <- c(rep(rollYquant0[1],half.width), rollYquant0, rep(tail(rollYquant0,1),half.width-1))

  ## remove consecutive BPs
  Yexceed <- abs(y)/abs(q.Y*rollYquant)
  BPcandi <- which(Yexceed > 1)
  # breakpoint candidates : absolute difference larger than 5 times of Q95 among +- half.width positions

  if (length(BPcandi) > 2) {
    BPcandi.info <- data.frame(BPs=BPcandi[2:length(BPcandi)],
                               BPs.sign=y[1:(length(BPcandi)-1)]*y[2:length(BPcandi)],
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
#' @import zoo
#' @import igraph
#' @import stringr
#' @noRd
get_BPsegment_v2 <- function(X, Y, Z=NULL, sampleID,
                             q.Y.BP1 = 5,
                             q.Y.BP2 = 10,
                             BP.ydepth = 0.08, BP.xdepth = 10) {
  require(zoo)
  if (is.numeric(sampleID)) {
    sam <- sampleID
  } else {
    sam <- which(colnames(X)==sampleID)
  }
  if (is.null(Z)) {
    Z <- t(apply(Y,1,function(x) pd.rate.hy(x,qrsc=T))) ########################
  }

  d = NROW(X)

  # Xdiff <- apply(X, 2, FUN=function(x) diff(x))
  # Ydiff <- apply(Y, 2, FUN=function(x) diff(x))
  xdiff <- diff(X[,sam])
  ydiff <- diff(Y[,sam])
  ymed <- apply(Y, 1, median) #######################################

  testdata <- data.frame(x=X[,sam], y=Y[,sam], z=Z[,sam])

  BPs.q10 <- detect_BPs(y=ydiff, q.Y=q.Y.BP2, half.width=250, prob.cutoff = 0.95)
  BPs.q05 <- detect_BPs(y=ydiff, q.Y=q.Y.BP1, half.width=250, prob.cutoff = 0.95)
  BPs.cutx <- which(abs(xdiff) > BP.xdepth)
  BPs.cuty <- which(abs(ydiff) > BP.ydepth)

  # BPs.cutx[order(-abs(xdiff[BPs.cutx]))<10]
  # xdiff[BPs.cutx][rank(-abs(xdiff[BPs.cutx]))<10]
  # BPs.cutx[rank(-abs(xdiff[BPs.cutx]))<10]

  BPs.level1 <- base::intersect(BPs.q05, BPs.cutx)
  BPs.level2 <- base::intersect(BPs.q10, BPs.cutx)
  BPs.level1.plus <- base::intersect(BPs.level1, BPs.cuty)

  BPs <- unique(c(BPs.level1.plus, BPs.level2))
  # ydiff[BPs]

  # BPs <- newBPs
  BPs.sam <- length(BPs)
  if (BPs.sam > 0) {
    if (BPs.sam %% 2 > 0) {
      BPs.sam <- BPs.sam + 1
      newBP <- c(c(1:length(ydiff))[-BPs])[which.max(abs(ydiff[-BPs]))]
      BPs <- c(BPs,newBP)

      # c(c(1:length(ydiff))[-BPs])[which.max(abs(ydiff[-BPs]))]
      # abs(ydiff)[4300:5500] %>% plot(4300:5500,.);abline(v=BPs,col="red")
      # abs(ydiff[-BPs])[4300:5500] %>% plot(4300:5500,.);abline(v=BPs -(1:3),col="red")
    }
    BPs <- sort(BPs)

    nseg <- length(BPs)
    nseg2 <- nseg + 1
    srt.BPs <- sort(BPs)
    segment.table <- data.frame(seg = 1:nseg2,
                                begin = c(1,srt.BPs+1),
                                end = c(srt.BPs, d))
    segbps <- sapply(1:nseg, FUN=function(i) segment.table$begin[i]:segment.table$end[i])
    segbps[[1]] <- c(segbps[[1]], segment.table$begin[nseg2]:segment.table$end[nseg2])

    ## Find a baseline
    Y.tmp <- Y
    Q1.zscore <- max.zscore <- mean.zscore.outwinreg <- rep(1000,nseg)
    med.outwinreg <- len.outwinreg <- rep(0,nseg)
    names(mean.zscore.outwinreg) <- 1:nseg
    for (js in 1:nseg) {
      outwinreg <- segbps[[js]]
      len.outwinreg[js] <- length(outwinreg)
      # if (length(outwinreg)>200) {
      newmed = median(Y[outwinreg,sam])
      med.outwinreg[js] <- newmed
      if (newmed > 0.05) {
        Y.tmp[,sam] = Y[,sam]/newmed
        Z.tmp = t(apply(Y.tmp,1,function(x) pd.rate.hy(x,qrsc=T)))
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
    # baseseg

    ## Make a table for each segment
    seg_len_vec = sapply(segbps,FUN=length)
    segtable <- segment.table[1:nseg, ] # for now, drop the last segment (this is the same as the first segment)
    segtable$mu.x <- sapply(segbps, FUN=function(ibps) mean(testdata$x[ibps]))
    segtable$sd.x <- sapply(segbps, FUN=function(ibps) sd(testdata$x[ibps]))
    segtable$mu.y <- sapply(segbps, FUN=function(ibps) mean(testdata$y[ibps]))
    segtable$sd.y <- sapply(segbps, FUN=function(ibps) sd(testdata$y[ibps]))
    segtable$mu.z <- sapply(segbps, FUN=function(ibps) mean(testdata$z[ibps]))
    segtable$sd.z <- sapply(segbps, FUN=function(ibps) sd(testdata$z[ibps]))
    segtable$cohort.y <- sapply(segbps, FUN=function(ibps) mean(ymed[ibps]))

    # baseseg <- which.min(abs(segtable$mu.y-1))
    segtable$baseline <- rep(0,nseg)
    segtable$baseline[baseseg] <- 1

    ## Find more base segments
    baseseg2 <- which(abs(segtable$mu.y - segtable$mu.y[baseseg]) < segtable$sd.y[baseseg])
    baseseg3 <- which(abs(segtable$mu.z - segtable$mu.z[baseseg]) < segtable$sd.z[baseseg])
    segtable$baseline[base::intersect(baseseg2,baseseg3)] <- 1
    segtable$state <- segtable$baseline

    i_nonbase = which(segtable$baseline==0)
    n_nonbase = length(i_nonbase)
    if(n_nonbase!=0){
      if(n_nonbase==1){
        segtable$state[which(segtable$baseline==0)] <-
          seq_len(length(which(segtable$baseline==0))) + 1
      }else if(n_nonbase>=2){

        # is_1st_above = segtable$mu.y[i_nonbase][1] > segtable$mu.y[baseseg]
        # if(is_1st_above){
        #   is_1st_grp = segtable$mu.y[i_nonbase] > segtable$mu.y[baseseg]
        # }else{
        #   is_1st_grp = segtable$mu.y[i_nonbase] < segtable$mu.y[baseseg]
        # }
        is_above = segtable$mu.y[i_nonbase] > segtable$mu.y[baseseg]
        i_nonbase_sp =
          split(
            i_nonbase,is_above,drop = TRUE)

        # find close copynumber segment groups
        i_grouped_lst =
          i_nonbase_sp |>
          lapply(\(i_tmp){
            # print(i_tmp)

          # })
        # for(i_tmp in i_nonbase_sp){
          # print(i_tmp)
          if(length(i_tmp)==1){
            out = list(i_tmp)
          }else if(length(i_tmp)==2){
            # baseseg2 <- which(abs(segtable$mu.y[i_tmp] - segtable$mu.y[baseseg]) < segtable$sd.y[baseseg])
            # baseseg3 <- which(abs(segtable$mu.z - segtable$mu.z[baseseg]) < segtable$sd.z[baseseg])

            psd_y = pooled_sd(segtable$sd.y[i_tmp[1]],segtable$sd.y[i_tmp[2]],seg_len_vec[i_tmp[1]],seg_len_vec[i_tmp[2]])
            is_simil_y = abs(diff(segtable$mu.y[i_tmp])) < psd_y
            psd_z = pooled_sd(segtable$sd.z[i_tmp[1]],segtable$sd.z[i_tmp[2]],seg_len_vec[i_tmp[1]],seg_len_vec[i_tmp[2]])
            is_simil_z = abs(diff(segtable$mu.z[i_tmp])) < psd_z

            if(is_simil_y&is_simil_z){
              # segtable$state[i_tmp] <- max(segtable$state) + 1 # same group
              out = list(i_tmp)  # same seg
            }else{
              # segtable$state[i_tmp] <- max(segtable$state) + seq_along(i_tmp) # diff group
              out = as.list(i_tmp) # different seg
            }
          }else if(length(i_tmp)>=3){
            i_tmp_cmb = combinat::combn(i_tmp,2)
            max_rel_d = i_tmp_cmb |>
              apply(2,\(j_tmp){
                psd_y = pooled_sd(segtable$sd.y[j_tmp[1]],segtable$sd.y[j_tmp[2]],seg_len_vec[j_tmp[1]],seg_len_vec[j_tmp[2]])
                simil_y = abs(diff(segtable$mu.y[j_tmp]))/psd_y
                psd_z = pooled_sd(segtable$sd.z[j_tmp[1]],segtable$sd.z[j_tmp[2]],seg_len_vec[j_tmp[1]],seg_len_vec[j_tmp[2]])
                simil_z = abs(diff(segtable$mu.z[j_tmp]))/psd_z
                return(c(max(simil_y,simil_z) ))
              })

            g = gen_dist_graph(i_tmp,i_tmp_cmb,max_rel_d,th_d=1)


            subn_grps =
              get_subcl(g) |>
              lapply(str_replace,"^v","") |>
              lapply(as.numeric)

            out = subn_grps

          }

            return(out)
        }) |> unlist(recursive = FALSE)

        # update group
        for( k_tmp in  i_grouped_lst[order(sapply(i_grouped_lst,min))]){
          segtable$state[k_tmp] <- max(segtable$state) + 1
        }

      }
    }


    segtable2 <- rbind(segtable,segtable[1,])
    segtable2[nseg2,1:3] <- segment.table[nseg2,1:3]
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


# gen_dist_graph(i_tmp,i_tmp_cmb,max_rel_d,th_d=1)
gen_dist_graph = function(i_tmp,i_tmp_cmb,max_rel_d,th_d){
  max_rel_s = 1/(max_rel_d+1)

  is_close = max_rel_d < th_d


  g <- make_empty_graph(n = length(i_tmp),directed = FALSE)
  V(g)$name <- paste0("v",i_tmp)

  if(any(is_close)){
    edges = paste0("v",i_tmp_cmb[,is_close] |> as.vector())
    g <- add_edges(g, edges)  # Connect node 1 to 2 and node 3 to 4
    E(g)$weight <- max_rel_s[is_close]
  }


  return(g)
}


clique_weight <- function(clique, graph) {
  subg <- induced_subgraph(graph, vids = clique)
  return(sum(E(subg)$weight))
}

get_subcl = function(g){
  cliques_list <- cliques(g)

  clique_weights <- sapply(cliques_list, clique_weight, graph = g)

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

  out_subnetworks = subnetworks |> lapply(\(x) induced_subgraph(g,vids=x)) |> lapply(V) %>% lapply(names)
  return(out_subnetworks)
}



# get_window = function(Y,sam) {
#   ips = which((Y[,sam]-1)[1:(nrow(Y)-1)]*(Y[,sam]-1)[2:nrow(Y)]<=0)
#   win = matrix(c(1,rep(ips,each=2),dim(Y)[1]),ncol=2,byrow=T)
#   win = win[which((win[,2] - win[,1])>10),]
#   return(win)
# }

get_window2 = function(y,cutoff=1,min.length=10) {
  d = length(y)
  ips = which((y-cutoff)[1:(d-1)]*(y-cutoff)[2:d]<=0)
  win = matrix(c(1,rep(ips,each=2),d),ncol=2,byrow=T)
  win = win[which((win[,2] - win[,1])>min.length),]
  return(win)
}

get_window = function(y,pos=NULL,cutoff=1,min.length=10) {
  if (is.null(pos)) {
    return(get_window2(y=y,cutoff=cutoff,min.length=min.length))
  } else {
    wincan = get_window2(y=y,cutoff=cutoff,min.length=min.length)
    windiff = wincan - pos
    win.tmp = wincan[which((windiff[,1]*windiff[,2])<=0),]
    return(win.tmp)
  }
}



t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color

  ## Get RGB values for named color
  rgb.val <- col2rgb(color)

  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)

  ## Save the color
  invisible(t.col)
}

compute_CN <- function(segment.table, base.state=1) {
  if (is.null(base.state)) base.state=1
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

get_BPsegment <- function(X, Y, Z=NULL, sampleID,
                          BP.ydepth = 0.08) {
  if (is.numeric(sampleID)) {
    sam <- sampleID
  } else {
    sam <- which(colnames(X)==sampleID)
  }
  if (is.null(Z)) {
    Z <- t(apply(Y,1,function(x) pd.rate.hy(x,qrsc=T)))
  }
  Xdiff <- apply(X, 2, FUN=function(x) diff(x))
  Ydiff <- apply(Y, 2, FUN=function(x) diff(x))
  prob.cutoff <- 0.9
  BPnums1 <- apply(Ydiff,2,FUN=function(x) which(abs(x)>BP.ydepth))
  BPcutoffs <- apply(Ydiff, 2, FUN=function(x) 5*quantile(x,probs=prob.cutoff))
  BPnums2 <- apply(Ydiff,2,FUN=function(x) which(abs(x)>(5*quantile(x,probs=prob.cutoff))))

  testdata <- data.frame(x=X[,sam], y=Y[,sam],z=Z[,sam])
  ymed <- apply(Y, 1, median)
  BPs.sam <- length(base::intersect(BPnums1[[sam]], BPnums2[[sam]]))
  if (BPs.sam > 0) {
    if (BPs.sam %% 2 > 0) {
      BPs.sam <- BPs.sam + 1
    }
    BPs <- order(abs(Ydiff[,sam]),decreasing=TRUE)[1:min(BPs.sam, 10)]
    # Xdiff[BPs,sam]

    nseg <- length(BPs)
    nseg2 <- nseg + 1
    srt.BPs <- sort(BPs)
    segment.table <- data.frame(seg = 1:nseg2,
                                begin = c(1,srt.BPs+1),
                                end = c(srt.BPs, d))
    segbps <- sapply(1:nseg, FUN=function(i) segment.table$begin[i]:segment.table$end[i])
    segbps[[1]] <- c(segbps[[1]], segment.table$begin[nseg2]:segment.table$end[nseg2])

    ## Make a table for each segment
    segtable <- segment.table[1:nseg, ] # for now, drop the last segment (this is the same as the first segment)
    segtable$mu.x <- sapply(segbps, FUN=function(ibps) mean(testdata$x[ibps]))
    segtable$sd.x <- sapply(segbps, FUN=function(ibps) sd(testdata$x[ibps]))
    segtable$mu.y <- sapply(segbps, FUN=function(ibps) mean(testdata$y[ibps]))
    segtable$sd.y <- sapply(segbps, FUN=function(ibps) sd(testdata$y[ibps]))
    segtable$mu.z <- sapply(segbps, FUN=function(ibps) mean(testdata$z[ibps]))
    segtable$sd.z <- sapply(segbps, FUN=function(ibps) sd(testdata$z[ibps]))
    segtable$cohort.y <- sapply(segbps, FUN=function(ibps) mean(ymed[ibps]))

    baseseg <- which.min(abs(segtable$mu.y-1))
    segtable$baseline <- rep(0,nseg)
    segtable$baseline[baseseg] <- 1

    ## Find more base segments
    baseseg2 <- which(abs(segtable$mu.y - segtable$mu.y[baseseg]) < segtable$sd.y[baseseg])
    baseseg3 <- which(abs(segtable$mu.z - segtable$mu.z[baseseg]) < segtable$sd.z[baseseg])
    segtable$baseline[base::intersect(baseseg2,baseseg3)] <- 1
    segtable$state <- segtable$baseline
    segtable$state[which(segtable$baseline==0)] <- seq(1:length(which(segtable$baseline==0))) + 1

    segtable2 <- rbind(segtable,segtable[1,])
    segtable2[nseg2,1:3] <- segment.table[nseg2,1:3]
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
} # end of function

