#' Get performance from the contingency table
#'
#' @param tbl 2 x 2 table with columns indicating predicted classes and rows indicating actual classes and the values indicating the number of cases. 1st row and column should indicate the target true event.
#'
#' @return performance metrics
#' @noRd
get_perf <- function(tbl){
    rs <- tbl |> apply(1,sum)
    cs <- tbl |> apply(2,sum)
    Recall <- Sensitivity <- tbl[1,1]/rs[1]
    Specificity <- tbl[2,2]/rs[2]
    Precision <- PPV <- tbl[1,1]/cs[1] # positive predictive value or precision
    NPV <- tbl[2,2]/cs[2]
    F1 <- (2*Precision*Recall)/(Precision+Recall)
    Accuracy <- sum(diag(tbl))/sum(tbl)

    out <- data.frame(
        Sensitivity = Sensitivity,
        Specificity = Specificity,
        PPV = PPV, # positive predictive value or precision
        NPV = NPV,
        Accuracy = Accuracy,
        F1 = F1
    )

    return(out)

}





#' Check if decimal places are trivially small, i.e. if it is integer
#'
#' @param x a numeric to test
#' @param th If `x - as.integer(x)` is less than th, this function returns TRUE
#'
#' @return logical indicating if input number is trivially small
#' @noRd
#' @examples
#'
#' is_zero_dec(1.3) # FALSE
#' is_zero_dec(1.01) # FALSE
#' is_zero_dec(1.3-1) # TRUE
#'
is_zero_dec <- function(x,th=10^(-10)){
    (x - as.integer(x)) < th
}




#' Internal function for get_A_B_Prop_Subcount
#'
#' @param A A count
#' @param A_str A name
#' @param B B count
#' @param B_str B name
#' @param target_prop desired proportions of A (= A/(A + B))
#' @param n_cycle_th maximum number of revision cycles
#'
#' @return list containint updated A and B that satisfies target proportion
#' @noRd
get_A_B_Prop_Subcount__internal <-
    function(A,A_str,B,B_str,target_prop,n_cycle_th = 10){
        n_cycle <- 1
        while( n_cycle <= n_cycle_th ){

            msg <- paste0("n_cycle : ",n_cycle)
            message(msg)

            A <- B*target_prop/(1-target_prop)
            if(is_zero_dec(A)){break}
            message(A_str)
            A <- as.integer(A)

            B <- A/target_prop - A
            if(is_zero_dec(B)){break}
            message(B_str)
            B <- as.integer(B)

            n_cycle <- n_cycle + 1
        }
        return(
            structure(  list(A,B),
                        names = c(A_str,B_str)
            )

        )
    }


#' Get an exact integer number pair with starting numbers A and B that satisfies target proportion
#'
#' @param A target number A
#' @param B target number B
#' @param target_prop target proportion the user want to achieve in the form of A/(A+B).
#' @param n_cycle_th number of adjustment cycle to make both A and B integer-valued
#'
#' @return list containint updated A and B that satisfies target proportion
#' @noRd
#'
#' @examples
#' N_alt_ori <- 50L
#' N_ctrl_ori <- 600L
#' get_A_B_Prop_Subcount(A=N_alt_ori,B=N_ctrl_ori,target_prop=0.01)
#' get_A_B_Prop_Subcount(A=N_alt_ori,B=N_ctrl_ori,target_prop=0.05)
#' get_A_B_Prop_Subcount(A=N_alt_ori,B=N_ctrl_ori,target_prop=0.1)
#' get_A_B_Prop_Subcount(A=N_alt_ori,B=N_ctrl_ori,target_prop=0.2)
get_A_B_Prop_Subcount <-
    function(A,B,target_prop,n_cycle_th = 10){
        A_str <- deparse(substitute(A))
        B_str <- deparse(substitute(B))

        ori_prop <- A/(A + B)
        if(ori_prop > target_prop){
            output <-
                get_A_B_Prop_Subcount__internal(A=A,A_str = A_str,B=B,B_str = B_str,target_prop = target_prop,n_cycle_th=n_cycle_th)
        }else if(ori_prop < target_prop){
            output <-
                get_A_B_Prop_Subcount__internal(A=B,A_str = B_str,B=A,B_str = A_str,target_prop = 1-target_prop,n_cycle_th=n_cycle_th) %>%
                rev
        }else{
            message("No change")
            output <- structure(list(A,B),
                                names = c(A_str,B_str))
        }
        return(output)
    }









generate_palettes<-FALSE
if(!generate_palettes){
    choi_discrete_palettes <-
        list(
            Rushmore1 = c("#E1BD6D", "#EABE94", "#0B775E", "#35274A", "#F2300F")
            ,Royal1 = c("#899DA4", "#C93312", "#FAEFD1", "#DC863B")
            ,Royal2 = c("#9A8822", "#F5CDB4", "#F8AFA8", "#FDDDA0", "#74A089")
            ,Zissou1 = c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00")
            ,Darjeeling1 = c("#FF0000", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6")
            ,Darjeeling2 = c("#ECCBAE", "#046C9A", "#D69C4E", "#ABDDDE", "#000000")
            ,Chevalier1 = c("#446455", "#FDD262", "#D3DDDC", "#C7B19C")
            ,FantasticFox1 = c("#DD8D29", "#E2D200", "#46ACC8", "#E58601", "#B40F20")
            ,Moonrise2 = c("#798E87", "#C27D38", "#CCC591", "#29211F")
            ,Moonrise3 = c("#85D4E3", "#F4B5BD", "#9C964A", "#CDC08C", "#FAD77B")
            ,Cavalcanti1 = c("#D8B70A", "#02401B", "#A2A475", "#81A88D", "#972D15")
            ,GrandBudapest1 = c("#F1BB7B", "#FD6467", "#5B1A18", "#D67236")
            ,GrandBudapest2 = c("#E6A0C4", "#C6CDF7", "#D8A499", "#7294D4")
            ,aurora = c("#BF616A", "#D08770", "#EBCB8B", "#A3BE8C", "#B48EAD")
            ,Accent = c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17", "#666666")
            ,Dark2 = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
            ,Pastel1 = c("#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2")
            ,Pastel2 = c("#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC", "#CCCCCC")
            ,Set1 = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
            ,Set2 = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")
            ,Set3 = c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")
        )
    choi_paired_palettes <-
        list(
            blue = c("#A6CEE3", "#1F78B4")
            ,green = c("#B2DF8A", "#33A02C")
            ,red = c("#FB9A99", "#E31A1C")
            ,orange = c("#FDBF6F", "#FF7F00")
            ,purple = c("#CAB2D6", "#6A3D9A")
            ,brown = c("#FFFF99", "#B15928")
        )

    choi_ordered_palettes <-
        list(
            lumina = c("#EDDAEB", "#AD8CAE", "#4F93B8", "#306489", "#222B4C")
            ,silver_mine = c("#4B644B", "#647D4B", "#E1E1E1", "#7D96AF", "#647D96")
            ,lake_superior = c("#7D4B19", "#C89664", "#C87d4B", "#4B647D", "#324B64", "#19324B")
            ,victory_bonds = c("#AF1900", "#C83200", "#E19600", "#193264", "#001964")
            ,halifax_harbor = c("#E1C8AF", "#C8AF96", "#AF967D", "#967D7D", "#644B64", "#4B324b")
            ,moose_pond = c("#4B3232", "#7D4B32", "#966432", "#AF7D32", "#E19632", "#E1AF4B", "#C8C896", "#4B4B4B")
            ,algoma_forest = c("#4B4B4B", "#967D4B", "#AFAF7D", "#C89632", "#647D64", "#96AFAF", "#7D96AF")
            ,red_mountain = c("#7D3232", "#7D4B4B", "#7D6464", "#AF967D", "#FAC87D", "#E1AF64", "#C8964B", "#32324B")
            ,afternoon_prarie = c("#486090", "#6078A8", "#7890A8", "#90A8C0", "#F0D8C0", "#D6BBCF", "#A8C0C0", "#C0D8D8", "#A8A890")
        )

    choi_continuous_palettes <-
        list(
            RdYlBu = c("#A50026", "#FFFFBF", "#313695")
            ,RdGy = c("#67001F", "#FFFFFF", "#1A1A1A")
            ,PRGn = c("#40004B", "#F7F7F7", "#00441B")
            ,PiYG = c("#8E0152", "#F7F7F7", "#276419")
            ,BrBG = c("#543005", "#F5F5F5", "#003C30")
            ,YlGB = c("yellow", "grey", "black")
            ,BlGn = c("#1a1334", "#26294a", "#01545a", "#017351", "#03c383", "#aad962")
            ,PuYl = c("#110141", "#710162", "#a12a5e", "#ed0345", "#ef6a32", "#fbbf45")
            ,BuGn = c("#3e71a8", "#577f9f", "#698e96", "#779d8d", "#84ad83", "#8fbd77", "#99cd6b", "#a2dd5c", "#aaee49", "#b2ff2e")
            ,rainbow = c("#D12600", "#DB6A00", "#B2FF2E", "#00AD00", "#9CCADE", "#005B94", "#1E2085", "#610052", "#953272")
        )



    candicol1 <-
        c("#FFFFCC", "#FFF2AE", "#FFFFE5", "#FFFFCC", "#FFEDA0", "#FFF5F0", "#FFF7F3", "#FFF7EC", "#FEE8C8", "#FFF5EB", "#FEE6CE", "aliceblue"
        ) # candidate colors for exonic regions

    candicol2 <-
        c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F"
        ) # candidiate colors for regions with shape changes

    candicol3 <-
        c("#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC", "#CCCCCC"
        ) # candidiate colors for regions with shape changes

    candicol <- c(candicol2,candicol3);
    exon.col <- candicol1[9]

    col_yarrr_info2  <- c("#006A40FF", "#F08892FF", "#75B41EFF", "#95828DFF", "#708C98FF", "#8AB8CFFF", "#007E7FFF", "#358359FF", "#8BA1BCFF", "#5A5895FF", "#F2990CFF", "#5A5895FF", "#E5BA3AFF", "#D86C4FFF")

}else{
    pkgs_undetected <- c()
    for(pkg in c("nord","wesanderson","RColorBrewer")){
        if (!requireNamespace(pkg, quietly = TRUE)) {
            pkgs_undetected <- c(pkgs_undetected,pkg)
        }
    }

    if(length(pkgs_undetected)>0){
        pkgs_undetected_str <- paste0("'",paste(pkgs_undetected,collapse="', '"),"'")
        stop(glue("The '{pkgs_undetected_str}' package is required but not installed. Please install it."))
    }


    if(!requireNamespace("nord"))
        choi_discrete_palettes <- c(wes_palettes[ c("Rushmore1","Royal1","Royal2","Zissou1",
                                                    "Darjeeling1","Darjeeling2",
                                                    "Chevalier1","FantasticFox1","Moonrise2",
                                                    "Moonrise3","Cavalcanti1",
                                                    "GrandBudapest1","GrandBudapest2")],
                                    nord_palettes["aurora"])
    choi_discrete_palettes$Accent <- RColorBrewer::brewer.pal(8,"Accent")
    choi_discrete_palettes$Dark2 <- RColorBrewer::brewer.pal(8,"Dark2")
    choi_discrete_palettes$Pastel1 <- RColorBrewer::brewer.pal(9,"Pastel1")
    choi_discrete_palettes$Pastel2 <- RColorBrewer::brewer.pal(8,"Pastel2")
    choi_discrete_palettes$Set1 <- RColorBrewer::brewer.pal(9,"Set1")
    choi_discrete_palettes$Set2 <- RColorBrewer::brewer.pal(8,"Set2")
    choi_discrete_palettes$Set3 <- RColorBrewer::brewer.pal(12,"Set3")

    choi_paired_palettes <- as.list(NULL)
    choi_paired_palettes$blue <- RColorBrewer::brewer.pal(12,"Paired")[c(1,2)]
    choi_paired_palettes$green <- RColorBrewer::brewer.pal(12,"Paired")[c(3,4)]
    choi_paired_palettes$red <- RColorBrewer::brewer.pal(12,"Paired")[c(5,6)]
    choi_paired_palettes$orange <- RColorBrewer::brewer.pal(12,"Paired")[c(7,8)]
    choi_paired_palettes$purple <- RColorBrewer::brewer.pal(12,"Paired")[c(9,10)]
    choi_paired_palettes$brown <- RColorBrewer::brewer.pal(12,"Paired")[c(11,12)]

    choi_ordered_palettes <- nord_palettes[c(5,7,8,9,10,11,12,14,16)]
    # lumina, silver_mine, lake_superior, victory_bonds, halifax_harbor,
    # moose_pond, algoma_forest, red_mountain, afternnon_prairie
    choi_continuous_palettes <- as.list(NULL)
    choi_continuous_palettes$RdYlBu <- RColorBrewer::brewer.pal(11,"RdYlBu")[c(1,6,11)]
    choi_continuous_palettes$RdGy <- RColorBrewer::brewer.pal(11,"RdGy")[c(1,6,11)]
    choi_continuous_palettes$PRGn <- RColorBrewer::brewer.pal(11,"PRGn")[c(1,6,11)]
    choi_continuous_palettes$PiYG <- RColorBrewer::brewer.pal(11,"PiYG")[c(1,6,11)]
    choi_continuous_palettes$BrBG <- RColorBrewer::brewer.pal(11,"BrBG")[c(1,6,11)]
    choi_continuous_palettes$YlGB <- c( "yellow","grey","black")
    choi_continuous_palettes$BlGn <- c( "#1a1334", "#26294a", "#01545a", "#017351", "#03c383",
                                        "#aad962")
    choi_continuous_palettes$PuYl <- c( "#110141", "#710162", "#a12a5e", "#ed0345", "#ef6a32",
                                        "#fbbf45")
    choi_continuous_palettes$BuGn <- c( "#3e71a8", "#577f9f", "#698e96", "#779d8d", "#84ad83",
                                        "#8fbd77", "#99cd6b", "#a2dd5c", "#aaee49", "#b2ff2e")
    choi_continuous_palettes$rainbow <- c(  rosso_corsa = "#D12600", spanish_orange = "#DB6A00",
                                            green_yellow = "#B2FF2E", green = "#00AD00", pale_cerulean = "#9CCADE",
                                            sea_blue = "#005B94", st_patricks_blue = "#1E2085", tyrian_purple = "#610052",
                                            amaranth_deep_purple = "#953272")






    candicol1 <- c( RColorBrewer::brewer.pal(9,"Pastel1")[6], # candidate colors for exonic regions
                    RColorBrewer::brewer.pal(8,"Pastel2")[6],
                    RColorBrewer::brewer.pal(9,"YlOrBr")[1],
                    RColorBrewer::brewer.pal(9,"YlOrRd")[1],
                    RColorBrewer::brewer.pal(9,"YlOrRd")[2],
                    RColorBrewer::brewer.pal(9,"Reds")[1],
                    RColorBrewer::brewer.pal(9,"RdPu")[1],
                    RColorBrewer::brewer.pal(9,"OrRd")[1],
                    RColorBrewer::brewer.pal(9,"OrRd")[2],
                    RColorBrewer::brewer.pal(9,"Oranges")[1],
                    RColorBrewer::brewer.pal(9,"Oranges")[2],
                    "aliceblue");
    candicol2 <- RColorBrewer::brewer.pal(12,"Set3") # candidiate colors for regions with shape changes
    candicol3 <- RColorBrewer::brewer.pal(8,"Pastel2") # candidiate colors for regions with shape changes
    candicol <- c(candicol2,candicol3);
    exon.col <- candicol1[9]


    col_yarrr_info2 <- piratepal(palette = "info2")



}











#' from SCISSOR package https://github.com/hyochoi/SCISSOR/blob/377789106d34200baa1d4d826952e2ad5313ba3a/R/compScales.R#L4
#' @noRd
compScales <- function( x,
                        rmZeroes=FALSE,maxRatio=NULL,precScale=1e-10){
    # Computes the scales sa and sb (above and below the median).
    # Assumes that x is an array of numbers.
    #
    x <- x[!is.na(x)] # we always take out NAs
    temp <- fastSplitSample(x)
    xa   <- temp$xa
    xb   <- temp$xb
    med  <- temp$med
    sall <- scale1StepM((x-med),precScale=precScale)
    if(rmZeroes){ # reduces breakdown value but yields fewer implosions
        xa <- xa[xa > precScale]
        xb <- xb[xb > precScale]
    }
    sa <- scale1StepM(xa,precScale=precScale)
    sb <- scale1StepM(xb,precScale=precScale)
    if(!is.null(maxRatio)){
        if(maxRatio < 2) stop("maxRatio must be at least 2")
        sa <- min(  c(max(sa,sall/maxRatio,na.rm = TRUE),sall*maxRatio),
                    na.rm = TRUE)
        sb <- min(  c(max(sb,sall/maxRatio,na.rm=TRUE),sall*maxRatio),
                    na.rm = TRUE)
    }
    return(list(sa=sa,sb=sb,med=med))
}


#' from SCISSOR package https://github.com/hyochoi/SCISSOR/blob/377789106d34200baa1d4d826952e2ad5313ba3a/R/fastSplitSample.R#L2
#' @noRd
fastSplitSample <- function(x){
    # Centers sample by median, and divides in 2 equal halves.
    # Assumes that NAs have already been removed.
    # This function has time complexity O(n).
    #
    med <- median(x) # takes only O(n) time
    x <- x - med # centering
    n <- length(x)
    h <- n %/% 2   # = integer part of n/2
    xa <- x[x > 0] # length(xa) <= h
    xb <- x[x < 0] # length(xa) <= h
    xa <- c(rep(0,(n - h - length(xa))),xa)
    xb <- c(rep(0,(n - h - length(xb))),abs(xb)) # abs() !
    return(list(xa=xa,xb=xb,med=med))
}


#' Projection Depth from R package SCISSOR : https://github.com/hyochoi/SCISSOR/blob/377789106d34200baa1d4d826952e2ad5313ba3a/R/pd.rate.hy.R
#' @noRd
pd.rate.hy <- function(x,qrsc=FALSE) {
    # projection depth
    m <- mad(x)
    if (m<1e-5) {
        return(rep(0,length(x)))
    } else {
        if (qrsc) {
            rsc<-compScales(x)
            y<-rep(0,length(x))

            above.ind<-which((x-rsc$med)>=0)
            below.ind<-which((x-rsc$med)<0)
            if (rsc$sa>1e-5) {
                y[above.ind]<-(x[above.ind]-rsc$med)/rsc$sa
            }
            if (rsc$sb>1e-5) {
                y[below.ind]<-(x[below.ind]-rsc$med)/rsc$sb
            }
            return(y)
        } else {
            return((x-median(x))/m)
        }
    }
}



#' from SCISSOR https://github.com/hyochoi/SCISSOR/blob/377789106d34200baa1d4d826952e2ad5313ba3a/R/rhoHuber.R#L4
#' @noRd
rhoHuber <- function(x,c=2.1){
    # x is a univariate sample
    # c is the tuning constant
    # output is rho(x)
    #
    rho <- (x/c)^2
    rho[rho > 1] <- 1
    1.54^2*rho
}

#' https://github.com/hyochoi/SCISSOR/blob/377789106d34200baa1d4d826952e2ad5313ba3a/R/scale1StepM.R#L4
#' @noRd
scale1StepM <- function(x,precScale=1e-10) {
    # Computes the first step of an algorithm for
    # a scale M-estimator using the given rho function.
    # The scatter is computed relative to zero.
    #
    x <- x[!is.na(x)] # we always take out NAs
    n <- length(x)
    if(n == 0) { return(0.0)
    } else {
        sigma0 <- 1.4826*median(abs(x))
        if(sigma0 < precScale) { return(0.0)
        } else {
            rho <- rhoHuber(x/sigma0)
            return(sigma0 * sqrt(sum(rho)*2/n))
        }
    }
}

#' from R package SCISSOR  : https://github.com/hyochoi/SCISSOR/blob/377789106d34200baa1d4d826952e2ad5313ba3a/R/yaxis.hy.R
#' @noRd
yaxis.hy <- function(mat){
    #  mat : d by n matrix
    tempmax <- max(mat) ;
    tempmin <- min(mat) ;
    templen <- tempmax-tempmin ;
    return(c(tempmin-0.002*templen, tempmax+0.002*templen)) ;
}





#' @noRd
.onLoad <- function(libname, pkgname) {
    if (is.null(getOption("repos")) || getOption("repos")["CRAN"] == "@CRAN@") {
        options(repos = c(CRAN = "https://cloud.r-project.org"))
    }
}
