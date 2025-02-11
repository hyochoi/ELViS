globalVariables(c(".", "%>%"))

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
#' @importFrom glue glue
#' @importFrom stringr str_to_title str_extract_all str_replace_all str_split str_replace
#'
#' @examples
#' data(mtrx_samtools_reticulate)
#' th <- 50
#' depth_hist(mtrx_samtools_reticulate,th=th,smry_fun=max)
#' depth_hist(mtrx_samtools_reticulate,th=th,smry_fun = quantile,prob=0.95)
#' depth_hist(mtrx_samtools_reticulate,th=th,smry_fun = mean)
depth_hist <- function(mtrx,th=50,title_txt=NULL,smry_fun=max,...){

    ## input checking
    stopifnot_mtrx_or_df(mtrx)
    stopifnot_numeric1(th)

    if(is.null(title_txt)){
        dot_info <- capture_params_glue(...)
        fun_name <- str_to_title(deparse(substitute(smry_fun)))
        if(length(dot_info)==0){
            title_txt <- (glue("{fun_name} Depth Distribution"))
        }else{
            title_txt <- (glue("{fun_name} {dot_info} Depth Distribution"))
        }

    }

    stopifnot_character1(title_txt)
    stopifnot(is(smry_fun,"function"))
    ## input checking done

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



#' @noRd
log10p1_trans <- trans_new(
    name = "log10p1",
    transform = function(x) log10(x + 1),
    inverse = function(x) 10^x - 1,
    breaks = function(x) log_breaks(base = 10)(x + 1)
)




#' @import ggplot2
#' @importFrom txdbmaker makeTxDbFromGFF
#' @importFrom GenomicFeatures transcriptsBy cdsBy genes
#' @rawNamespace import(GenomicRanges, except=c(subtract,intersect))
#' @rawNamespace import(IRanges, except=c(slice,mad,median,quantile,sd,intersect))
#' @rawNamespace import(BiocGenerics, except=c(mad,sd,combine,Position,normalize,path,image,density,setequal,intersect))
#' @noRd
get_gene_anno_plot_ori <- function(
        gff3_fn,
        space_length,
        annot_margin = 0.01,
        arrow_spacing = 0.05,
        geme_name_space = 0.5,
        col_pal =  col_yarrr_info2,
        exclude_genes = NULL
){

    # plot configure
    annot_half_height <- (1-geme_name_space)/2
    annot_halfmargin_bp <- as.integer(annot_margin*space_length/2)
    arrow_spacing_bp <- as.integer(arrow_spacing*space_length)

    # prep annotation
    gff_parsed <- parse_gff(gff3_fn,exclude_genes)
    txdb <- gff_parsed$txdb
    genes <- gff_parsed$gene
    cds <- gff_parsed$cds
    all_features <- gff_parsed$all_features

    genes_margin <- genes + annot_halfmargin_bp
    start(genes_margin) <- pmax(1,start(genes_margin))

    ylevels <-
        structure(
            rep(Inf,length(genes_margin)),
            names = names(genes)
        )

    ovlp_status <- genes_margin %>%
        findOverlaps() %>%
        sort %>%
        as.data.frame() %>%
        filter(.data$queryHits>.data$subjectHits)

    # non-overlapping -> level 1
    ylevels[!(seq_along(ylevels) %in% ovlp_status$queryHits)] <- 1

    for(ovlp in ovlp_status %>% group_split(.data$queryHits)){
        q_tmp <- ovlp$queryHits[1]
        s_tmp <- ovlp$subjectHits
        # take the lowest level among those not taken already by overlapping element
        ylevels[q_tmp] <- min(setdiff(seq_along(ylevels),ylevels[s_tmp]))

    }


    cds_plotdata <-
        names(cds) %>%
        lapply(\(gn){
            y_tmp <- ylevels[[gn]]
            cds_tmp <- cds[[gn]]
            s <- start(cds_tmp)
            e <-   end(cds_tmp)
            strnd <- as.vector(strand(cds_tmp))[1]
            cds_df <-
                data.frame(
                    Type    = "cds",
                    Gene    = gn,
                    Start   = s,
                    End     = e,
                    Strand  = strnd,
                    Y_Level = y_tmp
                )
            if(length(cds_tmp)==1){
                output <- cds_df
            }else{
                gap_tmp <- gaps(cds_tmp,start=min(s),end=max(e),ignore.strand=TRUE)
                gs <- start(gap_tmp)
                ge <-   end(gap_tmp)

                gap_df <-
                    data.frame(
                        Type    = "intron",
                        Gene    = gn,
                        Start   = gs,
                        End     = ge,
                        Strand  = strnd,
                        Y_Level = y_tmp
                    )
                output <-
                    rbind(
                        cds_df,
                        gap_df
                    )

            }

            return(output)

        }) %>%
        rbindlist()



    cds_plotdata__CDS  <- cds_plotdata%>% filter(.data$Type=="cds")
    cds_plotdata__INTRON <- cds_plotdata %>% filter(.data$Type=="intron")
    cds_plotdata__INTRON_arrow_pos <-
        seq_len(NROW(cds_plotdata__INTRON)) %>%
        lapply(\(i_tmp){
            df_tmp <-  cds_plotdata__INTRON[i_tmp,]
            arrow_pos <- seq(
                as.integer(df_tmp$Start/arrow_spacing_bp)+1,
                as.integer(df_tmp$End/arrow_spacing_bp)
            )*arrow_spacing_bp

            cbind(
                df_tmp,
                arrow_pos = arrow_pos
            )
        }) %>%
        rbindlist


    arrow_style <- arrow(type = "open", length = unit(0.05, "inch"))

    Gene_levels <- levels(cds_plotdata$Gene)
    if(is.null(Gene_levels)){
        Gene_levels <- sort(unique(cds_plotdata$Gene))
    }



    col_pal_fin <- make_col_pal_fin_gene(col_pal_gene=col_pal,Gene_levels=Gene_levels)





    annot_gene <-
        ggplot(cds_plotdata) +
        aes(col=.data$Gene,fill=.data$Gene) +
        geom_rect(
            data = cds_plotdata__CDS,
            mapping = aes(xmin = .data$Start,xmax=.data$End,ymin=.data$Y_Level-annot_half_height,ymax=.data$Y_Level+annot_half_height)
        )
    if(NROW(cds_plotdata__INTRON)!=0){
        annot_gene <-
            annot_gene +
            geom_linerange(
                data = cds_plotdata__INTRON,
                mapping = aes(xmin = .data$Start,xmax=.data$End,y=.data$Y_Level),
                linetype="dashed"
            ) +
            geom_segment(
                data = cds_plotdata__INTRON_arrow_pos,
                aes(x = .data$arrow_pos, y = .data$Y_Level, xend = .data$arrow_pos + 0.01*ifelse(.data$Strand=="+",1,-1), yend = .data$Y_Level),  # Calculate arrow direction
                arrow = arrow_style, size = 0.5
            )
    }

    annot_gene <-
        annot_gene +
        geom_text(
            data = cds_plotdata__CDS %>% group_by(.data$Gene) %>% slice_min(.data$Start),
            mapping = aes(x=.data$Start,y=.data$Y_Level+annot_half_height,label=.data$Gene),col="black",hjust=0,vjust=-0.3
        ) +
        ylim(1-annot_half_height,max(cds_plotdata$Y_Level)*(1.25))


    annot_gene_fin <-
        annot_gene +
        theme_void() +
        theme(legend.position = "none")

    annot_gene_fin <-
        annot_gene_fin +
        scale_color_manual(values = col_pal_fin) +
        scale_fill_manual(values = col_pal_fin)

    return(annot_gene_fin)
}



get_gene_anno_plot <- memoise(get_gene_anno_plot_ori)




#' Get a pile-up plot, internal function
#'
#' @param x target pile-up vector to plot
#' @param target_cn_table final copy number table(default : NULL)
#' @param baseline the copy number state index you want to set as baseline
#' @param col_cn_baseline color for baseline
#' @param col_pal_cn color palette for non-baseline copy number states
#' @param scale_plot_yaxis Determine whether to scale the plot so that the lower bound of the y-axis is set to the lesser value between 0 and the minimum data value. (Default : TRUE)
#'
#' @importFrom memoise memoise
#' @importFrom grDevices axisTicks axisTicks col2rgb colorRampPalette palette rainbow rgb
#' @importFrom graphics abline axis box image layout lines mtext par plot.new points polygon segments text title
#' @importFrom stats as.dendrogram cutree density dist hclust median quantile sd mad
#' @importFrom utils read.csv tail globalVariables combn
#' @return ggplot object of pile-up plot
#'
#' @noRd
plot_pileUp <-
    function(   x,target_cn_table=NULL,baseline=1,col_cn_baseline = "#708C98",col_pal_cn = col_yarrr_info2  [-5],
                scale_plot_yaxis=TRUE){

        baseline_target <- baseline
        plotdata <-
            data.frame(value = x) %>% mutate(pos=seq_len(n()))

        nonbase_states <-
            target_cn_table$state[target_cn_table$state!=baseline_target]

        nonbase_states_us <- unique(sort(nonbase_states))
        col_pal_cn_fin <- col_pal_cn
        if(length(col_pal_cn)>length(nonbase_states_us)){
            col_pal_cn_fin <- col_pal_cn[seq_len(length(nonbase_states_us))]
        }

        if(length(nonbase_states_us)>0){
            col_map_fin <- structure(
                c(col_cn_baseline,col_pal_cn_fin),
                names =
                    c("baseline",glue("s{nonbase_states_us}"))
            )

            target_cn_table <-
                target_cn_table %>%
                mutate(
                    CN_state_plotting =
                        ifelse(.data$state==baseline_target,"baseline",glue("s{.data$state}")) %>%
                        factor(levels=c("baseline",glue("s{sort(unique(.data$state[.data$state!=baseline_target]))}")))
                )

        }else{
            col_map_fin <- structure(
                c(col_cn_baseline),
                names =
                    c("baseline")
            )
            target_cn_table <-
                target_cn_table %>%
                mutate(
                    CN_state_plotting = factor("baseline",levels="baseline")
                )
        }




        gg_line <-
            ggplot() +
            geom_line(data=plotdata,aes(x=.data$pos,y=.data$value),linewidth=0.5,col="gray40")

        if(scale_plot_yaxis){
            gg_line <-
                gg_line + ylim(min(0,min(x)),max(x))
        }



        if(!is.null(target_cn_table)){
            gg_line <-
                gg_line +
                geom_linerange(data = target_cn_table,aes(  xmin=.data$begin,xmax=.data$end,y=.data$mu,group=.data$seg,
                                                            col=.data$CN_state_plotting
                ))  +
                geom_rect(
                    data = target_cn_table,aes( xmin=.data$begin,xmax=.data$end,ymin=.data$mu-.data$sd,ymax=.data$mu+.data$sd,group=.data$seg,
                                                fill=.data$CN_state_plotting),
                    alpha=0.5
                ) +
                scale_color_manual(values = col_map_fin) +
                scale_fill_manual(values = col_map_fin)  +
                labs(color="CN States",fill="CN States")
        }

        return(gg_line)

    }


#' Get a list of pile-up plots over multiple samples
#'
#' @param result analysis result
#' @param X_raw input raw depth matrix
#' @param target_indices sample indices to plot
#' @param plot_target target data type to plot (Default : `"x"`)
#' @param gff3_fn gene annotation file name
#' @param baseline the state index to set as baseline (Default : `1`)
#' @param annot_margin minimum of margin between gene annotations allowed. As a fraction of plotting area. (Default : `0.01`)
#' @param arrow_spacing gene annotation arrow spacing. As a fraction of plotting area. (Default : `0.05`)
#' @param geme_name_space the height of white space reserved for gene names in the annotation. (Default : `0.5`)
#' @param col_pal gene color palette
#' @param col_cn_baseline color for baseline (Default : `"#708C98"`)
#' @param col_pal_cn color palette for non-baseline copy number states
#' @param exclude_genes name of genes to exclude from the annotation track (Default : NULL)
#' @param annot_plot_ratio ratio of the annotation plot under the pileup plot
#'
#' @return a list of pile-up ggplot object
#' @export
#' @importFrom patchwork plot_layout
#' @examples
#'
#' # gff3 gene model file
#' package_name <- "ELViS"
#' gff3_fn <- system.file("extdata","HPV16REF_PaVE.gff",package = package_name)
#'
#' # loading precalculated depth matrix
#' data(mtrx_samtools_reticulate)
#'
#' # threshold
#' th <- 50
#'
#' # filtered matrix
#' base_resol_depth <- filt_samples(mtrx_samtools_reticulate,th=th,smry_fun=max)
#'
#' data(ELViS_toy_run_result)
#' result <- ELViS_toy_run_result
#'
#' # get line plots for shape-change samples
#' gg_lst_x <-
#' plot_pileUp_multisample(
#'   result = result,
#'   X_raw = base_resol_depth,
#'   plot_target = "x",
#'   gff3 = gff3_fn,
#'   baseline=1,
#'   exclude_genes = c("E6*","E1^E4","E8^E2"),
#'   target_indices = result$final_call$cnv_samples[seq_len(3)]
#'   )
#'
#' gg_lst_x[[1]]
#'
#'
#'
plot_pileUp_multisample <- function(
        result,
        X_raw,
        target_indices = NULL,
        plot_target = "x",
        gff3_fn,
        baseline=1,
        annot_margin = 0.01,
        arrow_spacing = 0.05,
        geme_name_space = 0.5,
        col_pal = col_yarrr_info2,
        col_cn_baseline = "#708C98",
        col_pal_cn = col_yarrr_info2  [-5],
        exclude_genes = NULL,
        annot_plot_ratio = 0.3
){

    stopifnot_ELViS_result(result)
    stopifnot_mtrx_or_df(X_raw)
    stopifnot_character1(plot_target)
    stopifnot_character1(gff3_fn)
    stopifnot_baseline(baseline,n_samples = NCOL(X_raw))
    stopifnot_numeric1(annot_margin)
    stopifnot_numeric1(arrow_spacing)
    stopifnot_numeric1(geme_name_space)
    stopifnot_character_ge1(col_pal)
    stopifnot_character1(col_cn_baseline)
    stopifnot_character_ge1(col_pal_cn)
    if(!is.null(exclude_genes)) stopifnot_character_ge1(exclude_genes)
    stopifnot_numeric1(annot_plot_ratio)

    baseline <- make_baseline_vec(baseline,L=NCOL(X_raw))



    if(is.null(target_indices)){
        target_indices <- result$final_call$cnv_samples
    }
    stopifnot_numeric_ge1(target_indices)

    # plot data setting
    if(plot_target=="x"){
        matrix_type<-"Raw Depth"
        mtrx_for_plotting <- X_raw
        cn_for_plotting <-
            result$final_output %>%
            transmute(
                sampleID=.data$sampleID,
                seg=.data$seg,
                begin=.data$begin,
                end=.data$end,
                mu = .data$mu.x,
                sd = .data$sd.x,
                baseline=.data$baseline,
                state=.data$state,
                cn=.data$cn
            )
    }else if(plot_target=="y"){
        matrix_type<-"Normalized Depth"
        mtrx_for_plotting <- result$new_Y_p2
        cn_for_plotting <-
            result$final_output %>%
            transmute(
                sampleID=.data$sampleID,
                seg=.data$seg,
                begin=.data$begin,
                end=.data$end,
                mu = .data$mu.y,
                sd = .data$sd.y,
                baseline=.data$baseline,
                state=.data$state,
                cn=.data$cn
            )
    }else if(plot_target=="z"){
        matrix_type<-"Robust Z-score"
        mtrx_for_plotting <-
            Z <- t(apply(result$new_Y_p2,1,function(x) pd_rate_hy(x,qrsc=TRUE)))
        cn_for_plotting <-
            result$final_output %>%
            transmute(
                sampleID=.data$sampleID,
                seg=.data$seg,
                begin=.data$begin,
                end=.data$end,
                mu = .data$mu.z,
                sd = .data$sd.z,
                baseline=.data$baseline,
                state=.data$state,
                cn=.data$cn
            )
    }


    plot_gene_anno <- FALSE
    if(!is.null(gff3_fn)){

        plot_gene_anno <- TRUE
        annot_gene <- get_gene_anno_plot(
            gff3_fn = gff3_fn,
            space_length = NROW(X_raw),
            annot_margin = annot_margin,
            arrow_spacing = arrow_spacing,
            geme_name_space = geme_name_space,
            col_pal = col_pal,
            exclude_genes = exclude_genes
        )

        annot_gene_fin <-
            annot_gene +
            xlim(0,NROW(mtrx_for_plotting)+1) +
            theme(plot.margin = unit(c(0, 5.5, 5.5, 5.5), "pt"))
    }



    YLIM <-
        c(
            min(0,mtrx_for_plotting[,target_indices]),
            max(mtrx_for_plotting[,target_indices])
        )

    output <-
        target_indices %>%
        lapply(\(sample_idx){
            x <- mtrx_for_plotting[,sample_idx]
            sample_ID <- colnames(X_raw)[sample_idx]
            target_cn_table <- cn_for_plotting %>% filter(.data$sampleID==sample_ID)


            gg_line <-
                plot_pileUp(
                    x=x,target_cn_table=target_cn_table,baseline=baseline[sample_idx],
                    col_cn_baseline=col_cn_baseline,col_pal_cn=col_pal_cn,
                    scale_plot_yaxis=FALSE
                )
            gg_line_fin <- gg_line + ylim(YLIM[1],YLIM[2])

            gg_line_fin <- gg_line_fin + xlim(0,NROW(mtrx_for_plotting)+1) +
                theme(plot.margin = unit(c(5.5, 5.5, 1, 5.5), "pt")) +
                ggtitle(paste0(matrix_type," ",sample_ID))

            if(plot_gene_anno){
                return(
                    gg_line_fin / annot_gene_fin + plot_layout(heights = c(1-annot_plot_ratio, annot_plot_ratio))

                )
            }else{
                return(gg_line_fin)
            }

        })

    return(output)

}









#' Plot heatmaps based on simple integrative clustering of multiple matrices
#'
#' @param X_raw Raw depth matrix
#' @param result Run result
#' @param gff3_fn gene annotation file name
#' @param exclude_genes name of genes to exclude from the annotation track (Default : NULL)
#' @param col_pal_gene color palette for gene colors
#' @param col_cn Color scheme for copy number heatmap (Default :`colorRamp2(c(0.5,1,1.5),c(muted("blue"),"white",muted("red")))`)
#' @param col_y Color scheme for normalized read depth(Y) heatmap (Default : `colorRamp2(c(0.5,1,2),c(muted("blue"),"white",muted("red")))`)
#' @param col_z Color scheme for Z-score heatmap (Default : `colorRamp2(c(-4,0,4),c(muted("blue"),"white",muted("red")))`)
#' @param col_x_scaled Color scheme for scaled raw depth(X) heatmap (Default : `"auto"`)
#' @param col_vl Color scheme for positional viral load heatmap (Default : `"auto"`)
#' @param baseline Vector of state numbers to use as baseline for each sample. If it is single integer, then the given state number is used for all samples. (Default : `1`)
#' @param matrices_to_plot Names and orders of the matrices to show as heatmap. Any permutation of `c("CN","Y","Z","X_Scaled","Viral_Load")` of any length is allowed. The vertical orders of stacked heatmaps follows the order of this vector. If set to `"all"`, `c("CN","Y","Z","X_Scaled","Viral_Load")` is used. (Default : `"all"`)
#' @param matrices_integ_cluster Names of the matrices to be used for integrative clustering for column orders. Any combination of `c("CN","Y","Z","X_Scaled","Viral_Load")` of length > 1 is allowed. If the length is less then 2, then it is ignored and the first matrix specified in `matrices_to_plot` argument is used for column ordering. The vertical orders of stacked heatmaps follows the order of this vector. If set to `"all"`, `c("CN","Y","Z","X_Scaled","Viral_Load")` is used. (Default : `"all"`)
#' @param total_aligned_base__host_and_virus Total aligned bases for each sample(i.e. from picard,gatk,qualimap). Used to calculate positional load of viral DNA. Makes sense if regions in host genome are also included in the target panel. Ignored if set to NULL. (Default : NULL)
#' @param return_data_matrices boolean whether to return the data matrices used. (Default : `FALSE`)
#'
#' @return A ComplexHeatmap Heatmap List object vertically stacked
#' @export
#'
#' @importFrom magrittr set_names set_colnames set_rownames %>%
#' @importFrom ComplexHeatmap Heatmap rowAnnotation HeatmapAnnotation %v% column_order
#' @importFrom scales alpha trans_new log_breaks muted viridis_pal hue_pal
#'
#' @examples
#'
#'
#' # gff3 gene model file
#' package_name <- "ELViS"
#' gff3_fn <- system.file("extdata","HPV16REF_PaVE.gff",package = package_name)
#'
#' # loading precalculated depth matrix
#' data(mtrx_samtools_reticulate)
#'
#' # threshold
#' th <- 50
#'
#' # filtered matrix
#' base_resol_depth <- filt_samples(mtrx_samtools_reticulate,th=th,smry_fun=max)
#'
#' # viral load data
#' data(total_aligned_base__host_and_virus)
#' viral_load <- (10^6)*(apply(base_resol_depth,2,\(x) sum(x)) )/total_aligned_base__host_and_virus
#'
#' # load ELViS run result
#' data(ELViS_toy_run_result)
#' result <- ELViS_toy_run_result
#'
#' # genes to exclude from plotting
#' exclude_genes <- c("E6*","E1^E4","E8^E2")
#'
#' # heatmap based on integrative clustering
#' integ_ht_result <- integrative_heatmap(
#'   X_raw = base_resol_depth,
#'   result = result,
#'   gff3_fn = gff3_fn,
#'   exclude_genes = exclude_genes,
#'   baseline=1,
#'   total_aligned_base__host_and_virus = total_aligned_base__host_and_virus
#'   )
#'
#' integ_ht_result
#'
#'
#'
integrative_heatmap <- function(
        X_raw,
        result,
        gff3_fn,
        exclude_genes,
        col_pal_gene = col_yarrr_info2  ,
        col_cn = colorRamp2(c(0.5,1,1.5),c(muted("blue"),"white",muted("red"))),
        col_y = colorRamp2(c(0.5,1,2),c(muted("blue"),"white",muted("red"))),
        col_z = colorRamp2(c(-4,0,4),c(muted("blue"),"white",muted("red"))),
        col_x_scaled = "auto",
        col_vl = "auto",
        baseline = 1,
        matrices_to_plot = "all",  # the kind and order of the matrix to show as heatmap
        matrices_integ_cluster = "all", # If only single matrix is specified, then it is ignored and the first matrix of `matrices_to_plot` will be used for clustering. If set to `"all"`, all the matrices available are used for integrative clustering.
        total_aligned_base__host_and_virus = NULL,
        return_data_matrices = FALSE
){

    stopifnot_mtrx_or_df(X_raw)
    stopifnot_ELViS_result(result)
    stopifnot_character1(gff3_fn)
    stopifnot_character_ge1(exclude_genes)
    stopifnot_character_ge1(col_pal_gene)
    stopifnot_ComplexHeatmap_col(col_cn)
    stopifnot_ComplexHeatmap_col(col_y)
    stopifnot_ComplexHeatmap_col(col_z)
    stopifnot_ComplexHeatmap_col(col_x_scaled)
    stopifnot_ComplexHeatmap_col(col_vl)
    stopifnot_baseline(baseline,n_samples = NCOL(X_raw))
    stopifnot_matrices_available(matrices_to_plot)
    stopifnot_matrices_available(matrices_integ_cluster)
    if(!is.null(total_aligned_base__host_and_virus)) stopifnot_numeric_ge1(total_aligned_base__host_and_virus)
    stopifnot_logical1(return_data_matrices)


    rnt_gene_name <- get_gene_rnt(
        gff3_fn = gff3_fn,
        space_length = NROW(X_raw),
        exclude_genes = exclude_genes,
        col_pal_gene = col_pal_gene,
        annotation_name_side = "top"
    )
    rnt_rm_gnm <- get_gene_rnt(
        gff3_fn = gff3_fn,
        space_length = NROW(X_raw),
        exclude_genes = exclude_genes,
        col_pal_gene = col_pal_gene,
        show_annotation_name = FALSE
    )



    viral_load <- (10^6)*(apply(X_raw,2,\(x) sum(x)) )/total_aligned_base__host_and_virus


    matrix_all_ori <- c("CN","Y","Z","X_Scaled","Viral_Load")
    if(is.null(total_aligned_base__host_and_virus)){
        matrix_all <- setdiff(matrix_all,"Viral_Load")
    }else{
        matrix_all <- matrix_all_ori
    }


    baseline <- make_baseline_vec(baseline,L=NCOL(X_raw))
    baseline_target <- baseline

    if(matrices_to_plot[1]=="all"){
        matrices_to_plot <- matrix_all
    }else{
        matrices_to_plot <- intersect(matrices_to_plot,matrix_all)
    }
    matrices_to_plot_original <- matrices_to_plot

    if(matrices_integ_cluster[1]=="all"){
        matrices_integ_cluster <- matrix_all
    }else{
        matrices_integ_cluster <- intersect(matrices_integ_cluster,matrix_all)
    }

    matrix_not_shown <- setdiff(matrices_integ_cluster,matrices_to_plot)
    if(length(matrix_not_shown)>0){
        warning(glue("These matrices were used for clustering but will not be shown in the figure : {paste(matrix_not_shown,collapse=', ')}"))
    }

    if(!inherits(col_vl, "function")){
        if( col_vl[1] =="auto"){
            col_vl <- viridis_pal(option="turbo")(11)
        }
    }

    X_scaled <- apply(X_raw,2,function(x) (10^6)*(x/sum(x)))
    X_scaled_log2p1 <- log2(X_scaled+1)

    if(!inherits(col_x_scaled, "function")){
        if(col_x_scaled[1]=="auto"){
            LQ_Smin <- X_scaled_log2p1 %>% apply(2,min) %>% quantile(0.20)
            MED <- X_scaled_log2p1 %>% median
            UQ_Smax <- X_scaled_log2p1 %>% apply(2,max) %>% quantile(0.95)

            cntr <- MED
            col_x_scaled <-
                viridis_pal(option="turbo")(11) %>%
                colorRamp2(
                    c(
                        (LQ_Smin-cntr)*((length(.)%/%2):1)/(length(.)%/%2),
                        0,
                        (UQ_Smax-cntr)*(seq_len(length(.)%/%2))/(length(.)%/%2)
                    ) +cntr
                    ,.)
        }
    }


    list_ht_name <-
        list(
            CN = "CN",
            Y = "Normalized\nRead Depth",
            Z = "Z-Score",
            X_Scaled = "X Scaled\nLog2p1",
            Viral_Load = "Viral Load\nLog2p1"
        )
    list_ht_rowtitle <-
        list(
            CN = "Copy Number",
            Y = "Normalized\nRead Depth(Y)",
            Z = "Z-Score",
            X_Scaled = "X Scaled",
            Viral_Load = "Positional\nViral Load"
        )
    list_col <-
        list(
            CN = col_cn,
            Y = col_y,
            Z = col_z,
            X_Scaled = col_x_scaled,
            Viral_Load = col_vl
        )
    list_rowant <-
        matrix_all_ori %>% {
            tmp <- .
            structure(
                ifelse(
                    tmp %in% matrices_to_plot[1],
                    list(rnt_gene_name),
                    list(rnt_rm_gnm)
                ),
                names=tmp
            )
        }





    list_matrix <- list()
    list_ht <- list()
    if(is.null(total_aligned_base__host_and_virus)){
        # remove positional viral load plot if overall viral load was not given
        matrices_to_plot <- setdiff(matrices_to_plot,"Viral_Load")
        matrices_integ_cluster <- setdiff(matrices_integ_cluster,"Viral_Load")
    }else{
        key <- "Viral_Load"
        list_matrix[[key]] <-
            log2(t(t((10^6)*X_raw)/total_aligned_base__host_and_virus)+1)

        list_ht[[key]] <-
            Heatmap(
                name=list_ht_name[[key]],
                row_title = list_ht_rowtitle[[key]],row_title_rot = 0,
                list_matrix[[key]],
                cluster_rows = FALSE,
                col=list_col[[key]],
                border=TRUE,
                show_column_names = FALSE,
                right_annotation = list_rowant[[key]]
            )
    }





    key <- "CN"
    list_matrix[[key]] <-
        result$final_output %>%
        group_by(.data$id) %>%
        mutate(cn = .data$cn-.data$cn[.data$state==unique(baseline_target[.data$id])][1]+1 ) %>%
        transmute(cn=.data$cn,id=.data$id,begin=.data$begin,end=.data$end) %>%
        apply(1,\(x) data.frame(pos=x[[3]]:x[[4]],sid=x[[2]],value=x[[1]])) %>% rbindlist %>%
        dcast(pos~sid,value.var = "value") %>%
        select(-pos)

    list_ht[[key]] <-
        Heatmap(
            name=list_ht_name[[key]],
            row_title = list_ht_rowtitle[[key]],row_title_rot = 0,
            list_matrix[[key]],
            cluster_rows = FALSE,
            col=list_col[[key]],
            border=TRUE,
            show_column_names = FALSE,
            right_annotation = list_rowant[[key]]
        )


    key <- "Y"
    list_matrix[[key]] <- result$new_Y_p2

    list_ht[[key]] <-
        Heatmap(
            name=list_ht_name[[key]],
            row_title = list_ht_rowtitle[[key]],row_title_rot = 0,
            list_matrix[[key]],
            cluster_rows = FALSE,
            col=list_col[[key]],
            border=TRUE,
            show_column_names = FALSE,
            right_annotation = list_rowant[[key]]
        )

    key <- "Z"
    list_matrix[[key]] <- t(apply(result$new_Y_p2,1,function(x) pd_rate_hy(x,qrsc=TRUE)))

    list_ht[[key]] <-
        Heatmap(
            name=list_ht_name[[key]],
            row_title = list_ht_rowtitle[[key]],row_title_rot = 0,
            list_matrix[[key]],
            cluster_rows = FALSE,
            col=list_col[[key]],
            border=TRUE,
            show_column_names = FALSE,
            right_annotation = list_rowant[[key]]

        )



    key <- "X_Scaled"
    list_matrix[[key]] <- X_scaled_log2p1

    list_ht[[key]] <-
        Heatmap(
            name=list_ht_name[[key]],
            row_title = list_ht_rowtitle[[key]],row_title_rot = 0,
            list_matrix[[key]],
            cluster_rows = FALSE,
            col=list_col[[key]],
            border=TRUE,
            show_column_names = FALSE,
            right_annotation = list_rowant[[key]]
        )


    final_matrices_used_for_clustering <- matrices_to_plot[1]
    if(length(matrices_integ_cluster)>1){
        final_matrices_used_for_clustering <- matrices_integ_cluster
        Top_matrix_name <- matrices_to_plot[1]

        key <- "Integrative"
        matrices_to_plot[1] <- key
        list_matrix[[key]] <- list_matrix[[Top_matrix_name]]
        list_ht_name[[key]] <- list_ht_name[[Top_matrix_name]]
        list_ht_rowtitle[[key]] <- list_ht_rowtitle[[Top_matrix_name]]
        list_col[[key]] <- list_col[[Top_matrix_name]]
        list_rowant[[key]] <- list_rowant[[Top_matrix_name]]


        dist_list <-
            list_matrix[matrices_integ_cluster] %>%
            lapply(\(mtrx){
                # get sample distance matrix
                mtrx %>% as.matrix %>% t %>% dist()
            })

        ltm_quantile <- Reduce(max,dist_list %>% lapply(\(dst) mean(dst<mean(dst)) ))
        q_p <- 0.75
        q_p <- min(1,max(q_p,ltm_quantile*1.1))
        norm_dist_list <- dist_list %>% lapply(\(dst) dst/quantile(dst,q_p) )

        integ_dendro <-
            {Reduce(`+`,norm_dist_list)/length(matrices_integ_cluster)} %>%
            hclust(method="complete") %>%
            as.dendrogram

        list_ht[[key]] <-
            Heatmap(
                name=list_ht_name[[key]],
                row_title = list_ht_rowtitle[[key]],row_title_rot = 0,
                list_matrix[[key]],column_dend_reorder = TRUE,
                cluster_rows = FALSE,
                col=list_col[[key]],
                border=TRUE,
                show_column_names = FALSE,
                cluster_columns = integ_dendro,
                right_annotation = list_rowant[[key]]
            )

    }

    list_ht_reordered <- list_ht[matrices_to_plot]
    htlst <- Reduce(`%v%`,list_ht_reordered)

    output <-
        list(
            Heatmap = htlst,
            matrices_shown = matrices_to_plot_original,
            matrices_clustered = final_matrices_used_for_clustering,
            column_order = column_order(list_ht_reordered[[1]])
        )

    if(return_data_matrices){
        output <-
            list_matrix[names(list_matrix) != "Integrative"] %>%
            lapply(as.matrix) %>%
            lapply(set_colnames,colnames(X_raw))
    }

    return(
        output
    )

}

parse_gff <- function(gff3_fn,exclude_genes){
    txdb <- makeTxDbFromGFF(gff3_fn, format = "gff3")
    all_features <- transcriptsBy(txdb, by = "gene")
    genes <- genes(txdb) %>% sort
    cds <- cdsBy(txdb, by = "gene")


    if(!is.null(exclude_genes)){
        genes <- genes[!(genes$gene_id %in% exclude_genes)]
        cds <- cds[!(names(cds) %in% exclude_genes)]
    }
    return(
        list(
            txdb = txdb,
            all_features = all_features,
            genes = genes,
            cds = cds
        )
    )
}


#' @noRd
make_col_pal_fin_gene <- function(col_pal_gene,Gene_levels){
    # if color palette can cover all genes, use the scale
    if(length(col_pal_gene)>=length(Gene_levels)){
        #
        if(is.null(names(col_pal_gene))||length(setdiff(Gene_levels,names(col_pal_gene)))!=0){ #pass if col palette designate all colors
            col_pal_fin <- structure(
                col_pal_gene[seq_len(length(Gene_levels))],
                names = Gene_levels
            )
        }else{
            col_pal_fin <- col_pal_gene
        }
    }else{
        col_pal_fin <-
            structure(
                hue_pal()(length(Gene_levels))
                ,names = Gene_levels
            )
    }

    return(
            col_pal_fin
    )
}


#' Row annotation that shows positions of CDSs of genes.
#'
#' @param gff3_fn gene annotation file name
#' @param space_length length of target region
#' @param exclude_genes name of genes to exclude from the annotation track (Default : NULL)
#' @param col_pal_gene palette for gene colors
#'
#' @return ComplexHeatmap RowAnnotation
#' @noRd
get_gene_rnt_ori <- function(
        gff3_fn,
        space_length,
        exclude_genes,
        col_pal_gene = col_yarrr_info2  ,
        ...
){
    gff_parsed <- parse_gff(gff3_fn,exclude_genes)
    txdb <- gff_parsed$txdb
    genes <- gff_parsed$gene
    cds <- gff_parsed$cds

    Gene_levels <- unique(genes$gene_id)

    col_pal_fin <- make_col_pal_fin_gene(col_pal_gene=col_pal_gene,Gene_levels=Gene_levels)

    cds_df <-
        cds[Gene_levels] %>%
        lapply(\(gr){

            target_indice <-
                data.frame(
                    s = start(gr),
                    e = end(gr)
                ) %>%
                apply(1,\(x) seq(x[1],x[2])) %>%
                unlist
            as.character(seq_len(space_length) %in% target_indice)
        }) %>%
        as.data.frame() %>%
        set_colnames(Gene_levels)

    rnt <-
        rowAnnotation(
            df = cds_df,
            col =
                structure(
                    lapply(col_pal_fin,\(x) c("TRUE"=x,"FALSE"="ivory")),
                    names = names(col_pal_fin)
                ),
            border=TRUE,
            show_legend = FALSE,
            ...
        )

    return(rnt)

}

get_gene_rnt <- memoise(get_gene_rnt_ori)


make_baseline_vec <- function(baseline,L){
    if(length(baseline)==1){
        baseline <- rep(baseline,L)
    }
    return(baseline)
}



#' Gene Copy Number Heatmap
#'
#' @param X_raw Raw depth matrix
#' @param result Run result
#' @param gff3_fn gene annotation file name
#' @param gene_ref The name of the gene to set as reference for relative gene dosage heatmap
#' @param baseline Vector of state numbers to use as baseline for each sample. If it is single integer, then the given state number is used for all samples. (Default : `1`)
#' @param exclude_genes name of genes to exclude from the annotation track (Default : NULL)
#' @param col_cn relative gene dosage color palette. (Default : `colorRamp2(c(0.5,1,1.5),c(muted("blue"),"white",muted("red")))`)
#' @param heatmap_height heatmap height specified using unit function. (Default : `unit(1.5,"in")`)
#'
#' @return a ComplexHeatmap Heatmap List object
#' @export
#'
#' @examples
#'
#'
#'
#' # gff3 gene model file
#' package_name <- "ELViS"
#' gff3_fn <- system.file("extdata","HPV16REF_PaVE.gff",package = package_name)
#'
#' # loading precalculated depth matrix
#' data(mtrx_samtools_reticulate)
#'
#' # threshold
#' th <- 50
#'
#' # filtered matrix
#' base_resol_depth <- filt_samples(mtrx_samtools_reticulate,th=th,smry_fun=max)
#'
#' # viral load data
#' data(total_aligned_base__host_and_virus)
#' viral_load <- (10^6)*(apply(base_resol_depth,2,\(x) sum(x)) )/total_aligned_base__host_and_virus
#'
#' # load ELViS run result
#' data(ELViS_toy_run_result)
#' result <- ELViS_toy_run_result
#'
#' # genes to exclude from plotting
#' exclude_genes <- c("E6*","E1^E4","E8^E2")
#'
#' # heatmap of gene dosage
#' gene_ref<-"E7"
#'
#' gene_cn <-
#'   gene_cn_heatmaps(
#'     X_raw = base_resol_depth,
#'     result = result,
#'     gff3_fn = gff3_fn,
#'     baseline = 1,
#'     gene_ref = gene_ref,
#'     exclude_genes = exclude_genes
#'     )
#'
#' gene_cn
#'
#'
#'
#'
gene_cn_heatmaps <-
    function(
        X_raw,
        result,
        gff3_fn,
        gene_ref,
        baseline=1,
        exclude_genes,
        col_cn = colorRamp2(c(0.5,1,1.5),c(muted("blue"),"white",muted("red"))),
        heatmap_height = unit(1.5,"in")
    ){


        stopifnot_mtrx_or_df(X_raw)
        stopifnot_ELViS_result(result)
        stopifnot_character1(gff3_fn)
        stopifnot_character1(gene_ref)
        stopifnot_character_ge1(exclude_genes)
        stopifnot_ComplexHeatmap_col(col_cn)
        stopifnot_unit1(heatmap_height)


        baseline <- make_baseline_vec(baseline,L=NCOL(X_raw))
        baseline_target <- baseline

        gff_parsed <- parse_gff(gff3_fn,exclude_genes)
        txdb <- gff_parsed$txdb
        genes <- gff_parsed$gene
        cds <- gff_parsed$cds

        cds_ul <- unlist(cds,use.names = TRUE)
        mcols(cds_ul)$gene_id <-
            rep(names(cds),vapply(cds,length,0))


        cds_ul_sort <- cds_ul %>% sort

        genes_levels <- unique(genes$gene_id)

        chr <- as.character(seqnames(cds_ul_sort))[1]

        CN_info <-
            result$final_output %>%
            group_by(id) %>%
            mutate(cn = .data$cn-.data$cn[.data$state==unique(baseline_target[.data$id])][1]+1 ) %>%
            ungroup() %>%
            group_split(id) %>%
            lapply(\(df){
                gr <-
                    df %>%
                    transmute(chr=chr,start=.data$begin,end=.data$end) %>%
                    makeGRangesFromDataFrame()


                findOverlaps(cds_ul_sort,gr) %>%
                    as.data.frame() %>%
                    transmute(
                        gene_id = cds_ul_sort$gene_id[.data$queryHits] %>% factor(levels=genes_levels),
                        cn = df$cn[.data$subjectHits]
                    ) %>%
                    group_by(.data$gene_id) %>%
                    summarise(Min_CN=min(.data$cn),Max_CN=max(.data$cn)) %>%
                    mutate(id=df$id[1])
            }) %>%
            rbindlist


        minCN_mtrx <-
            CN_info %>%
            dcast(gene_id~id,value.var="Min_CN") %>%
            as.data.frame %>%
            set_rownames(.$gene_id) %>%
            select(-matches("^gene_id$")) %>%
            set_colnames(colnames(X_raw)) %>%
            as.matrix

        relCN_mtrx <-
            t(t(minCN_mtrx)/minCN_mtrx[gene_ref,])

        ht_gene <-
            Heatmap(
                name="Gene CN",
                row_title = "Intact\nGene CN",
                minCN_mtrx,
                col=col_cn,
                cluster_rows=FALSE,
                show_column_names=FALSE,
                heatmap_height = unit(1.5,"in"),
                row_title_rot = 0,
                border=TRUE
            )

        ht_gene_rel <-
            Heatmap(
                row_title=glue("Dosage Relative\nto {gene_ref}"),
                name = "Relative\nDosage",
                relCN_mtrx,
                col=col_cn,
                cluster_rows=FALSE,
                show_column_names=FALSE,
                heatmap_height = unit(1.5,"in"),
                row_title_rot = 0,
                border=TRUE
            )

        output <-
            list(
                Heatmaps =
                    list(
                        intact_gene_cn = ht_gene,
                        rel_dosage = ht_gene_rel
                    ),
                Matrices =
                    list(
                        intact_gene_cn = ht_gene,
                        rel_dosage = ht_gene_rel
                    )
            )

    }
