
stopifnot_mtrx_or_df <- function(mtrx){
    stopifnot(is(mtrx,"matrix")|is(mtrx,"data.frame"))
}


stopifnot_numeric1 <- function(th){
    stopifnot(is(th,"numeric")&(length(th)==1))
}

stopifnot_numericL <- function(num,L){
    stopifnot_classL(num,class="numeric",L=L)
}

stopifnot_numeric_ge1 <- function(i){
    stopifnot_class_ge1(i,"numeric")
}

stopifnot_integer1 <- function(i){
    stopifnot(is(i,"integer")&(length(i)==1))
}
stopifnot_integer_ge1 <- function(i){
    stopifnot_class_ge1(i,"integer")
}

stopifnot_character1 <- function(tx){
    stopifnot(is(tx,"character")&(length(tx)==1))
}

stopifnot_character_ge1 <- function(tx){
    stopifnot_class_ge1(tx,"character")
}

stopifnot_characterL <- function(tx,L){
    stopifnot_classL(tx,class="character",L=L)
}

stopifnot_listL <- function(lst,L){
    stopifnot_classL(lst,class="list",L=L)
}

stopifnot_list_ge1 <- function(lst,L){
    stopifnot_class_ge1(lst,class="list")
}

stopifnot_classL <- function(input,class,L){
    stopifnot_numeric1(L)
    stopifnot_character1(class)
    stopifnot(is(input,class)&(length(input)==L))
}


stopifnot_class_ge1 <- function(input,class){
    stopifnot_character1(class)
    stopifnot(is(input,class)&(length(input)>0))
}

stopifnot_logical1 <- function(bool){
    stopifnot(is(bool,"logical")&(length(bool)==1))
}

stopifnot_unit1 <- function(u){
    stopifnot(is(u,"unit")&(length(u)==1))
}

stopifnot_probs <- function(probs){
    stopifnot(
        is(probs,"numeric")&
            (length(probs)==1)&
            (0<=probs)&(probs<=1))
}

stopifnot_percent <- function(percent){
    stopifnot(
        is(percent,"numeric")&
            (length(percent)==1)&
            (0<=percent)&(percent<=100))
}


stopifnot_coordinate1 <- function(coord){
    stopifnot_character1(coord)
    stopifnot(grepl("[^:-]+:[0-9,]+-[0-9,]+",coord))
}

stopifnot_baseline <- function(baseline,n_samples){
    stopifnot_numeric1(n_samples)
    stopifnot_numeric_ge1(baseline)
    if(length(baseline)>1){
        stopifnot(length(baseline)==n_samples)
    }
}

stopifnot_win <- function(win){
    stopifnot(is(win,"matrix")&(NCOL(win)==2))
}


stopifnot_normalized_data <- function(normalized_data){
    stopifnot(is(normalized_data,"list"))
    stopifnot(all(c("sample_Ids","Y","Z") %in% names(normalized_data)))
    stopifnot_mtrx_or_df(normalized_data$Y)
    stopifnot_mtrx_or_df(normalized_data$Z)
    stopifnot_characterL(normalized_data$sample_Ids,L=NCOL(normalized_data$Y))
}

stopifnot_refupate_data <- function(refupate_data){
    stopifnot(is(refupate_data,"list"))
    stopifnot(all(c("sample_Ids","Y_recentered","Z_recentered") %in% names(refupate_data)))
    stopifnot_mtrx_or_df(refupate_data$Y_recentered)
    stopifnot_mtrx_or_df(refupate_data$Z_recentered)
    stopifnot_characterL(refupate_data$sample_Ids,L=NCOL(refupate_data$Y_recentered))
    stopifnot(all(c("segment.K_initial") %in% names(refupate_data)))
    stopifnot_numericL(refupate_data$segment.K_initial,L=NCOL(refupate_data$Y_recentered))
}
stopifnot_rescued_data <- function(rescued_data){
    stopifnot_normalized_data(rescued_data)
    stopifnot(all(c("segment.K","clust.list") %in% names(rescued_data)))
    stopifnot_numericL(rescued_data$segment.K,L=NCOL(rescued_data$Y))
    stopifnot_listL(rescued_data$clust.list,L=NCOL(rescued_data$Y))
    stopifnot(is(rescued_data$clust.list[[1]],"segmentation")|is.null(rescued_data$clust.list[[1]]))
}


stopifnot_segment.table <- function(segment.table){
    stopifnot(is(segment.table,"data.frame"))
    stopifnot(all(colnames(segment.table)[seq_len(3)] == c("seg","begin","end")) )
    if(!is_coordinate_only(segment.table)){
        stopifnot(all(colnames(segment.table)[4:12] ==
                        c("mu.x", "sd.x", "mu.y", "sd.y", "mu.z", "sd.z", "cohort.y", "baseline", "state")) )
    }
}

is_coordinate_only <- function(segment.table) {
    return(NCOL(segment.table)<=3)
}

stopifnot_ELViS_result <- function(result){
    stopifnot(is(result,"list"))
    stopifnot( all(c("is_reduced_output", "final_output", "final_call","new_Y_p2" ) %in% names(result) ))
    stopifnot_logical1(result$is_reduced_output)
    stopifnot(is(result$final_call,"list"))
    stopifnot( all(c("segmented_samples", "cnv_samples" ) %in% names(result$final_call) ))
    stopifnot_numeric_ge1(result$final_call$segmented_samples)
    stopifnot_numeric_ge1(result$final_call$cnv_samples)
    stopifnot_mtrx_or_df(result$new_Y_p2)
    stopifnot(is(result$final_output,"data.frame"))
    stopifnot(all(c("id","sampleID","seg","begin","end","mu.x", "sd.x", "mu.y", "sd.y", "mu.z", "sd.z", "cohort.y", "baseline", "state","cn") %in% colnames(result$final_output)))
}


stopifnot_ComplexHeatmap_col <- function(col){
    stopifnot(
        is.null(col)|is.character(col)|is.function(col)
    )
}


stopifnot_matrices_available <- function(matrix_names){
    stopifnot(
        any(
            matrix_names %in%
                c("all","CN","Y","Z","X_Scaled","Viral_Load")
        )
    )
}
