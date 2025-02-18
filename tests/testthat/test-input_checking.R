expect_true_wrapper_func_stopifnot <- function(expr){
    expr1 <- substitute(expr)  # Capture the unevaluated expression
    res <- tryCatch(expr = is.null(eval(expr)),error = \(e) FALSE)
    return(expect_true(res))
}

clust_obj <-
    segmentation(data.frame(x=1:10,y=10:1),lmin=5,Kmax=2)


segment.table_coordonly =
    data.frame(
        seg = 1:3
        ,begin = c(10,50,100)
        ,end = c(10,50,100)
    )

segment.table_full =
    cbind(
        segment.table_coordonly,
        data.frame(
            mu.x = c(100,200,300)
            ,sd.x = c(10,10,20)
            ,mu.y = c(1,1,1.5)
            ,sd.y = c(0.1,0.1,0.15)
            ,mu.z = c(-0.1,+0.1,+0.3)
            ,sd.z = c(0.1,+0.1,+0.3)
            ,cohort.y = c(0.9,1,1.6)
            ,baseline = c(1,1,0)
            ,state = c(1,1,2)
        )

    )

data(ELViS_toy_run_result)
result = ELViS_toy_run_result


test_that("input_checking", {

    # running without failing itself ensures these functions are working in the desired way

    # stopifnot_mtrx_or_df
    expect_true_wrapper_func_stopifnot(stopifnot_mtrx_or_df(matrix(1:12,nrow = 4)))

    # stopifnot_numeric1,stopifnot_numericL,stopifnot_numeric_ge1
    expect_true_wrapper_func_stopifnot(stopifnot_numeric1(1))
    expect_true_wrapper_func_stopifnot(stopifnot_numericL(1:3,L=3))
    expect_true_wrapper_func_stopifnot(stopifnot_numeric_ge1(1:3))

    # stopifnot_integer1,stopifnot_integer_ge1
    expect_true_wrapper_func_stopifnot(stopifnot_integer1(1L))
    expect_true_wrapper_func_stopifnot(stopifnot_integer_ge1(1L:3L))

    # stopifnot_character1,stopifnot_characterL,stopifnot_character_ge1
    expect_true_wrapper_func_stopifnot(stopifnot_character1("a"))
    expect_true_wrapper_func_stopifnot(stopifnot_characterL(c("c","d"),L=2))
    expect_true_wrapper_func_stopifnot(stopifnot_character_ge1(c("e","f")))

    # stopifnot_listL
    expect_true_wrapper_func_stopifnot(stopifnot_listL(list(1,"2",TRUE),L=3))

    # stopifnot_classL
    expect_true_wrapper_func_stopifnot(stopifnot_classL(list(1,"2",TRUE),class="list",L=3))

    # stopifnot_class_ge1
    expect_true_wrapper_func_stopifnot(stopifnot_class_ge1(list(1,"2"),class="list"))

    # stopifnot_list_ge1
    expect_true_wrapper_func_stopifnot(stopifnot_list_ge1(list(1,"2")))

    # stopifnot_unit1
    expect_true_wrapper_func_stopifnot(stopifnot_unit1(unit(1,unit="in")))

    # stopifnot_logical1
    expect_true_wrapper_func_stopifnot(stopifnot_logical1(TRUE))

    # stopifnot_probs
    expect_true_wrapper_func_stopifnot(stopifnot_probs(0.5))

    # stopifnot_percent
    expect_true_wrapper_func_stopifnot(stopifnot_percent(50))

    # stopifnot_coordinate1
    expect_true_wrapper_func_stopifnot(stopifnot_coordinate1("chr1:123,456-789,123"))

    # stopifnot_baseline
    expect_true_wrapper_func_stopifnot(stopifnot_baseline(rep(1,5),n_samples = 5))

    # stopifnot_win
    expect_true_wrapper_func_stopifnot(stopifnot_win(matrix(1:10,ncol=2)))

    # stopifnot_normalized_data
    expect_true_wrapper_func_stopifnot(stopifnot_normalized_data(
        list(
            sample_Ids = paste0("s",1:5)
            ,Y = matrix(1:20,ncol=5)
            ,Z = t(scale(t(matrix(1:20,ncol=5))))
            )
    ))

    # stopifnot_refupate_data
    expect_true_wrapper_func_stopifnot(stopifnot_refupate_data(
        list(
            sample_Ids = paste0("s",1:5)
            ,Y_recentered = matrix(1:20,ncol=5)
            ,Z_recentered = t(scale(t(matrix(1:20,ncol=5))))
            ,segment.K_initial = c(2,2,3,1,2)
        )
    ))

    # stopifnot_rescued_data
    expect_true_wrapper_func_stopifnot(stopifnot_rescued_data(
        list(
            sample_Ids = paste0("s",1:5)
            ,Y = matrix(1:20,ncol=5)
            ,Z = t(scale(t(matrix(1:20,ncol=5))))
            ,segment.K = c(2,2,3,1,2)
            ,clust.list = rep(list(clust_obj),5)
        )
    ))

    # stopifnot_segment.table
    expect_true_wrapper_func_stopifnot(
        stopifnot_segment.table(segment.table_full)
    )

    # is_coordinate_only
    expect_true(
        is_coordinate_only(segment.table_coordonly)
    )

    # stopifnot_ELViS_result
    expect_true_wrapper_func_stopifnot(
        stopifnot_ELViS_result(result)
    )

    # stopifnot_ComplexHeatmap_col
    expect_true_wrapper_func_stopifnot(
        stopifnot_ComplexHeatmap_col(c("blue","white","red"))
    )

    # stopifnot_matrices_available
    expect_true_wrapper_func_stopifnot(
        stopifnot_matrices_available( c("all","CN","Y","Z","X_Scaled","Viral_Load"))
    )



})
