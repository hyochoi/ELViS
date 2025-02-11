test_that("misc functions used for plotting", {

    # make_baseline_vec
    expect_equal(make_baseline_vec(1,3),c(1,1,1))

    # make_col_pal_fin_gene
    expect_equal(
        make_col_pal_fin_gene(NULL,c("gene1","gene2","gene3"))
        ,c(gene1="#F8766D",gene2="#00BA38",gene3="#619CFF")
        )


})
