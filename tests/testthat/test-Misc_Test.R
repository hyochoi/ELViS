precScale <- 1e-10
x <- 1:10
fsp <- fastSplitSample(x)
x2 <- x - fsp$med

N_alt_ori <- 50L
N_ctrl_ori <- 600L

test_that("Misc functions", {

    # is_zero_dec
    expect_true(is_zero_dec(10))
    expect_true(is_zero_dec(10.0000000000001))
    expect_false(is_zero_dec(10.1))

    # get_perf
    expect_equal(
        get_perf(table(c(0,1,0,1,1,0,1),c(0,0,1,1,0,0,0)))
        ,data.frame(
            Sensitivity=0.666666666666667,Specificity=0.25,PPV=0.4,NPV=0.5,Accuracy=0.428571428571429,F1=0.5
            ,row.names = 0
        )
    )

    # get_A_B_Prop_Subcount
    expect_equal(
        get_A_B_Prop_Subcount(A=N_alt_ori,B=N_ctrl_ori,target_prop=0.01)
        ,list(N_alt_ori=6,N_ctrl_ori=594)
    )

    # get_A_B_Prop_Subcount__internal
    expect_equal(
        get_A_B_Prop_Subcount__internal(A=N_alt_ori,A_str = "N_alt_ori",B=N_ctrl_ori,B_str = "N_ctrl_ori",target_prop = 0.01,n_cycle_th=10)
        ,list(N_alt_ori=6,N_ctrl_ori=594)
    )

    # pd_rate_hy
    expect_equal(pd_rate_hy(x),c(-1.21408336705787,-0.944287063267233,-0.674490759476595,-0.404694455685957,-0.134898151895319,0.134898151895319,0.404694455685957,0.674490759476595,0.944287063267233,1.21408336705787))

    # compScales
    expect_equal(
        compScales(1:10),
        list(sa=2.97881408,sb=2.97881408,med=5.5)
    )

    # fastSplitSample
    expect_equal(fsp,
                list(
                    xa = c(0.5, 1.5, 2.5, 3.5, 4.5)
                    ,xb = c(4.5, 3.5, 2.5, 1.5, 0.5)
                    ,med = 5.5
                    )
                )

    # scale1StepM
    expect_equal(
        scale1StepM(x2,precScale),2.97881408
    )

    # rhoHuber
    expect_equal(
        rhoHuber(x2/ (1.4826*median(abs(x2))) ),
        c(0.792683595920693,0.479524644445851,0.24465543083972,0.0880759551022992,0.0097862172335888,0.0097862172335888,0.0880759551022992,0.24465543083972,0.479524644445851,0.792683595920693)
    )

    # yaxis_hy
    expect_equal(
        yaxis_hy(matrix(1:10,ncol=2))
        ,c(0.982,10.018)
    )

})


