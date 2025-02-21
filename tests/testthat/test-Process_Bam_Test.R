package_name <- "ELViS"
coord_or_target_virus_name <- "gi|333031|lcl|HPV16REF.1|"
is_virus <- TRUE
# get bam file pathes
ext_path <- system.file("extdata",package = package_name)
bam_files <- list.files(ext_path,full.names = TRUE,pattern = "bam$")
# number of threads to use
N_cores <- 1L
# get read depth matrix
tmpdir <- tempdir()

os_name <- Sys.info()["sysname"]



#### Function to be tested ####
if( os_name != "Windows" ){
mtrx_samtools_basilisk <-
    get_depth_matrix(
        bam_files = bam_files,coord_or_target_virus_name = coord_or_target_virus_name,is_virus = is_virus
        ,mode = "samtools_basilisk"
        ,N_cores = N_cores
        ,min_mapq = 30
        ,tmpdir=tempdir()
        ,condaenv = "env_samtools"
    )
}

max_depth <-  1e5
min_mapq <- 30
min_base_quality <- 0

#### Function to be tested ####
mtrx_Rsatools <-
    get_depth_matrix_Rsamtools(
        bam_files = bam_files,
        coord_or_target_virus_name = coord_or_target_virus_name,is_virus = is_virus,N_cores = N_cores,max_depth = max_depth,min_mapq=min_mapq
        ,min_base_quality=min_base_quality
    )





if(is_virus){
    # viral genome size
    bam <- Rsamtools::BamFile(bam_files[1])
    bam_header <- Rsamtools::scanBamHeader(bam)
    virus_genome_size <- bam_header$targets[[coord_or_target_virus_name]]

    chr <- coord_or_target_virus_name;start <- 1;end <- virus_genome_size
}else{
    coord_lst <- coord_to_lst(coord)

    chr <- coord_lst$chr;start <- coord_lst$start;end <- coord_lst$end
}


target_grng <-
    data.frame(chr=chr,start=start,end=end) %>%
    makeGRangesFromDataFrame()
paramScanBam <- Rsamtools::ScanBamParam(which=target_grng)

paramPileup <-
    Rsamtools::PileupParam(
        min_base_quality=min_base_quality,
        max_depth=max_depth,
        min_mapq=min_mapq,
        min_nucleotide_depth=0,
        distinguish_strands=FALSE,
        distinguish_nucleotides=FALSE,
        ignore_query_Ns=FALSE,
        include_deletions=FALSE,
        include_insertions=FALSE,
        left_bins=NULL,
        query_bins=NULL,
        cycle_bins=NULL
    )

#### Function to be tested ####
mtrx_Rsatools_core <-
    get_depth_matrix_Rsamtools_core(
        bam_files,
        N_cores = N_cores,max_depth = max_depth,min_mapq=min_mapq,
        min_base_quality=min_base_quality,
        chr,start,end
    )

#### Function to be tested ####
depth_Rsatools <- get_depth_Rsamtools(bam_files[1],paramScanBam,paramPileup)


if( os_name != "Windows" ){

condaenv <- "env_samtools"
condaenv_samtools_version <- "1.21"
modules = NULL

#### Function to be tested ####
envs <- get_envs_samtools_basilisk(condaenv_samtools_version,condaenv)

bash_script_base <-
    get_bash_script_base(modules=modules,envs=envs)

samtools <- NULL
region <- coord_or_target_virus_name
bam_fn <- bam_files[1]
vec_i <- 1
depth_bed_fn <- glue("{tmpdir}/{UUIDgenerate()}_{vec_i}.depth.bed")

#### Function to be tested ####
samtools_run_status <-
    run_samtools(
        bash_script_base = bash_script_base,
        command = glue("depth -a -r '{region}' --min-MQ {min_mapq} --min-BQ {min_base_quality} -g 256 {bam_fn}"),
        output_name = depth_bed_fn,
        samtools = samtools,
        depth_count_only = TRUE
    )

#### Function to be tested ####
depth_satools <-
    get_depth_samtools(
        vec_i
        ,bam_files
        ,coord_or_target_virus_name
        ,is_virus
        ,min_mapq=min_mapq
        ,min_base_quality=min_base_quality
        ,tmpdir=tmpdir
        ,bash_script_base
        ,samtools=NULL
    )

#### Function to be tested ####
mtrx_satools <-
    get_depth_matrix_samtools(
        bam_files=bam_files
        ,coord_or_target_virus_name=coord_or_target_virus_name
        ,is_virus = is_virus
        #common options
        ,N_cores = N_cores,min_mapq=min_mapq,min_base_quality=min_base_quality
        #samtools specific options
        ,modules=modules,envs=envs
        ,tmpdir=tmpdir
        ,samtools=samtools
    )

}

#### True values ####
{

true_depth <-
    c(57,58,59,61,64)
true_mtrx <-
    matrix(c(57,58,59,61,64,48,49,51,55,59),ncol=2)
true_mtrx_dimnames <-
    structure(
        true_mtrx
        ,dimnames = list(NULL, c("Control_100X_1102.bam","Control_100X_1119.bam"))
    )

if( os_name != "Windows" ){

samtools_env <- BasiliskEnvironment(
    envname=condaenv
    ,pkgname="ELViS"
    ,channels = c("conda-forge","bioconda")
    ,packages=c(glue("samtools=={condaenv_samtools_version}"))
)

env_dir <- obtainEnvironmentPath(samtools_env)

true_env <- c(
    PATH = file.path(env_dir,"bin"),
    LD_LIBRARY_PATH = file.path(env_dir,"lib")
)

}

}

test_that("Bam Processing_Common", {

    # get_depth_matrix_Rsamtools
    testthat::expect_equal(
        mtrx_Rsatools[1:5,1:2],
        true_mtrx
    )

    # get_depth_matrix_Rsamtools_core
    testthat::expect_equal(
        mtrx_Rsatools_core[1:5,1:2],
        true_mtrx
    )

    # get_depth_Rsamtools
    testthat::expect_equal(
        depth_Rsatools[1:5],
        true_depth
    )

})


if( os_name != "Windows" ){

test_that("Bam Processing Linux and MacOS", {

    # get_depth_matrix
    testthat::expect_equal(
        mtrx_samtools_basilisk[1:5,1:2],
        true_mtrx_dimnames
    )

    # get_depth_matrix_samtools
    testthat::expect_equal(
        mtrx_satools[1:5,1:2],
        true_mtrx
    )

    # get_depth_samtools
    testthat::expect_equal(
        depth_satools[1:5],
        true_depth
    )

    # get_envs_samtools_basilisk
    testthat::expect_equal(
        envs,
        true_env
    )


    # run_samtools
    expect_equal(samtools_run_status,0L)


})

}

false_if_error <- function(expr){
    expr1 <- substitute(expr)  # Capture the unevaluated expression
    res <- tryCatch(expr = eval(expr),error = \(e) FALSE)
    return(res)
}


test_that("Misc Functions",{
    # sanity_check
    expect_false(false_if_error(sanity_check("rm -rf")))
    expect_true(false_if_error(sanity_check("#!/bin/bash")))

    # detect_dollar_unusual
    expect_true(detect_dollar_unusual("bash $"))
    expect_false(detect_dollar_unusual("${3}"))

    # detect_unquoted_pipe
    expect_true(detect_unquoted_pipe("abc|d"))
    expect_false(detect_unquoted_pipe("'|'"))

    # get_bash_script_base
    expect_equal(get_bash_script_base(),"#!/bin/bash\nset -e")

    # check_mode_os
    expect_equal(check_mode_os("Rsamtools"),"Rsamtools")

    # coord_to_lst
    expect_equal(
        coord_to_lst("chr1:1,234-5,678,912"),
        list(
            chr = "chr1",
            start = 1234,
            end = 5678912
        )
    )

    # coord_to_grng
    expect_equal(
        coord_to_grng("chr1:1,234-5,678,912"),
        data.frame(chr="chr1",start=1234,end=5678912) %>%
            makeGRangesFromDataFrame()
    )



})


