#' Generate a read depth matrix of positions x samples from input BAM files list
#'
#' @param bam_files Vector containing bam file names in character
#' @param mode Mode of read depth calculation. Either of `c("samtools_basilisk","samtools_custom","Rsamtools")` are acceptable. If run on Windows OS, it will coerced to `"Rsamtools"` (Default : `"samtools_basilisk"`)
#' @param coord_or_target_virus_name The name of the target virus. This should be equal to the name of the sequence in the FASTA file reads are aligned to.
#' @param is_virus logical indicating if the coord_or_target_virus_name is for viral genome(TRUE) or non-viral genome(FALSE) (default : TRUE)
#' @param N_cores Number of cores to use for parallel processing (Default : min(10,available cores))
#' @param min_mapq Minimum MAPQ. (Default : 30)
#' @param min_base_quality Minimum basecall quality score (Default : 0)
#' @param max_depth (Rsamtools) Maximum read depth. (Default : 1e5)
#' @param modules (samtools) Environment modulefile name. (Default : NULL)
#' @param envs (samtools) Environmental variables for samtools. (Default : NULL)
#' @param tmpdir (samtools) Temporary file directory (Default : `tempdir()`)
#' @param samtools (samtools) Absolute path to samtools executable (Default : NULL)
#' @param condaenv (samtools_basilisk) Name of the conda environment in which samtools are installed. If no environment with this name is available, one will be created. (Default : `"env_samtools"`)
#' @param condaenv_samtools_version (samtools_basilisk) The version of samtools to install in the conda environment using basilisk (Default : "1.21")
#'
#'
#' @return a matrix of positions x samples containing base-resolution raw read depth
#' @importFrom  parallel mclapply    detectCores
#' @rawNamespace import(data.table, except=c(between,first,last,shift,yearmon,yearqtr))
#' @importFrom basilisk BasiliskEnvironment obtainEnvironmentPath
#' @export
#'
#' @examples
#'
#' package_name <- "ELViS"
#'
#' # The name of the target virus
#' # in the reference sequence FASTA file used for alignment.
#' # Can be check by samtools view -H input.bam
#' target_virus_name <- "gi|333031|lcl|HPV16REF.1|"
#'
#' # get bam file pathes
#' ext_path <- system.file("extdata",package = package_name)
#' bam_files <- list.files(ext_path,full.names = TRUE,pattern = "bam$")
#'
#' # number of threads to use
#' N_cores <- 1L
#'
#' # get read depth matrix
#' tmpdir <- tempdir()
#'
#' mtrx_samtools_basilisk <-
#'  get_depth_matrix(
#'   bam_files = bam_files,coord_or_target_virus_name = target_virus_name,is_virus = TRUE
#'  ,mode = "samtools_basilisk"
#'  ,N_cores = N_cores
#'  ,min_mapq = 30
#'  ,tmpdir=tempdir()
#'  ,condaenv = "env_samtools"
#'  )
#'
#'
get_depth_matrix <-
    function(
        bam_files,mode="samtools_basilisk"
        ,coord_or_target_virus_name
        ,is_virus = TRUE
        #common options
        ,N_cores = detectCores()
        ,min_mapq=30
        ,min_base_quality=0
        #Rsamtools specific options
        ,max_depth = 1e5
        #samtools specific options - custom and basilisk
        ,modules=NULL
        ,envs=NULL
        ,tmpdir=tempdir()
        ,samtools=NULL # absolute path to samtools
        #basilisk specific options
        ,condaenv = "env_samtools"
        ,condaenv_samtools_version = "1.21"
    ){

        stopifnot_character_ge1(bam_files)
        stopifnot_character1(mode)
        stopifnot_character1(coord_or_target_virus_name)
        stopifnot_logical1(is_virus)
        stopifnot_numeric1(N_cores)
        stopifnot_numeric1(min_mapq)
        stopifnot_numeric1(min_base_quality)
        stopifnot_numeric1(max_depth)
        if(!is.null(modules)) stopifnot_character_ge1(modules)
        if(!is.null(envs)) stopifnot_character_ge1(envs)
        stopifnot_character1(tmpdir)
        if(!is.null(samtools)) stopifnot_character1(samtools)
        stopifnot_character1(condaenv)
        stopifnot_character1(condaenv_samtools_version)



        # custom samtools not available for windows
        mode <- check_mode_os(mode)

        if(mode == "samtools_basilisk"){
            envs <- get_envs_samtools_basilisk(condaenv_samtools_version,condaenv)
        }

        if(mode == "Rsamtools"){

            # check if Rsamtools is installed and install if not
            check_Rsamtools_installation()

            out_mtrx <-
                get_depth_matrix_Rsamtools(
                    bam_files = bam_files,
                    coord_or_target_virus_name = coord_or_target_virus_name,is_virus = is_virus,N_cores = N_cores,max_depth = max_depth,min_mapq=min_mapq
                    ,min_base_quality=min_base_quality
                )
        }else if(mode %in% c("samtools_custom","samtools_basilisk")){
            if(!dir.exists(tmpdir)){ dir.create(tmpdir) }
            out_mtrx <-
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

        colnames(out_mtrx) <- basename(bam_files)
        if(
            grepl("Error",out_mtrx[1,1],ignore.case = TRUE)
        ){
            stop(out_mtrx[1,1])
        }

        return(out_mtrx)
    }


#' @noRd
get_depth_matrix_Rsamtools <-
    function(
        bam_files,
        coord_or_target_virus_name,is_virus,N_cores = detectCores(),max_depth = 1e5,min_mapq=30,
        min_base_quality=0
    ){

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


        depth_mtrx <-
            get_depth_matrix_Rsamtools_core(
                bam_files,
                N_cores=N_cores,max_depth=max_depth,min_mapq=min_mapq,
                min_base_quality=min_base_quality,
                chr=chr,start=start,end=end
            )

        return(depth_mtrx)
    }

#' @noRd
get_depth_Rsamtools <-
    function(bam_fn,scanBamParam,pileupParam){
        Rsamtools::pileup(bam_fn, scanBamParam=scanBamParam,pileupParam = pileupParam)$count
    }


















#' system call sanity check
#'
#' @param s string or partial string that are to be used for system call
#'
#' @return return TRUE if input string is clear and FALSE if not.
#' @noRd
sanity_check <- function(s){

    if(is.null(s)){
        return(TRUE)
    }

    dangerous_patterns <- c(
        "rm\\s+-rf",    # Recursive delete
        "mv\\s+.*\\s+/", # File move
        "[;&><`]",     # Command chaining and shell control
        "sudo",         # Privileged command execution
        "chmod",        # File permission changes
        "chown"         # File ownership changes
    )

    if (any(grepl(paste(dangerous_patterns, collapse = "|"), s))) {
        stop(glue("Dangerous command detected! : \n{s}"))
    }else if(grepl("[$]",s)){
        if(detect_dollar_unusual(s)){
            stop(glue("Dangerous command detected! : \n{s}"))
        }
    }else{

        if(detect_unquoted_pipe(s)){
            stop(glue("Dangerous command detected! : \n{s}"))
        }

    }
    return(TRUE)
}

#' check if dollar signs were improperly used if any
#' @noRd
detect_dollar_unusual <- function(s){
    all_dollars <- str_extract_all(s,"[$]",simplify = TRUE)
    allowed_patterns <- str_extract_all(s,'[$][{][^{}]+[}]|"[$][{][^{}]+[}]"',simplify = TRUE)
    if(length(all_dollars)!=length(allowed_patterns)){
        return(TRUE)
    }else{
        return(FALSE)
    }
}

#' check if there are any unquoted pipes
#' @noRd
detect_unquoted_pipe <- function(s) {
    # Split the text by quotes and retain quoted and unquoted parts
    parts <- strsplit(s, "'")[[1]]

    if( (length(parts)%%2==0)&(!grepl("'$",s)) ){
        return(TRUE)
    }

    # Initialize flag to alternate between quoted and unquoted sections
    is_quoted <- TRUE

    # Iterate over parts to check for dangerous patterns only in unquoted parts
    for (part in parts) {
        # Alternate the flag (quoted/unquoted)
        is_quoted <- !is_quoted

        # If the part is unquoted, check for dangerous patterns
        if (!is_quoted) {
            # Check for [ or ] or | in the unquoted part
            if (grepl("[|]", part)) {
                # Dangerous pattern detected!
                return(TRUE)
            }
        }
    }

    return(FALSE)
}






#' Function to prepare initializing commands for bash script
#'
#' @param module character vector containing the names of modules to load
#' @param env list of environment variables. Names should match the variable name and the elements are the strings to be inputed to the variable.
#' @noRd
get_bash_script_base <- function(modules=NULL,envs=NULL,debug_mode=FALSE){

    all_strings <- c(modules,envs,names(envs))
    all_strings %>% lapply(sanity_check)


    if(debug_mode){
        script <-
            "#!/bin/bash
set -ex"
    }else{
        script <-
            "#!/bin/bash
set -e"
    }



    if(!is.null(modules)){
        module_string <- glue("module load {paste(modules,collapse=' ')}") |> as.character()
        script <-   glue("{script}
                    {module_string}") |> as.character()
    }

    if(!is.null(envs)){
        if(is.null(names(envs))){
            stop("env should be named list. ex) list(PATH=c('/bin','/usr/bin','/home/user/bin'))")
        }


        env_string <-
            names(envs) |>
            vapply(\(env_names){
                env_colon_concat <- paste(envs[[env_names]],collapse = ":")
                glue('export {env_names}="${{{env_names}}}":{env_colon_concat}')
            },""
            ) |>
            paste(collapse="\n")


        script <-   glue("{script}
                    {env_string}") |> as.character()
    }

    sanity_check(script)
    return(script)


}

#' Run samtools using system function. Sanity check is included
#'
#' @param bash_script_base Character. Output of the function get_bash_script_base
#' @param command Character. samtools commands to use
#' @param output_name Character. The name of the output file. Do not save to a file if NULL. (Default : NULL)
#' @param samtools Character. Samtools path. Default samtools in the PATH is used if NULL (Default : NULL)
#'
#' @noRd
run_samtools <- function(
        bash_script_base,
        command,
        output_name=NULL,
        samtools=NULL,
        depth_count_only = TRUE
){

    sanity_check(command)
    sanity_check(samtools)
    sanity_check(output_name)



    if(is.null(samtools)){
        message("The path to samtools not provided.")
        samtools_path <-
            system2(
                command = "which"
                ,args = "samtools"
                ,env = paste0(bash_script_base,"\n")
                ,stdout=TRUE)

        if( !( samtools_path |> attr("status") |> is.null() ) ){
            stop("samtools is not in the PATH. Please provide paths to samtools and required environment variables if necessary.")
        }
        message(glue("Default samtools is used : {samtools_path}"))

        samtools <- samtools_path
    }

    args <- command
    output <- ""

    if(depth_count_only){
        args <- c(args,"| cut -f 3")
    }

    if(!is.null(output_name)){
        output <- output_name
    }


    system2(
        command = samtools
        ,args = args
        ,env = paste0(bash_script_base,"\n")
        ,stdout = output
    )

}



#' @noRd
get_depth_matrix_samtools <-
    function(
        bam_files
        ,coord_or_target_virus_name
        ,is_virus
        #common options
        ,N_cores = min(10,detectCores()),min_mapq=30,min_base_quality=0
        #samtools specific options
        ,modules=NULL,envs=NULL
        ,tmpdir=tempdir()
        ,samtools=NULL
    ){

        bash_script_base <-
            get_bash_script_base(modules=modules,envs=envs)

        depth_mtrx <-
            mclapply(   mc.cores = N_cores,
                        X = seq_along(bam_files)
                        ,FUN = get_depth_samtools
                        ,bam_files = bam_files
                        ,coord_or_target_virus_name = coord_or_target_virus_name
                        ,is_virus = is_virus
                        ,min_mapq=min_mapq
                        ,min_base_quality=min_base_quality
                        ,tmpdir=tmpdir
                        ,bash_script_base=bash_script_base
                        ,samtools=samtools
            ) %>%
            do.call(cbind, .)
        return(depth_mtrx)
    }




#' @importFrom uuid UUIDgenerate
#' @noRd
get_depth_samtools <-
    function(
        vec_i
        ,bam_files
        ,coord_or_target_virus_name
        ,is_virus
        ,min_mapq=30
        ,min_base_quality=0
        ,tmpdir=tempdir()
        ,bash_script_base
        ,samtools=NULL
    ){

        bam_fn <- bam_files[vec_i]

        depth_bed_fn <- glue("{tmpdir}/{UUIDgenerate()}_{vec_i}.depth.bed")

        if(is_virus){
            region <- coord_or_target_virus_name
        }else{
            coord_lst <- coord_to_lst(coord_or_target_virus_name)
            region <- glue("{coord_lst$chr}:{coord_lst$start}-{coord_lst$end}")
        }


        run_samtools(
            bash_script_base = bash_script_base,
            command = glue("depth -a -r '{region}' --min-MQ {min_mapq} --min-BQ {min_base_quality} -g 256 {bam_fn}"),
            output_name = depth_bed_fn,
            samtools = samtools,
            depth_count_only = TRUE
        )

        depth_bed <- unlist(read.csv(depth_bed_fn,sep = "\t",header=FALSE),use.names = FALSE)
        sanity_check(depth_bed_fn)
        system2(command = "rm",args = depth_bed_fn)

        return(depth_bed)
    }

check_mode_os <- function(mode){
    os_name <- Sys.info()["sysname"]

    if( os_name == "Windows" ){
        if(mode != "Rsamtools"){
            warning(glue("mode={mode} is not supported for {os_name}. Changing mode to Rsamtools..."))
            mode <- "Rsamtools"
        }
    }

    if( !(mode %in% c("samtools_basilisk","samtools_custom","Rsamtools")) ){
        stop(glue("mode='{mode}' is not an allowed argument. Available arguments are 'samtools_basilisk','samtools_custom', and 'Rsamtools'"))
    }
    return(mode)
}

get_envs_samtools_basilisk <- function(condaenv_samtools_version,condaenv){
    # check if basilisk is installed and install if not
    if (!requireNamespace("basilisk", quietly = TRUE)) {
        stop("R Package 'basilisk' does not exist. Please install it by following instructions in 'https://www.bioconductor.org/packages/release/bioc/html/basilisk.html'")
    }

    # samtools version sanity check
    if(grepl("[^0-9.]",condaenv_samtools_version)){
        stop("Invalid samtools version number. Please find correct version number refering to 'https://anaconda.org/bioconda/samtools'.")
    }


    # Load samtools conda environment
    samtools_env <- BasiliskEnvironment(
        envname=condaenv
        ,pkgname="ELViS"
        ,channels = c("conda-forge","bioconda")
        ,packages=c(glue("samtools=={condaenv_samtools_version}"))
    )

    env_dir <- obtainEnvironmentPath(samtools_env)

    envs <- c(
        PATH = file.path(env_dir,"bin"),
        LD_LIBRARY_PATH = file.path(env_dir,"lib")
    )

    return(envs)
}

check_Rsamtools_installation <- function(){
    if (!requireNamespace("Rsamtools", quietly = TRUE)) {
        stop("R Package 'Rsamtools' does not exist. Please install it by executing following command.

if (!requireNamespace('BiocManager', quietly = TRUE))
    utils::install.packages('BiocManager')

BiocManager::install('Rsamtools')")
    }
}



get_depth_matrix_Rsamtools_core <-
    function(
        bam_files,
        N_cores = detectCores(),max_depth = 1e5,min_mapq=30,
        min_base_quality=0,
        chr,start,end
    ){

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





        depth_mtrx <-
            mclapply(   X = bam_files,
                        mc.cores = N_cores,
                        FUN = get_depth_Rsamtools,
                        scanBamParam=paramScanBam,pileupParam = paramPileup) %>%
            do.call(cbind, .)
        return(depth_mtrx)
    }





#' Convert coordinate string to list of chr,start and end
#'
#' @param coord string in the form of "chr1:123-456" or "chr1:1,234-5,678,912"
#'
#' @return a list of 3 elements. Chromosome name, start position and end position.
#' @export
#'
#' @examples
#' coord_to_lst("chr1:123-456")
#' coord_to_lst("chr1:1,234-5,678,912")
#'
coord_to_lst <- function(coord){
    stopifnot_coordinate1(coord)
    coord_lst <-
        coord %>%
        str_replace_all(",","") %>%
        str_split(":|-",simplify = TRUE) %>%
        as.list() %>%
        structure(names=c("chr","start","end")) %>%
        within({start <- as.numeric(start);end <- as.numeric(end)})

    return(coord_lst)
}

#' Convert coordinate string to grng object
#'
#' @param coord string in the form of "chr1:123-456" or "chr1:1,234-5,678,912"
#'
#' @return GRanges object corresponding to the input coordinate string
#' @export
#'
#' @examples
#' coord_to_grng("chr1:123-456")
#' coord_to_grng("chr1:1,234-5,678,912")
#'
coord_to_grng <- function(coord){
    stopifnot_coordinate1(coord)
    coord_lst <- coord_to_lst(coord)

    grng <-
        data.frame(chr=coord_lst$chr,start=coord_lst$start,end=coord_lst$end) %>%
        makeGRangesFromDataFrame()

    return(grng)
}


