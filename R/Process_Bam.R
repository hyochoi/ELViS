

#' Check if conda installed via reticulate
#' @noRd
check_and_install_conda <- function() {
  # Load the reticulate package

  # Check if Conda is installed
  conda_path <- tryCatch(reticulate::conda_binary(), error = function(e) NULL)

  if (is.null(conda_path)) {
    message("Conda is not installed. Installing Miniconda...")
    reticulate::install_miniconda()  # Install Miniconda
    conda_path <- reticulate::conda_binary()
  } else {
    message("Conda is already installed at: ", conda_path)
  }

  return(invisible(conda_path))
}


#' Generate a read depth matrix of positions x samples from input BAM files list
#'
#' @param bam_files Vector containing bam file names in character
#' @param mode Mode of read depth calculation. Either of `c("samtools_reticulate","samtools_custom","Rsamtools")` are acceptable. If run on Windows OS, it will coerced to `"Rsamtools"` (Default : `"samtools_reticulate"`)
#' @param target_virus_name The name of the target virus. This should be equal to the name of the sequence in the FASTA file reads are aligned to
#' @param N_cores Number of cores to use for parallel processing (Default : min(10,available cores))
#' @param min_mapq Minimum MAPQ. (Default : 30)
#' @param min_base_quality Minimum basecall quality score (Default : 0)
#' @param max_depth (Rsamtools) Maximum read depth. (Default : 1e5)
#' @param modules (samtools) Environment modulefile name. (Default : NULL)
#' @param envs (samtools) Environmental variables for samtools. (Default : NULL)
#' @param tmpdir (samtools) Temporary file directory (Default : `"tmpdir"`)
#' @param remove_tmpdir remove temporary directory after processing is done. Setting to TRUE is not recommended if you are running multiple processes. (Default : FALSE)
#' @param samtools (samtools) Absolute path to samtools executable (Default : NULL)
#' @param condaenv (samtools_reticulate) Name of the conda environment in which samtools are installed. If no environment with this name is available, one will be created. (Default : `"env_samtools"`)
#' @param conda (samtools_reticulate) Path to conda executable. Set to "auto" to let the reticulate automatically decide which conda to use. (Default : "auto")
#'
#'
#' @return a matrix of positions x samples containing base-resolution raw read depth
#' @import parallel
#' @rawNamespace import(data.table, except=c(between,first,last,shift,yearmon,yearqtr))
#' @import reticulate
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
#' N_cores <- 2L
#'
#' # get read depth matrix
#' tmpdir<-"./tmpdir"
#' dir.create(tmpdir,recursive = TRUE)
#'
#' mtrx_samtools_reticulate <-
#'  get_depth_matrix(
#'   bam_files = bam_files,target_virus_name = target_virus_name
#'  ,mode = "samtools_reticulate"
#'  ,N_cores = N_cores
#'  ,min_mapq = 30
#'  ,tmpdir="tmpdir"
#'  ,condaenv = "env_samtools"
#'  ,conda = "auto")
#'
#'  unlink(tmpdir,recursive=TRUE)
#'
#'
get_depth_matrix <-
  function(
    bam_files,mode="samtools_reticulate"
    ,target_virus_name
    #common options
    ,N_cores = detectCores()
    ,min_mapq=30
    ,min_base_quality=0
    #Rsamtools specific options
    ,max_depth = 1e5
    #samtools specific options - custom and reticulate
    ,modules=NULL
    ,envs=NULL
    ,tmpdir="tmpdir"
    ,remove_tmpdir = FALSE
    ,samtools=NULL # absolute path to samtools
    #reticulate specific options
    ,condaenv = "env_samtools"
    ,conda = "auto"
  ){

    # custom samtools not available for windows
    os_name <- Sys.info()["sysname"]
    if( os_name == "Windows" ){
      if(mode != "Rsamtools"){
        warning(glue("mode={mode} is not supported for {os_name}. Changing mode to Rsamtools..."))
        mode <- "Rsamtools"
      }
    }


    if(mode == "samtools_reticulate"){

      # check if reticulate is installed and install if not
      if (!requireNamespace("reticulate", quietly = TRUE)) {
        stop("R Package 'reticulate' does not exist. Please install it by executing following command.\n\nutils::install.packages('reticulate')")
      }

      # Check if conda was installed and install if not
      if(conda == "auto"){
        conda <- check_and_install_conda()
      }

      # Check if environment exist and create one if not
      conda_path <- reticulate::conda_binary(conda = conda)
      env_dir <- paste(dirname(dirname(conda_path)),"envs",condaenv,sep="/")

      if(!reticulate::condaenv_exists(envname = condaenv, conda = conda)){
        reticulate::conda_create(envname = condaenv,channel = c("conda-forge","bioconda"),packages = "samtools",conda = conda)
      }else if(!dir.exists(env_dir)){
        reticulate::conda_create(envname = condaenv,channel = c("conda-forge","bioconda"),packages = "samtools",conda = conda)
      }

      envs <- c(
        PATH = paste(env_dir,"bin",sep="/"),
        LD_LIBRARY_PATH = paste(env_dir,"lib",sep="/")
      )

    }

    if(mode == "Rsamtools"){

      # check if Rsamtools is installed and install if not
      if (!requireNamespace("Rsamtools", quietly = TRUE)) {
        stop("R Package 'reticulate' does not exist. Please install it by executing following command.

if (!requireNamespace('BiocManager', quietly = TRUE))
    utils::install.packages('BiocManager')

BiocManager::install('Rsamtools')")
      }
      # require("Rsamtools")


      out_mtrx <-
        get_depth_matrix_Rsamtools(
          bam_files = bam_files,
          target_virus_name = target_virus_name,N_cores = N_cores,max_depth = max_depth,min_mapq=min_mapq
          ,min_base_quality=min_base_quality
        )
    }else if(mode %in% c("samtools_custom","samtools_reticulate")){
      dir.create(tmpdir)
      out_mtrx <-
        get_depth_matrix_samtools(
          bam_files=bam_files
          ,target_virus_name=target_virus_name
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


    remove_tmpdir

    return(out_mtrx)
  }


#' @noRd
get_depth_matrix_Rsamtools <-
  function(
    bam_files,
    target_virus_name,N_cores = detectCores(),max_depth = 1e5,min_mapq=30,
    min_base_quality=0
  ){

    # viral genome size
    bam <- Rsamtools::BamFile(bam_files[1])
    bam_header <- Rsamtools::scanBamHeader(bam)
    virus_genome_size <- bam_header$targets[[target_virus_name]]

    target_grng <-
      data.frame(chr=target_virus_name,start=1,end=virus_genome_size) %>%
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
      mclapply(X = bam_files,
               mc.cores = N_cores,
               FUN = get_depth_Rsamtools,
               scanBamParam=paramScanBam,pileupParam = paramPileup) %>%
      do.call(cbind, .)
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
  for(s in all_strings){  sanity_check(s)  }


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
    script <- glue("{script}
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


    script <- glue("{script}
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
    ,target_virus_name
    #common options
    ,N_cores = min(10,detectCores()),min_mapq=30,min_base_quality=0
    #samtools specific options
    ,modules=NULL,envs=NULL
    ,tmpdir="tmpdir"
    ,samtools=NULL
  ){

    bash_script_base <-
      get_bash_script_base(modules=modules,envs=envs)

    depth_mtrx <-
      mclapply(mc.cores = N_cores,
               X = seq_along(bam_files)
               ,FUN = get_depth_samtools
               ,bam_files = bam_files
               ,target_virus_name=target_virus_name
               ,min_mapq=min_mapq
               ,min_base_quality=min_base_quality
               ,tmpdir=tmpdir
               ,bash_script_base=bash_script_base
               ,samtools=samtools
      ) %>%
      do.call(cbind, .)
    return(depth_mtrx)
  }




#' @import uuid
#' @noRd
get_depth_samtools <-
  function(
    vec_i
    ,bam_files
    ,target_virus_name
    ,min_mapq=30
    ,min_base_quality=0
    ,tmpdir="tmpdir"
    ,bash_script_base
    ,samtools=NULL
  ){

    bam_fn <- bam_files[vec_i]

    depth_bed_fn <- glue("{tmpdir}/{UUIDgenerate()}_{vec_i}.depth.bed")


    run_samtools(
      bash_script_base = bash_script_base,
      command = glue("depth -a -r '{target_virus_name}' --min-MQ {min_mapq} --min-BQ {min_base_quality} -g 256 {bam_fn}"),
      output_name = depth_bed_fn,
      samtools = samtools,
      depth_count_only = TRUE
    )

    depth_bed <- unlist(read.csv(depth_bed_fn,sep = "\t",header=FALSE),use.names = FALSE)
    sanity_check(depth_bed_fn)
    system2(command = "rm",args = depth_bed_fn)

    return(depth_bed)
  }




#' Generate a read depth matrix of positions x samples from input BAM files list
#'
#' @param bam_files Vector containing bam file names in character
#' @param mode Mode of read depth calculation. Either of `c("samtools_reticulate","samtools_custom","Rsamtools")` are acceptable. If run on Windows OS, it will coerced to `"Rsamtools"` (Default : `"samtools_reticulate"`)
#' @param coord Target region coordinate. This should be a string in the form of "chr1:123-456" or "chr1:1,234-5,678,912"
#' @param N_cores Number of cores to use for parallel processing (Default : min(10,available cores))
#' @param min_mapq Minimum MAPQ. (Default : 30)
#' @param min_base_quality Minimum basecall quality score (Default : 0)
#' @param max_depth (Rsamtools) Maximum read depth. (Default : 1e5)
#' @param modules (samtools) Environment modulefile name. (Default : NULL)
#' @param envs (samtools) Environmental variables for samtools. (Default : NULL)
#' @param tmpdir (samtools) Temporary file directory (Default : `"tmpdir"`)
#' @param samtools (samtools) Absolute path to samtools executable (Default : NULL)
#' @param condaenv (samtools_reticulate) Name of the conda environment in which samtools are installed. If no environment with this name is available, one will be created. (Default : `"env_samtools"`)
#' @param conda (samtools_reticulate) Path to conda executable. Set to "auto" to let the reticulate automatically decide which conda to use. (Default : "auto")
#'
#' @return a matrix of positions x samples containing base-resolution raw read depth
#' @import parallel
#' @noRd
get_depth_matrix_gp <-
  function(
    bam_files,mode="samtools_custom"
    ,coord
    #common options
    ,N_cores = detectCores(),min_mapq=30,min_base_quality=0
    #Rsamtools specific options
    ,max_depth = 1e5
    #samtools specific options - custom and reticulate
    ,modules=NULL,envs=NULL
    ,tmpdir="tmpdir"
    ,samtools=NULL # absolute path to samtools
    #reticulate specific options
    ,condaenv = "env_samtools"
    ,conda = "auto"
  ){
    os_name <- Sys.info()["sysname"]
    if( os_name == "Windows" ){
      if(mode != "Rsamtools"){
        warning(glue("mode={mode} is not supported for {os_name}. Changing mode to Rsamtools..."))
        mode <- "Rsamtools"
      }
    }


    if(mode == "samtools_reticulate"){

      if (!requireNamespace("reticulate", quietly = TRUE)) {
        stop("R Package 'reticulate' does not exist. Please install it by executing following command.\n\nutils::install.packages('reticulate')")
      }



      # Check if environment exist and create one if not
      if(!reticulate::condaenv_exists(envname = condaenv, conda = conda)){
        reticulate::conda_create(envname = condaenv,channel = c("conda-forge","bioconda"),packages = "samtools")
      }
      reticulate::use_condaenv(condaenv = condaenv, conda = conda)

      #get paths
      conda_path <- reticulate::conda_binary()
      env_dir <- paste(dirname(dirname(conda_path)),"envs",condaenv,sep="/")
      envs <- c(
        PATH = paste(env_dir,"bin",sep="/"),
        LD_LIBRARY_PATH = paste(env_dir,"lib",sep="/")
      )

    }

    if(mode == "Rsamtools"){

      # check if Rsamtools is installed and install if not
      if (!requireNamespace("Rsamtools", quietly = TRUE)) {
        stop("R Package 'reticulate' does not exist. Please install it by executing following command.

if (!requireNamespace('BiocManager', quietly = TRUE))
    utils::install.packages('BiocManager')

BiocManager::install('Rsamtools')")
      }

      out_mtrx <-
        get_depth_matrix_Rsamtools_gp(
          bam_files = bam_files,
          coord = coord,N_cores = N_cores,max_depth = max_depth,min_mapq=min_mapq
          ,min_base_quality=min_base_quality
        )
    }else if(mode %in% c("samtools_custom","samtools_reticulate")){
      dir.create(tmpdir)
      out_mtrx <-
        get_depth_matrix_samtools_gp(
          bam_files=bam_files
          ,coord=coord
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
get_depth_matrix_Rsamtools_gp <-
  function(
    bam_files,
    coord,N_cores = detectCores(),max_depth = 1e5,min_mapq=30,
    min_base_quality=0
  ){


    coord_lst <-
      coord %>%
      str_replace_all(",","") %>%
      str_split(":|-",simplify = TRUE) %>%
      as.list() %>%
      structure(names=c("chr","start","end")) %>%
      within({start <- as.numeric(start);end <- as.numeric(end)})


    target_grng <-
      data.frame(chr=coord_lst$chr,start=coord_lst$start,end=coord_lst$end) %>%
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
      mclapply(X = bam_files,
               mc.cores = N_cores,
               FUN = get_depth_Rsamtools,
               scanBamParam=paramScanBam,pileupParam = paramPileup) %>%
      do.call(cbind, .)
    return(depth_mtrx)
  }






#' @noRd
get_depth_matrix_samtools_gp <-
  function(
    bam_files
    ,coord
    #common options
    ,N_cores = detectCores(),min_mapq=30,min_base_quality=0
    #samtools specific options
    ,modules=NULL,envs=NULL
    ,tmpdir="tmpdir"
    ,samtools=NULL
  ){

    bash_script_base <-
      get_bash_script_base(modules=modules,envs=envs)


    depth_mtrx <-
      mclapply(X = seq_along(bam_files)
               ,mc.cores = N_cores
               ,FUN = get_depth_samtools_gp
               ,bam_files = bam_files
               ,coord=coord
               ,min_mapq=min_mapq
               ,min_base_quality=min_base_quality
               ,tmpdir=tmpdir
               ,bash_script_base=bash_script_base
               ,samtools=samtools
      ) %>%
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
  coord_lst <- coord_to_lst(coord)

  grng <-
    data.frame(chr=coord_lst$chr,start=coord_lst$start,end=coord_lst$end) %>%
    makeGRangesFromDataFrame()

  return(grng)
}



#' @import uuid
#' @noRd
get_depth_samtools_gp <-
  function(
    vec_i
    ,bam_files
    ,coord
    ,min_mapq=30
    ,min_base_quality=0
    ,tmpdir="tmpdir"
    ,bash_script_base
    ,samtools=NULL
  ){

    coord_lst <-
      coord %>%
      str_replace_all(",","") %>%
      str_split(":|-",simplify = TRUE) %>%
      as.list() %>%
      structure(names=c("chr","start","end")) %>%
      within({start <- as.numeric(start);end <- as.numeric(end)})


    bam_fn <- bam_files[vec_i]

    depth_bed_fn <- glue("{tmpdir}/{UUIDgenerate()}_{vec_i}.depth.bed")


    run_samtools(
      bash_script_base,
      command <- glue("depth -a -r '{coord_lst$chr}:{coord_lst$start}-{coord_lst$end}' --min-MQ {min_mapq} --min-BQ {min_base_quality} -g 256 {bam_fn}"),
      output_name <- depth_bed_fn,
      samtools <- samtools,
      depth_count_only <- TRUE
    )

    depth_bed <- unlist(read.csv(depth_bed_fn,sep = "\t",header=FALSE),use.names = FALSE)
    sanity_check(depth_bed_fn)
    system2(command = "rm",args = depth_bed_fn)

    return(depth_bed)
  }
