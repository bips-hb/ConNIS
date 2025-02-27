#' Results based on subsamples
#'
#' 'subsampleResults' draws *m* subsamples of the observerd IS and calculates
#' for each subsample the results for a user selected analysis method. Will be
#' used by `instabilities`.
#'
#' @importFrom parallel parLapply
#'
#' @param method Method that will be applied to subsamples. Available are
#' `"Binomial"`, `"ConNIS"`, `"Geometric"` and `"Tn5Gaps"`.
#' @param ins.positions Numeric vector of position of observed insertions sites.
#' @param gene.names The names of the genes.
#' @param gene.starts Starting position within the genome of each gene.
#' @param gene.stops Ending position within the genome of each gene.
#' @param num.ins.per.gene Number of unique insertion sites within each gene.
#' @param genome.length Length of the genome.
#' @param weightings A sequence of weight values that are applied to the
#' selected method.
#' @param m Number of sub samples.
#' @param d Proportion of the original IS used as sub samples.
#' @param use.parallel Should the calcutltions be done in parallel?
#' @param parallelization.type Which method should be used for the
#' parallelization? Available are `"mclapply"` and `"parLapply"` of the
#' `"parallel"` package.
#' @param numCores Number of cores to juse for parallel calculation. Should not
#' exit the number of available (logical) cores of the system.
#' @param cluster.type If `"parLapply"` is used as parallelization.type, a
#' cluster type like `"PSOCK"` needs tp be specified. Check which types are
#' supported on your system. See als `?makeCluster` of the `parallel` package.
#' @param seed Seed for subsampling. If `NULL` no seed is set and results might
#' not be reproducable.
#' @param rng Which random number generator (RNG) should be used. If NULL, the
#' default RNG of the system is used. NOTE: for parallelisation.
#' `"L'Ecuyer-CMRG"` should be used. For details see the manual of the
#' `parallel` package.
#'
#' @returns A list of `tibble`s. Each `tibble` is based on one subsample and
#' contains the results of each genes of each tuning/weight value.
#'
#' @examples
#' # generate random insertion sites
#' set.seed(1)
#' random_is <- sort(sample(1:10000, 2000))
#'
#' # generate dummy gene names
#' genes <- paste("gene_", 1:30)
#'
#' # generate radom start points and stop points of the dummy genes
#' set.seed(2)
#' x <- sort(sample(1:10000, 60))
#' starts <- x[seq(1,60, 2)]
#' stops <- x[seq(2,60, 2)]
#'
#' genome_length <- 10000
#'
#' # run supsumple procedure with ConNIS
#' subsampleResults(method="ConNIS",
#'                   ins.positions = random_is,
#'                   gene.names = genes,
#'                   gene.starts = starts,
#'                   gene.stops = stops,
#'                   genome.length = genome_length,
#'                   weightings = c(0.3, 0.6),
#'                   m = 2,
#'                   d = 0.5,
#'                   use.parallel = FALSE,
#'                   seed = 1)
#'
#' @export

subsampleResults <-
  function(method=c("Binomial", "ConNIS", "Geometric", "Tn5Gaps"),
           ins.positions,
           gene.names,
           gene.starts,
           gene.stops,
           num.ins.per.gene = NULL,
           genome.length,
           weightings = 1,
           m = 100,
           d = 0.5,
           use.parallel=FALSE,
           parallelization.type=c("mclapply", "parLapply"),
           numCores=3,
           cluster.type= NULL,
           seed = NULL,
           rng = NULL){

    system_RNG <-RNGkind()[1]

    # catch errors:
    if(d<=0 | d >= 1 | length(d)!=1){
      stop("Error: d has to be a single value >0 and <1.")
    }

    # print warnings if RNG is not set appropriately
    if(is.null(rng) & use.parallel==FALSE){

      warning(paste("You did not specify an RNG. Will use your defaul one: ", system_RNG, sep=""))

    }else if(is.null(seed) | is.null(rng)){

      warning("You set either no RNG or no seed: results might not be reproducable.")

    }else if(system_RNG != "L'Ecuyer-CMRG" & use.parallel==TRUE){

      warning(paste("You want to parallize the subsampling but you did not use the RNG >>L'Ecuyer-CMRG<<: results might not be reproducable.", sep=""))

    }


    if(use.parallel==TRUE){

      if(numCores >= detectCores()){
        warning("You use at least as manes cores as there are (logical) cores available. Might slow down calculations.")
      }


      if(parallelization.type == "parLapply"){

        cl <- makeCluster(type = cluster.type, numCores)
        clusterSetRNGStream(cl = cl, iseed = seed)

        clusterExport(cl,)

      }else if(parallelization.type == "mclapply"){

        RNGkind(rng)
        set.seed(seed)

        subsample_results <- mclapply(X = 1:m, mc.cores = numCores, FUN =  function(i){

          ins_positions_subsample <-

              sort(sample(ins.positions, size=length(ins.positions)*d, replace = F))


          if(method == "Binomial"){

            results_tunings <- lapply(weightings, function(w){
              result_w <- Binomial(ins.positions = ins_positions_subsample,
                       gene.names,
                       gene.starts,
                       gene.stops,
                       num.ins.per.gene,
                       genome.length,
                       weighting=w)

              result_w$subsample_number <- i
              result_w
            })

          }else if(method == "ConNIS"){

            results_tunings <- lapply(weightings, function(w){
              result_w <- ConNIS(ins.positions = ins_positions_subsample,
                     gene.names,
                     gene.starts,
                     gene.stops,
                     num.ins.per.gene,
                     genome.length,
                     weighting=w)

              result_w$subsample_number <- i
              result_w
            })

          }else if(method == "Geometric"){

            results_tunings <- lapply(weightings, function(w){
              result_w <- Geometric(ins.positions = ins_positions_subsample,
                     gene.names,
                     gene.starts,
                     gene.stops,
                     num.ins.per.gene,
                     genome.length,
                     weighting=w)

              result_w$subsample_number <- i
              result_w
            })

          }else if(method == "Tn5Gaps"){

            results_tunings <- lapply(weightings, function(w){
              result_w <- Tn5Gaps(ins.positions = ins_positions_subsample,
                        gene.names,
                        gene.starts,
                        gene.stops,
                        genome.length,
                        weighting=w)

              result_w$subsample_number <- i
              result_w
            })

          }

          do.call(rbind, results_tunings)

        })

      }

      subsample_results

    }else{

      subsample_results <- lapply(X = 1:m, FUN =  function(i){

        ins_positions_subsample <-

          sort(sample(ins.positions, size=length(ins.positions)*d, replace = F))


        if(method == "Binomial"){

          results_tunings <- lapply(weightings, function(w){
            result_w <- Binomial(ins.positions = ins_positions_subsample,
                     gene.names,
                     gene.starts,
                     gene.stops,
                     num.ins.per.gene,
                     genome.length,
                     weighting=w)

            result_w$subsample_number <- i
            result_w
          })

        }else if(method == "ConNIS"){

          results_tunings <- lapply(weightings, function(w){
            result_w <- ConNIS(ins.positions = ins_positions_subsample,
                   gene.names,
                   gene.starts,
                   gene.stops,
                   num.ins.per.gene,
                   genome.length,
                   weighting=w)

            result_w$subsample_number <- i
            result_w
          })

        }else if(method == "Geometric"){

          results_tunings <- lapply(weightings, function(w){
            result_w <- Geometric(ins.positions = ins_positions_subsample,
                      gene.names,
                      gene.starts,
                      gene.stops,
                      num.ins.per.gene,
                      genome.length,
                      weighting=w)

            result_w$subsample_number <- i
            result_w
          })

        }else if(method == "Tn5Gaps"){

          results_tunings <- lapply(weightings, function(w){
            result_w <- Tn5Gaps(ins.positions = ins_positions_subsample,
                    gene.names,
                    gene.starts,
                    gene.stops,
                    genome.length,
                    weighting=w)

            result_w$subsample_number <- i
            result_w
          })

        }

        do.call(rbind, results_tunings)

      })

      subsample_results

    }
    subsample_results
  }



