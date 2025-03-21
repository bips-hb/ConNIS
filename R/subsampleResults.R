#' Results based on subsamples
#'
#' 'subsampleResults' draws \eqn{m} subsamples of the observerd IS and
#' calculates for each subsample the results for a user selected analysis
#' method.
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
#' @param weights A sequence of weight values that are applied to the
#' selected method.
#' @param m Number of sub samples.
#' @param d Proportion of the original IS used as sub samples.
#' @param use.parallelization Should the calculations be done in parallel?
#' @param parallelization.type Which method should be used for the
#' parallelization? Available are `"mclapply"` (fork mechanism) and
#' `"parLapply"` (creation of sockets) of the `"parallel"` package.
#' @param numCores Number of cores to juse for parallel calculation. Should not
#' exit the number of available (logical) cores of the system.
#' @param cluster.type If `"parLapply"` is used as parallelization.type, a
#' cluster type like `"PSOCK"` needs tp be specified. Check which types are
#' supported on your system. See als `?makeCluster` of the `parallel` package.
#' @param seed Seed for subsampling. If `NULL` no seed is set and results might
#' not be reproducable.
#' @param set.rng Which random number generator (RNG) should be used. If NULL, the
#' default RNG of the system is used. NOTE: if use.parallelizationization is `TRUE` the
#' `"L'Ecuyer-CMRG"` value for RNG should be used. For details see the manual of
#' the `parallel` package.
#' @param keep.new.RNG If the RNG is changed, should it kept for the rest of the
#' session? Default is `FALSE` (recommended).
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
#' # run subsumple procedure with ConNIS (sequential)
#' subsampleResults(method="ConNIS",
#'                   ins.positions = random_is,
#'                   gene.names = genes,
#'                   gene.starts = starts,
#'                   gene.stops = stops,
#'                   genome.length = genome_length,
#'                   weights = c(0.3, 0.6),
#'                   m = 2,
#'                   d = 0.5,
#'                   use.parallelization = FALSE,
#'                   seed = 1)
#'
#' # Example with parallelization (mclapply()); use detectCores()-1 to set
#' # workers
#' \dontrun{
#' subsampleResults(method="ConNIS",
#'                   ins.positions = random_is,
#'                   gene.names = genes,
#'                   gene.starts = starts,
#'                   gene.stops = stops,
#'                   genome.length = genome_length,
#'                   weights = c(0.3, 0.6),
#'                   m = 2,
#'                   d = 0.5,
#'                   use.parallelization = T,
#'                   parallelization.type = "mclapply",
#'                   set.rng = "L'Ecuyer-CMRG",
#'                   numCores=max(1,detectCores()-1),
#'                   seed = 1)
#' }
#'
#' @export

subsampleResults <- function(method="ConNIS",
                             ins.positions,
                             gene.names,
                             gene.starts,
                             gene.stops,
                             num.ins.per.gene = NULL,
                             genome.length,
                             weights = 1,
                             m = 100,
                             d = 0.5,
                             use.parallelization=FALSE,
                             parallelization.type="mclapply",
                             numCores=3,
                             cluster.type = NULL,
                             seed = NULL,
                             set.rng = NULL,
                             keep.new.RNG=FALSE){

  current_RNG <-RNGkind()[1]
  proceed <- 1

  # catch errors:
  if(d<=0 | d >= 1 | length(d)!=1){
    stop("d has to be a single value >0 and <1.")
  }

  if(use.parallelization == TRUE & parallelization.type == "parLapply" & is.null(cluster.type)){
    stop("Need to set a valide cluster type for using parLapply.")
  }

  # print warnings if RNG or seed is not set appropriately
  if(is.null(seed)){
    proceed <- menu(c("Yes", "No"),
                    title = "No seed set, results might not be reproducable. Do you want to proceed?")
  }

  if(is.null(set.rng) & use.parallelization == FALSE){
    warning("No RNG but a seed was set. Used current RNG this session: ", current_RNG)
  }

  if(is.null(set.rng)){
    if(use.parallelization == TRUE & parallelization.type == "mclapply"){
      proceed <- menu(c("Yes", "No"),
                      title = "Selected 'mclapply' for subsampling but RNG has not been set to >>L'Ecuyer-CMRG<<: results might not be reproducable. Do you want to proceed?")
    }
  }

  if(!is.null(set.rng)){
    if(set.rng != "L'Ecuyer-CMRG" & use.parallelization == TRUE & parallelization.type == "mclapply"){
      proceed <- menu(c("Yes", "No"),
                      title = "Selected 'mclapply' for subsampling but RNG has not been set to >>L'Ecuyer-CMRG<<: results might not be reproducable. Do you want to proceed?")
    }
  }

  if(proceed==1){

    if(use.parallelization==TRUE){

      if(numCores >= detectCores()){
        warning("You use at least as manes cores as there are (logical) cores available. Might slow down calculations.")
      }


      if(parallelization.type == "parLapply"){

        # make sure that the cluster is closed in any case
        on.exit(stopCluster(cl))

        # make cluster
        cl <- makeCluster(type = cluster.type, numCores)

        # set seed for RNG
        clusterSetRNGStream(cl = cl, iseed = seed)

        clusterExport(cl = cl,
                      envir = environment(),
                      varlist = list(
                        "method",
                        "ins.positions",
                        "gene.names",
                        "gene.starts",
                        "gene.stops",
                        "num.ins.per.gene",
                        "genome.length",
                        "weights",
                        "m",
                        "d",
                        "prob_seq_misses",
                        "freq_seq_misses",
                        "pgumbel",
                        "Binomial",
                        "ConNIS",
                        "Geometric",
                        "Tn5Gaps"
                      ))

        clusterEvalQ(cl, {
          library(tibble)
          library(gmp)}
        )

        subsample_results <- parLapply(cl = cl, X = 1:m, fun =  function(i){

          # draw subsamples of IS
          ins_positions_subsample <-
            sort(sample(ins.positions, size=length(ins.positions)*d, replace = F))


          if(method == "Binomial"){

            results_tunings <- lapply(weights, function(w){
              result_w <- Binomial(ins.positions = ins_positions_subsample,
                                   gene.names,
                                   gene.starts,
                                   gene.stops,
                                   num.ins.per.gene,
                                   genome.length,
                                   weight=w)

              result_w$subsample_number <- i
              result_w
            })

          }else if(method == "ConNIS"){

            results_tunings <- lapply(weights, function(w){
              result_w <- ConNIS(ins.positions = ins_positions_subsample,
                                 gene.names,
                                 gene.starts,
                                 gene.stops,
                                 num.ins.per.gene,
                                 genome.length,
                                 weight=w)

              result_w$subsample_number <- i
              result_w
            })

          }else if(method == "Geometric"){

            results_tunings <- lapply(weights, function(w){
              result_w <- Geometric(ins.positions = ins_positions_subsample,
                                    gene.names,
                                    gene.starts,
                                    gene.stops,
                                    num.ins.per.gene,
                                    genome.length,
                                    weight=w)

              result_w$subsample_number <- i
              result_w
            })

          }else if(method == "Tn5Gaps"){

            results_tunings <- lapply(weights, function(w){
              result_w <- Tn5Gaps(ins.positions = ins_positions_subsample,
                                  gene.names,
                                  gene.starts,
                                  gene.stops,
                                  genome.length,
                                  weight=w)

              result_w$subsample_number <- i
              result_w
            })

          }

          do.call(rbind, results_tunings)

        })


      }else if(parallelization.type == "mclapply"){

        # set RNG fpr mclapply; NOTE: other than "L'Ecuyer-CMRG" might not give
        # reproducable results
        RNGkind(set.rng)
        set.seed(seed)

        subsample_results <- mclapply(X = 1:m, mc.cores = numCores, FUN =  function(i){

          # draw subsamples of IS
          ins_positions_subsample <-
            sort(sample(ins.positions, size=length(ins.positions)*d, replace = F))


          if(method == "Binomial"){

            results_tunings <- lapply(weights, function(w){
              result_w <- Binomial(ins.positions = ins_positions_subsample,
                                   gene.names,
                                   gene.starts,
                                   gene.stops,
                                   num.ins.per.gene,
                                   genome.length,
                                   weight=w)

              result_w$subsample_number <- i
              result_w
            })

          }else if(method == "ConNIS"){

            results_tunings <- lapply(weights, function(w){
              result_w <- ConNIS(ins.positions = ins_positions_subsample,
                                 gene.names,
                                 gene.starts,
                                 gene.stops,
                                 num.ins.per.gene,
                                 genome.length,
                                 weight=w)

              result_w$subsample_number <- i
              result_w
            })

          }else if(method == "Geometric"){

            results_tunings <- lapply(weights, function(w){
              result_w <- Geometric(ins.positions = ins_positions_subsample,
                                    gene.names,
                                    gene.starts,
                                    gene.stops,
                                    num.ins.per.gene,
                                    genome.length,
                                    weight=w)

              result_w$subsample_number <- i
              result_w
            })

          }else if(method == "Tn5Gaps"){

            results_tunings <- lapply(weights, function(w){
              result_w <- Tn5Gaps(ins.positions = ins_positions_subsample,
                                  gene.names,
                                  gene.starts,
                                  gene.stops,
                                  genome.length,
                                  weight=w)

              result_w$subsample_number <- i
              result_w
            })

          }

          do.call(rbind, results_tunings)

        })

      }

      subsample_results

    }else{

      # run subsample pricedure sequentally -  might be slow !!!
      set.seed(seed)

      subsample_results <- lapply(X = 1:m, FUN =  function(i){

        # draw subsamples of IS
        ins_positions_subsample <-
          sort(sample(ins.positions, size=length(ins.positions)*d, replace = F))

        if(method == "Binomial"){

          results_tunings <- lapply(weights, function(w){
            result_w <- Binomial(ins.positions = ins_positions_subsample,
                                 gene.names,
                                 gene.starts,
                                 gene.stops,
                                 num.ins.per.gene,
                                 genome.length,
                                 weight=w)

            result_w$subsample_number <- i
            result_w
          })

        }else if(method == "ConNIS"){

          results_tunings <- lapply(weights, function(w){
            result_w <- ConNIS(ins.positions = ins_positions_subsample,
                               gene.names,
                               gene.starts,
                               gene.stops,
                               num.ins.per.gene,
                               genome.length,
                               weight=w)

            result_w$subsample_number <- i
            result_w
          })

        }else if(method == "Geometric"){

          results_tunings <- lapply(weights, function(w){
            result_w <- Geometric(ins.positions = ins_positions_subsample,
                                  gene.names,
                                  gene.starts,
                                  gene.stops,
                                  num.ins.per.gene,
                                  genome.length,
                                  weight=w)

            result_w$subsample_number <- i
            result_w
          })

        }else if(method == "Tn5Gaps"){

          results_tunings <- lapply(weights, function(w){
            result_w <- Tn5Gaps(ins.positions = ins_positions_subsample,
                                gene.names,
                                gene.starts,
                                gene.stops,
                                genome.length,
                                weight=w)

            result_w$subsample_number <- i
            result_w
          })

        }

        do.call(rbind, results_tunings)

      })

      subsample_results

    }

    # set the old RNG
    if(keep.new.RNG == FALSE){
      RNGkind(current_RNG)
    }

    # return subsample results
    subsample_results

  }else{
    print("Did not run subsampleResults()")
  }

}

