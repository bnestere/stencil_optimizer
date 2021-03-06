library("tidyverse")
library("optparse")
library("rjson")
library("rlist")

source("lib/opts.r")
source("lib/model.r")

options(warn = -1)

opt_parser = OptionParser(option_list=stencil_clioptions)
opt = parse_args(opt_parser)

print(opt$`in`)

if(is.null(opt$`in`)) {
  stop('no input stencil file')
}

# ----------------------------------------------------------------------------
# -----  Parse the input stencil code to extract domain-level 
# ----- props and opportunity indicators 
# ----------------------------------------------------------------------------
pcg_cmd <- 'pcg'
pcg_args <- c(paste0("-pin=", opt$`in`), "FindDomainProps.pt")
system2(pcg_cmd, pcg_args, wait=TRUE)

stencil_json <- "stencil_props.json"
exe_conf_json <- "exe_config.json"
sten_props <- fromJSON(file=stencil_json)
exe_conf <- fromJSON(file=exe_conf_json)

# ----------------------------------------------------------------------------
# -----  Shape the data that builds the regression model ---------------------
# ----------------------------------------------------------------------------
stencil_csv <- "stencilbench.semenuk.cfgd.csv"
dat <- read.csv(stencil_csv) %>%
  mutate(div = case_when(
                         divsLup > 0 ~ as.numeric(divsLup),
                         multsLup-addsLup==1 ~ as.numeric(0),
                         TRUE~as.numeric(multsLup-addsLup-1))) %>%
  mutate(var_coefs = case_when(
                               optId=='TDAC' ~ F,
                               (ldMats)>0 ~ T,
                               TRUE~F)) %>%
#  filter(optId != 'R + Omp') %>%
#  filter(optId != "TDAC") %>%
  mutate(optId = case_when(optId == 'orig' & (multsLup-addsLup)>=2 ~ paste0("Multiplicative Inversion"),
                           startsWith(paste0(optId), 'R') ~ paste0('Rectangular Tiling'),
                           startsWith(paste0(optId), 'Pd') ~ paste0('Partial Diamond Tiling'),
                           startsWith(paste0(optId), 'PD') ~ paste0('Partial Diamond Tiling'),
                           startsWith(paste0(optId), 'Fd') ~ paste0('Full Diamond Tiling'),
                           startsWith(paste0(optId), 'FD') ~ paste0('Full Diamond Tiling'),
                           startsWith(paste0(optId), 'TDA') ~ paste0('Cache-oblvious Divide and Conquer'),
                           T ~ as.character(optId))
         ) %>%
  mutate(FlopsPerLup = addsLup+subsLup+multsLup+divsLup) %>%
  mutate(InputSize = case_when(D == 1 ~ xmax,
                               D == 2 ~ xmax*ymax,
                               D == 3 ~ xmax*ymax*zmax)) %>%
  mutate(TotalLups = tmax*InputSize) %>%
  mutate(GFLOPS = round((((FlopsPerLup * TotalLups)/duration)/1000000000), digits=2)) %>%
  group_by(D, OA, TO, div, var_coefs, tmax, InputSize, optId) %>%
  mutate(mean_GFLOPS=mean(GFLOPS)) %>%
  ungroup() %>%

  group_by(D, OA, TO, div, var_coefs, tmax, InputSize) %>%
  #mutate(Original=max(case_when(optId=='orig' ~ mean_GFLOPS, TRUE ~ 0)),
  mutate(Original=max(case_when(optId=='OMP' & nThreads==1 ~ mean_GFLOPS, TRUE ~ 0)),
         Speedup=GFLOPS/Original) %>%
  ungroup() %>%
  filter(optId != 'orig') %>%
  group_by(optId, D, OA, TO, var_coefs, div, nThreads, tmax, InputSize) %>%
  summarize(mspeedup=mean(Speedup), mgflops=mean(GFLOPS), stdev=sd(Speedup)) %>%
  ungroup() 

print(dat)

dat <- dat %>%
  mutate(nts = nThreads^2)


print("Pre Build the Regression Model")
# ----------------------------------------------------------------------------
# -----  Build the Regression Model 
# ----------------------------------------------------------------------------
fitter <- dat %>%
  group_by(optId) %>%
  nest() %>%
  mutate(fit_wf = map(data, ~lm(mspeedup ~ D * OA * TO * var_coefs * tmax * (nThreads + nts) * InputSize, data=.)))

print("Built the Regression Model")

ncorrect <- 0
nwrong <- 0
nclose <- 0

modelInputSize <- 0
floatSize <- 8
if(sten_props$dim == 1) {
  modelInputSize <- (exe_conf$xmax)
}
if(sten_props$dim == 2) {
  modelInputSize <- (exe_conf$xmax * exe_conf$ymax)
}
if(sten_props$dim == 3) {
  modelInputSize <- (exe_conf$xmax * exe_conf$ymax * exe_conf$zmax)
}

max_n_threads <- exe_conf$nThreads
max_nts <- max_n_threads^2

model_input <-data.frame(D = sten_props$dim,
                         OA = sten_props$oa,
                         TO = sten_props$to,
                         var_coefs = as.logical(sten_props$dd),
                         nThreads = max_n_threads,
                         nts = exe_conf$nThreads^2,
                         tmax = exe_conf$nTimesteps,
                         InputSize = modelInputSize)

print("model input is")
print(model_input)
                         
filter_invalid_opts <- function(opts_list, indicators) {

  opts <- opts_list
  if(!indicators$isLinear) {
    if(length(opts) == 0) return(NULL)
    if("Rectangular Tiling" %in% opts) {
      if(length(opts) == 3) {
        opts <- NULL
      } else {
        opts <- list.remove(opts, c('fd','r','pd'))
      }
    }
  }

  if(!indicators$hasConstDiv) {
    if(length(opts) == 0) return(NULL)
    if("Multiplicative Inversion" %in% opts) {
      if(length(opts) == 1) {
        opts <- NULL
      } else {
        opts <- list.remove(opts, c('mi'))
      }
    }
  }
  
  if(!indicators$isSerial) {
    if(length(opts) == 0) return(NULL)
    if("OMP" %in% opts) {
      if(length(opts) == 1) {
        opts <- NULL
      } else {
        opts <- list.remove(opts, c('omp'))
      }
    }
  }

  opts
}


predict_strategy <- function(model_input, inp_sten_props, max_n_threads, fitter) {
  valid_opts <- list(fd="Full Diamond Tiling", omp="OMP", pd="Partial Diamond Tiling", r="Rectangular Tiling", mi="Multiplicative Inversion")
  valid_opts <- filter_invalid_opts(valid_opts, inp_sten_props)
  
  opt_strat = c()
  best_n_threads <- max_n_threads
  
  n_thread_counts <- c(1,2,4,8,16)
  highest_perf <- 0
  # ----------------------------------------------------------------------------
  # ----- Use the prediction model to specialize an optimization strategy ------
  # ----------------------------------------------------------------------------
  while(length(valid_opts) > 0) {
  
    max_pred <- 0
    best_opt <- ''
  
    for(i in seq_along(1:nrow(fitter))) {
      orow <- fitter[i,]
      optId <- orow[['optId']]
    
      if(!(optId %in% valid_opts)) next
    
      print("Showing orow")
      print(orow)
  
      for(nt in n_thread_counts) {
        if(nt > max_n_threads) break
  
        model_input$nThreads <- nt
        model_input$nts <- nt^2
        print(paste("Going to test nthreads", nt))
  
        pval <- round(predict(orow[['fit_wf']][[1]], model_input), 2)
        pred_opt <- orow[['optId']]
    
        if(pval > max_pred & !startsWith(optId, 'Cach')) {
          max_pred <- pval
          best_opt <- pred_opt 
  
          if(pval > highest_perf) {
            best_n_threads <- nt
            highest_perf <- pval
          }
        }
    
        if(startsWith(optId, 'Cach'))  pval <- pval / 3
  
        print(paste0('Predicted speedup (nt:',nt,'): ', pred_opt, ': ',pval, 'x'))
      }
    }
  
    print('')
    #print(paste0('Recommended optimization: ', best_opt))
  
    opt_strat <- c(opt_strat, best_opt)
  
    if(best_opt == "Full Diamond Tiling" || 
       best_opt == "Partial Diamond Tiling" || 
       best_opt == "Rectangular Tiling") {
  
      inp_sten_props$isLinear <- FALSE
    }
  
    if(best_opt == "OMP") {
      inp_sten_props$isSerial <- FALSE
    }
  
    if(best_opt == "Multiplicative Inversion") {
      inp_sten_props$hasConstDiv <- FALSE
    }
  
    valid_opts <- filter_invalid_opts(valid_opts, inp_sten_props)
  }

  list(opt_strat, best_n_threads)
}

#valid_opts <- filter_invalid_opts(valid_opts, sten_props)
#
#opt_strat = c()
#
#n_thread_counts <- c(1,2,4,8,16)
#best_n_threads <- max_n_threads
#highest_perf <- 0
## ----------------------------------------------------------------------------
## ----- Use the prediction model to specialize an optimization strategy ------
## ----------------------------------------------------------------------------
#while(length(valid_opts) > 0) {
#
#  max_pred <- 0
#  best_opt <- ''
#
#  for(i in seq_along(1:nrow(fitter))) {
#    orow <- fitter[i,]
#    optId <- orow[['optId']]
#  
#    if(!(optId %in% valid_opts)) next
#  
#    print("Showing orow")
#    print(orow)
#
#    for(nt in n_thread_counts) {
#      if(nt > max_n_threads) break
#
#      model_input$nThreads <- nt
#      model_input$nts <- nt^2
#      print(paste("Going to test nthreads", nt))
#
#      pval <- round(predict(orow[['fit_wf']][[1]], model_input), 2)
#      pred_opt <- orow[['optId']]
#  
#      if(pval > max_pred & !startsWith(optId, 'Cach')) {
#        max_pred <- pval
#        best_opt <- pred_opt 
#
#        if(pval > highest_perf) {
#          best_n_threads <- nt
#          highest_perf <- pval
#        }
#      }
#  
#      if(startsWith(optId, 'Cach'))  pval <- pval / 3
#
#      print(paste0('Predicted speedup (nt:',nt,'): ', pred_opt, ': ',pval, 'x'))
#    }
#  }
#
#  print('')
#  print(paste0('Recommended optimization: ', best_opt))
#
#  opt_strat <- c(opt_strat, best_opt)
#
#  if(best_opt == "Full Diamond Tiling" || 
#     best_opt == "Partial Diamond Tiling" || 
#     best_opt == "Rectangular Tiling") {
#
#    sten_props$isLinear <- FALSE
#  }
#
#  if(best_opt == "OMP") {
#    sten_props$isSerial <- FALSE
#  }
#
#  if(best_opt == "Multiplicative Inversion") {
#    sten_props$hasConstDiv <- FALSE
#  }
#
#  valid_opts <- filter_invalid_opts(valid_opts, sten_props)
#}

pred_info <- predict_strategy(model_input, sten_props, exe_conf$nThreads, fitter)
opt_strat <- pred_info[1]
best_n_threads <- pred_info[2]


# ----------------------------------------------------------------------------
# -----  Optimize the code ---------------------------------------------------
# ----------------------------------------------------------------------------
tile_id <- ""
tile_flag <- "--tile"
if("Rectangular Tiling" %in% opt_strat) {
  tile_flag <- paste(tile_flag, "--nodiamond-tile")
  tile_id <- 'rect'
} else if("Partial Diamond Tiling" %in% opt_strat) {
  tile_flag <- paste(tile_flag, "--diamond-tile")
  tile_id <- 'pd'
} else if("Full Diamond Tiling" %in% opt_strat) {
  tile_flag <- paste(tile_flag, "--full-diamond-tile")
  tile_id <- 'fd'
}
dim <- sten_props$dim
dim_id <- paste0(c(sten_props$dim, 'd'), collapse='')
tile_file <- paste0(c('tile',dim_id,tile_id,'sizes'),collapse='.')
tile_folder <- 'tile_sizes'
tile_path <- paste0(c(tile_folder, tile_file), collapse='/')
print(paste("using blocking file", tile_path))

cpargs <- c(tile_path, "tile.sizes")
system2('cp', args=cpargs)


parallel_flag <- ""
#if("Rectangular Tiling" %in% opt_strat) {
if("OMP" %in% opt_strat) {
  parallel_flag <- "--parallel"
}

# Multiplicative Inversion is performed previously through POET 



stage_file <- "poet_stage.c"
margs <- c(stage_file, tile_flag, parallel_flag, "--pet", "-o out.c")

command <- paste("/home/brandon/opt_devel/pluto/target/bin/polycc", paste0(margs, collapse=" "))
print(command)

system2('polycc', args=margs)

print("\n\n-------------------------------------------------")
print("-------------------------------------------------")
print("Successfully applied the following optimizations:")
print(" ")
print(opt_strat)
print(paste("Recommended Num Threads:", best_n_threads))
print("-------------------------------------------------")
print("-------------------------------------------------")

print("Model Input")
print(model_input)


#model_input <-data.frame(D = sten_props$dim,
#                         OA = sten_props$oa,
#                         TO = sten_props$to,
#                         var_coefs = as.logical(sten_props$dd),
#                         nThreads = max_n_threads,
#                         nts = exe_conf$nThreads^2,
#                         tmax = exe_conf$nTimesteps,
#                         InputSize = modelInputSize)

dims <- c(1, 2, 3)
oas <- c(2,4)
tos <- c(0,1,2)
var_coefs <- c(TRUE, FALSE)
tmaxs <- c(16, 32, 64, 128, 256)
input_sizes <- c(500000, 1000000, 2000000, 4000000)
max_threads <- 16

msten_props <- data.frame( isLinear = TRUE,
                          isSerial = TRUE,
                          hasConstDiv = FALSE)

#for(dim in dims) {
#  for(oa in oas) {
#    for(to in tos) {
#      for(vc in var_coefs) {
#        for(tmax in tmaxs) {
#          for(is in input_sizes) { 
#
#            input <- data.frame(D = dim,
#                         OA = oa,
#                         TO = to,
#                         var_coefs = vc,
#                         nThreads = max_threads,
#                         nts = max_threads^2,
#                         tmax = tmax,
#                         InputSize = is)
#            pred_info <- predict_strategy(input, msten_props, max_threads, fitter)
#            opt_strat <- pred_info[1]
#            best_n_threads <- pred_info[2]
#
#            if(best_n_threads != max_threads) {
#              print("")
#              print("Model Input")
#              print(input)
#              print("Predictions")
#              print(opt_strat)
#              print(paste("Recommended Num Threads:", best_n_threads))
#            }
#
#
#          }
#        }
#      }
#    }
#  }
#}
