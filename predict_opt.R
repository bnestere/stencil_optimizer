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

pcg_cmd <- 'pcg'
pcg_args <- c(paste0("-pin=", opt$`in`), "FindDomainProps.pt")
system2(pcg_cmd, pcg_args, wait=TRUE)

stencil_json <- "stencil_props.json"
exe_conf_json <- "exe_config.json"
sten_props <- fromJSON(file=stencil_json)
exe_conf <- fromJSON(file=exe_conf_json)

stencil_csv <- "stencilbench.csv"
dat <- read.csv(stencil_csv) %>%
  mutate(div = case_when(
                         divsLup > 0 ~ as.numeric(divsLup),
                         multsLup-addsLup==1 ~ as.numeric(0),
                         TRUE~as.numeric(multsLup-addsLup-1))) %>%
  mutate(var_coefs = case_when(
                               optId=='TDAC' ~ F,
                               (ldMats)>0 ~ T,
                               TRUE~F)) %>%
  filter(optId != 'R + Omp') %>%
  filter(optId != "TDAC") %>%
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
  mutate(Original=max(case_when(optId=='orig' ~ mean_GFLOPS, TRUE ~ 0)),
         Speedup=GFLOPS/Original) %>%
  ungroup() %>%
  filter(optId != 'orig') %>%
  group_by(optId, D, OA, TO, var_coefs, div, nThreads, tmax, InputSize) %>%
  summarize(mspeedup=mean(Speedup), mgflops=mean(GFLOPS), stdev=sd(Speedup)) %>%
  ungroup() 

fitter <- dat %>%
  group_by(optId) %>%
  nest() %>%
  mutate(fit_wf = map(data, ~lm(mspeedup ~ D * OA * TO * var_coefs * tmax * nThreads * InputSize, data=.)))

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

model_input <-data.frame(D = sten_props$dim,
                         OA = sten_props$oa,
                         TO = sten_props$to,
                         var_coefs = as.logical(sten_props$dd),
                         nThreads = exe_conf$nThreads,
                         tmax = exe_conf$nTimesteps,
                         InputSize = modelInputSize)

print("model input is")
print(model_input)
                         
filter_invalid_opts <- function(opts_list, indicators) {

  opts <- opts_list
  if(!sten_props$isLinear) {
    if(length(opts) == 0) return(NULL)
    if("Rectangular Tiling" %in% opts) {
      if(length(opts) == 3) {
        opts <- NULL
      } else {
        opts <- list.remove(opts, c('fd','r','pd'))
      }
    }
  }

  if(!sten_props$hasConstDiv) {
    if(length(opts) == 0) return(NULL)
    if("Multiplicative Inversion" %in% opts) {
      if(length(opts) == 1) {
        opts <- NULL
      } else {
        opts <- list.remove(opts, c('mi'))
      }
    }
  }
  
  if(!sten_props$isSerial) {
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

valid_opts <- list(fd="Full Diamond Tiling", omp="OMP", pd="Partial Diamond Tiling", r="Rectangular Tiling", mi="Multiplicative Inversion")
valid_opts <- filter_invalid_opts(valid_opts, sten_props)

opt_strat = c()

while(length(valid_opts) > 0) {

  max_pred <- 0
  best_opt <- ''

  for(i in seq_along(1:nrow(fitter))) {
    orow <- fitter[i,]
    optId <- orow[['optId']]
  
    if(!(optId %in% valid_opts)) next
  
    pval <- round(predict(orow[['fit_wf']][[1]], model_input), 2)
    pred_opt <- orow[['optId']]
  
    if(pval > max_pred & !startsWith(optId, 'Cach')) {
      max_pred <- pval
      best_opt <- pred_opt 
    }
  
    if(startsWith(optId, 'Cach'))  pval <- pval / 3

    print(paste0('Predicted speedup: ', pred_opt, ': ',pval, 'x'))
  }

  print('')
  print(paste0('Recommended optimization: ', best_opt))

  opt_strat <- c(opt_strat, best_opt)

  if(best_opt == "Full Diamond Tiling" || 
     best_opt == "Partial Diamond Tiling" || 
     best_opt == "Rectangular Tiling") {

    sten_props$isLinear <- FALSE
  }

  if(best_opt == "OMP") {
    sten_props$isSerial <- FALSE
  }

  if(best_opt == "Multiplicative Inversion") {
    sten_props$hasConstDiv <- FALSE
  }

  valid_opts <- filter_invalid_opts(valid_opts, sten_props)
}

tile_flag <- "--tile"
if("Rectangular Tiling" %in% opt_strat) {
  tile_flag <- paste(tile_flag, "--nodiamond-tile")
} else if("Partial Diamond Tiling" %in% opt_strat) {
  tile_flag <- paste(tile_flag, "--diamond-tile")
} else if("Full Diamond Tiling" %in% opt_strat) {
  tile_flag <- paste(tile_flag, "--full-diamond-tile")
}

parallel_flag <- ""
if("Rectangular Tiling" %in% opt_strat) {
  parallel_flag <- "--parallel"
}

# Multiplicative Inversion is performed previously through POET 


stage_file <- "poet_stage.c"
stage_file <- "poet_stage.c"
margs <- c(stage_file, tile_flag, parallel_flag, "--pet", "-o", "out.c")
command <- paste("/home/brandon/opt_devel/pluto/target/bin/polycc", margs)


system2('polycc', args=margs)

print("\n\n-------------------------------------------------")
print("-------------------------------------------------")
print("Successfully applied the following optimizations:")
print(" ")
print(opt_strat)
print("-------------------------------------------------")
print("-------------------------------------------------")
