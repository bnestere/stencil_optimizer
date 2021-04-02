library("tidyverse")
library("ggpubr")
library("ggpmisc")
library("optparse")
library("rlang")
library("cliplot")
library("hablar")
library("cowplot")
library("broom")
library("colorblindr")
library("yardstick")
library("firatheme")
library("lares")
library("moderndive")
library("skimr")
library("purrr")
library("parsnip")
library("workflows")
library("recipes")
library("pracma")

source("lib/opts.r")
source("lib/model.r")

opt <- cli_process(domain_opts=stencil_clioptions, 
                   default_infile='../stencilbench.semenuk.cfgd.csv', 
                   default_outfile='stenplot.png', 
                   default_x='NumThreads', default_y='MFLOPS', 
                   default_groupby='optId', default_facetvar='eqId')
gbsym <- rlang::sym(opt$groupby)
xsym <- rlang::sym(opt$xaxis)
ysym <- rlang::sym(opt$yaxis)
facetvarsym <- rlang::sym(opt$facetvar)

if(is.null(opt$infile)) {
  stop('no infile')
}

model_input <- data.frame(D = 2, OA = 4, TO = 0, tmax = 32, nThreads = 16, var_coefs=T, div=1)

print("beginning")

dat <- read.csv(opt$infile) %>%
  mutate(div = case_when(
                         divsLup > 0 ~ as.numeric(divsLup),
                         multsLup-addsLup==1 ~ as.numeric(0),
                         TRUE~as.numeric(multsLup-addsLup-1))) %>%
  mutate(var_coefs = case_when(
                               optId=='TDAC' ~ F,
                               (ldMats)>0 ~ T,
                               TRUE~F)) %>%
  mutate(optId = case_when((multsLup-addsLup)>=2 ~ paste0("MI + ",optId),
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
  #filter(optId != 'orig') %>%
  filter(!(optId == 'OMP' & nThreads == 1)) %>%
  group_by(optId, D, OA, TO, var_coefs, div, nThreads, tmax, InputSize) %>%
  summarize(mspeedup=mean(Speedup), stdev=sd(Speedup)) %>%
  ungroup() 

print("made it here")

bestopts <- dat %>%
  group_by(D, OA, TO, var_coefs, div, nThreads, tmax, InputSize) %>%
  #group_by(D, OA, TO, var_coefs, div, tmax, InputSize) %>%
  filter(mspeedup == max(mspeedup)) %>%
  summarize(bestopt = optId, speedup = mspeedup) %>%
  #summarize(bestopt = optId, speedup = mspeedup, bestNThreads = nThreads) %>%
  ungroup()

dat <- dat %>%
  mutate(nts = nThreads^2) 

idx_train = sample(1:nrow(dat), floor(0.95*nrow(dat)))
idx_test = setdiff(1:nrow(dat), idx_train)

dat_train = dat[idx_train,]
dat_test = dat[idx_test,]

print(bestopts)

  #mutate(fit_wf = map(data, ~lm(mspeedup ~ OA * TO * tmax * (nThreads + nts) * div * var_coefs * InputSize, data=.)))
predict_nthreads <- function(ntfitter, row, search_opt, search_dim) {
#  if(nThreads == 1) { 
#    pred_opt <- 'R'
#  #} else if(row[['var_coefs']] == FALSE & nThreads == 16) {
#  #  pred_opt <- 'TDAC'
#  } else {

  nt_wf <- ntfitter %>%
    filter(optId == search_opt & D == search_dim)

#  print("Predict nthreads")
#  print(nt_wf)

  max_val = 0
  max_nt = 0

  

  for(i in c(1,2,4,8,16)) {
    #pval <- predict(nt_wf[1,][['fit_wf']][[1]], list(nThreads=i, nts = i^2))
    prow <- list(nThreads=i,
                 nts=i^2,
                 OA=row[['OA']],
                 TO=row[['TO']],
                 tmax=row[['tmax']],
                 div=row[['div']],
                 var_coefs=row[['var_coefs']],
                 InputSize=row[['InputSize']])

    pval <- predict(nt_wf[1,][['fit_wf']][[1]], prow)

#    print(paste(search_opt, "with nt", i, "should have perf", pval))

    if(pval > max_val) {
      max_val <- pval
      max_nt <- i
    }
  }

  list(pnt=max_nt,pspeedup=max_val)

#  4
}

predict_best_opt <- function(fitter, row) {
  pred_opt <- ''
  nThreads <- row[['nThreads']]

#  if(nThreads == 1) { 
#    pred_opt <- 'R'
#  #} else if(row[['var_coefs']] == FALSE & nThreads == 16) {
#  #  pred_opt <- 'TDAC'
#  } else {

  if(1) {

    max_val = 0
    for(i in seq_along(1:nrow(fitter))) {
      orow <- fitter[i,]
      optId <- orow[['optId']]
      if(strcmp(optId,'TDAC') | strcmp(optId,'MI + TDAC')) {
        next 
      }

      pval <- predict(orow[['fit_wf']][[1]], row)


      if(pval > max_val) {
        max_val <- pval
        pred_opt <- orow[['optId']]
      }
    }

    #print(row)
    #print(pred_opt)

  }

  #if(nThreads > 1 & row[['D']] == 1) {
  #  if(strcmp(pred_opt, 'R') | strcmp(pred_opt,'MI + R')) {
  #    pred_opt <- paste0(pred_opt, ' + Omp')
  #  }
  #}

  if(row[['div']] > 0) {
    if(!startsWith(pred_opt, 'MI +')) {
      pred_opt <- paste0('MI + ',pred_opt)
    }
  } else {
    if(startsWith(pred_opt, 'MI')) {
      pred_opt <- sub('MI + ', '', pred_opt)
    }
  }

  pred_opt
}


fitter <- dat_train %>%
  group_by(optId) %>%
  nest() %>%
  mutate(fit_wf = map(data, ~lm(mspeedup ~ D * OA * TO * tmax * nThreads * div * var_coefs * InputSize, data=.)))

nt_fitter <- dat_train %>%
  group_by(optId, D) %>%
  nest() %>%
  mutate(fit_wf = map(data, ~lm(mspeedup ~ OA * TO * tmax * (nThreads + nts) * div * var_coefs * InputSize, data=.)))
  #mutate(fit_wf = map(data, ~lm(mspeedup ~ nThreads + nts, data=.)))

ncorrect <- 0
nwrong <- 0
nclose <- 0

#for(i in seq_along(1:nrow(dat_test))) {
for(i in seq_along(1:576)) {
  row = dat_test[i,]
  tD = row[['D']]
  tOA = row[['OA']]
  tTO = row[['TO']]
  ttmax = row[['tmax']]
  tnThreads = row[['nThreads']]
  tdiv = row[['div']]
  tvar_coefs = row[['var_coefs']]
  tInputSize = row[['InputSize']]

  best_opt <- subset(bestopts, D == tD & OA == tOA & TO == tTO & tmax == ttmax & nThreads == tnThreads & div == tdiv & var_coefs == tvar_coefs & InputSize == tInputSize, select=(bestopt))
  #best_opt <- subset(bestopts, D == tD & OA == tOA & TO == tTO & tmax == ttmax & div == tdiv & var_coefs == tvar_coefs & InputSize == tInputSize, select=(bestopt))

  pbo <- predict_best_opt(fitter, row)
  #pnt <- predict_nthreads(nt_fitter, row, pbo, row[['D']])
  prediction <- predict_nthreads(nt_fitter, row, pbo, row[['D']])
  pnt <- prediction[['pnt']]
  pred_speedup <- prediction[['pspeedup']]

  pred_opt_config <- list(optId=pbo,
                          nThreads=pnt,
                          D=tD,
                          OA=tOA,
                          TO=tTO,
                          tmax=ttmax,
                          nThreads=tnThreads,
                          div=tdiv,
                          var_coefs=tvar_coefs,
                          InputSize=tInputSize)

#  print(paste("Opt", pbo, "should run with", pnt, "threads"))

  # Override if we find a real value in our dataset
  real_speedup <- pred_speedup

  real_val_subset <- subset(dat, optId==pbo & D == tD & OA == tOA & TO == tTO & tmax == ttmax & nThreads == pnt & div == tdiv & var_coefs == tvar_coefs & InputSize == tInputSize, select=(mspeedup))
  if(nrow(real_val_subset) > 0) {
#    print(paste("We have a real value from subset of length", length(real_val_subset),'and is null?'))
#    print(real_val_subset)
    real_speedup <- as.numeric(real_val_subset[[1]])
  } else {
    print("No real data value found for row")
    print(paste("Length is", length(real_val_subset),'and is null?', is.null(real_val_subset)))
    print(row)
  }

  #print("real val subset")
  #print(real_val_subset)

#  if(is.null(real_speedup)) {
#    print(paste("Could not find value for row with prediction pbo:",pbo,"and pnt:",pnt))
#    print(row)
#    real_speedup <- pred_speedup
#  }
  #pred_speedup <- dat %>%
  #  filter(optId==pbo, D == tD & OA == tOA & TO == tTO & tmax == ttmax & nThreads == pnt & div == tdiv & var_coefs == tvar_coefs & InputSize == tInputSize)

  best_speedup <- as.numeric(subset(bestopts, D == tD & OA == tOA & TO == tTO & tmax == ttmax & nThreads == tnThreads & div == tdiv & var_coefs == tvar_coefs & InputSize == tInputSize, select=(speedup))[[1]])


#  print(paste("Testing row with real speedup:", real_speedup, "and best speedup", best_speedup))
#  print(row)
  diff <- abs((real_speedup - best_speedup)/best_speedup)
  if(real_speedup <= best_speedup & diff > 0.05) {

    print(paste("test failed"))
    print(row)
    print(paste("predicted:", pbo, ", best:",best_opt))
    print(paste("Real speedup:", real_speedup, "Predicted Speedup: ",pred_speedup, "Best Speedup:", best_speedup))
    #print(paste('predicted:',pred_speedup,'best:',best_speedup, "diff:",diff))

    nwrong <- nwrong + 1

    #print(pred_speedup)
    #print(best_speedup)

    if(diff < 0.1) {
      nclose <- nclose + 1
      print("It is close though")
    }
  } else {
    ncorrect <- ncorrect + 1
  }



#  if(!strcmp(pbo,best_opt)) {
#
#
#    #pred_speedup <- dat %>% filter(optId==pbo & D==tD & OA==tOA & TO == tTO & tmax == ttmax & nThreads == tnThreads & div == tdiv & var_coefs == tvar_coefs & InputSize == tInputSize) %>%
#    pred_speedup <- dat %>% filter(optId==pbo & D==tD & OA==tOA & TO == tTO & tmax == ttmax & div == tdiv & var_coefs == tvar_coefs & InputSize == tInputSize) %>%
#      group_by(optId, D, OA, TO, var_coefs, div, tmax, InputSize) %>%
#      #group_by(optId, D, OA, TO, var_coefs, div, nThreads, tmax, InputSize) %>%
#      #group_by(across(c(-mspeedup, -nThreads))) %>%
#      #slice_max(mspeedup) %>%
#      filter(mspeedup == max(mspeedup) & nThreads <= tnThreads) %>%
#      ungroup() %>%
#      #arrange(optId, D, OA, TO, tmax, var_coefs, div, nThreads, tmax, InputSize, mspeedup, stdev) 
#      select(mspeedup,stdev)
#    print("Tmp pred_speedup is")
#    print(pred_speedup)
#
#    #best_speedup <- subset(bestopts, D == tD & OA == tOA & TO == tTO & tmax == ttmax & div == tdiv & var_coefs == tvar_coefs & InputSize == tInputSize, select=(speedup))
#    best_speedup <- subset(bestopts, D == tD & OA == tOA & TO == tTO & tmax == ttmax & nThreads == tnThreads & div == tdiv & var_coefs == tvar_coefs & InputSize == tInputSize, select=(speedup))
#
#
#
#    #print(pred_speedup[[2]])
#    #print(pred_speedup[[1]])
#    #print(best_speedup[[1]])
#
#    if(nrow(pred_speedup) == 0) next
#
#    pred_speedup <- as.numeric(pred_speedup[[1]])
#    best_speedup <- as.numeric(best_speedup[[1]])
#    
#    diff <- abs((pred_speedup - best_speedup)/best_speedup)
#    if(diff > 0.05) {
#
#      print(paste("test failed"))
#      print(row)
#      print(paste("predicted:", pbo, ", best:",best_opt))
#      print(paste('predicted:',pred_speedup,'best:',best_speedup, "diff:",diff))
#
#      nwrong <- nwrong + 1
#
#      #print(pred_speedup)
#      #print(best_speedup)
#
#      if(diff < 0.1) {
#        nclose <- nclose + 1
#        print("It is close though")
#      }
#    } else {
#      ncorrect <- ncorrect + 1
#    }
#
#  } else {
#    ncorrect <- ncorrect + 1
#  }
#
#  nchecks <- nwrong + ncorrect
#  if(nchecks %% 10 == 0) {
#    print(paste0(nchecks,'/',nrow(dat_test)))
#    print(paste('prediction accuracy',ncorrect/nchecks,'; closeness accuracy', nclose/nchecks))
#    print(paste('ncorrect',ncorrect,'; nchecks',nchecks))
#  }
  #print('row best op real/predicated')
  #print(row)
  #print(paste(best_opt, pbo))
}

nchecks <- ncorrect+nwrong
print(paste('prediction accuracy',ncorrect/nchecks,'; closeness accuracy', nclose/nchecks))
print(paste('ncorrect',ncorrect,'; nchecks',nchecks))
