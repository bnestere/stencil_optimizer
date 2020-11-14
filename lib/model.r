resources <- tribble(
  ~rid, ~rcapacity, ~rbandwidth, ~rnum,
  "L1",  32000,    201, 0,
  "L2",  1024000,  131, 0.25,
  #"L2",  1024000,  95.6,
  "LLC", 22528000, 39.4, 0.5,
  #"LLC", 22528000, 24.7,
  #"Memory", 100000000000, 3.6
  #"Memory", 100000000000, 7
  "Memory", 100000000000, 20.4, 0.75
  #"Memory", 100000000000, 18
  #"Memory", 100000000000, 13
)

bytes <- function(value, ds_entry) {
  value * ds_entry[['fpSize']]
}

get_ts_size_1d <- function(., xname='xmax', yname='ymax', zname = 'zmax') {
  # Original
  .[['TotalMats']] * .[[xname]]

  #Considering cache replacement policy
  #.[[xname]] - .[['LongDistanceMats']]
  #.[[xname]]
}
get_ts_size_2d <- function(., xname='xmax', yname='ymax', zname = 'zmax') {
  .[['TotalMats']] * .[[xname]] * .[[yname]]
}
get_ts_size_3d <- function(., xname='xmax', yname='ymax', zname = 'zmax') {
  .[['TotalMats']] * .[[xname]] * .[[yname]] * .[[zname]] 
}

get_reuse_1d <- function(., xname='xmax', yname='ymax', zname = 'zmax') {
  .[[xname]] + .[['TotalMats']]
}
get_reuse_2d <- function(., xname='xmax', yname='ymax', zname = 'zmax') {
  .[['Oa']] * .[[xname]] + (.[['TotalMats']]*.[[xname]])

  # Make this more general for when we add coefficient matrices
  #.[['Oa']] * .[[xname]] + (.[['Oa']] * (.[['TotalMats']]*.[[xname]]))
}
get_reuse_3d <- function(., xname='xmax', yname='ymax', zname = 'zmax') {
  (.[['Dim']] - 1) * .[['Oa']] * .[[xname]] * .[[yname]] + (.[['Oa']] * (.[['TotalMats']]*.[[xname]]*.[[yname]]))
}

get_ts_size <- function(entry, xname, yname, zname) {
  #print("entry is")
  #print(entry)
  if(entry[['Dim']] == 1) {
    get_ts_size_1d(entry, xname=xname)
  } 
  else if(entry[['Dim']] == 2) {
    get_ts_size_2d(entry, xname=xname, yname=yname)
  } 
  else if(entry[['Dim']] == 3) {
    get_ts_size_3d(entry, xname=xname, yname=yname, zname=zname)
  } 
  else {
    stop("unrecognized dimension")
  }
}

get_reuse_size <- function(entry, xname, yname, zname) {
  if(entry[['Dim']] == 1) {
    get_reuse_1d(entry, xname=xname, yname=yname, zname=zname)
  } 
  else if(entry[['Dim']] == 2) {
    get_reuse_2d(entry, xname=xname, yname=yname, zname=zname)
  } 
  else if(entry[['Dim']] == 3) {
    get_reuse_3d(entry, xname=xname, yname=yname, zname=zname)
  } 
  else {
    stop("unrecognized dimension")
  }
}

check_bottleneck_generic <- function(resource, entry, size_func, reuse_func) {
  size <- size_func(entry) * entry[['fpSize']]
  if(size < resource[['rcapacity']]) {
    return(TRUE)
  }
  FALSE 
}

print_if_2d <- function(...,entry) {
  print(paste(..., sep=' '))
}

block_rsrc_pred <- function(resource, entry) {
  # Does everything fit within the resource?
  ts_size <- bytes(get_ts_size(entry, xname='bx', yname='by', zname='bz'), entry)
  if(ts_size < resource[['rcapacity']]) {
    #print_if_2d('entry has ts within capacity', entry[['Dim']], ts_size, resource[['rid']], resource[['rcapacity']], entry=entry)
    return(TRUE)
  }

#  reuse_size <- bytes(get_reuse_size(entry, xname='bx', yname='by', zname='bz'), entry)
#  if(reuse_size < resource[['rcapacity']]) {
#    print_if_2d('2d entry has reuse within capacity', reuse_size, resource[['rid']], resource[['rcapacity']], entry=entry)
#    return(TRUE)
#  }
#  print(paste('ts:', ts_size, sep=' '))
#  print(paste('reuse :', reuse_size, sep=' '))
  FALSE
}

resource_pred <- function(resource, entry) {
  if(entry[['Dim']] == 1) {
    return(check_bottleneck_generic(resource, entry, get_ts_size_1d, get_reuse_1d))
  } 
  else if(entry[['Dim']] == 2) {
    return(check_bottleneck_generic(resource, entry, get_ts_size_2d, get_reuse_2d))
  } 
  else if(entry[['Dim']] == 3) {
    return(check_bottleneck_generic(resource, entry, get_ts_size_3d, get_reuse_3d))
  } 
  else {
    stop("unrecognized dimension")
  }
}

get_reuse_distance_direct <- function(dim, oa, to, xmax, ymax, fpsize) {
  if(dim == 1) {
    v <- 1 * to
  } 
  else if(dim == 2) {
    v <-xmax * oa-1 *to
  } 
  else {
    v <- xmax * ymax * (oa*2-1) *to
  }

  v * fpsize
}

get_reuse_distance <- function(entry, xmax_name='xmax', ymax_name='ymax') {
  if(entry[['Dim']] == 1) {
    v <- 1 * entry[['To']]
  } 
  else if(entry[['Dim']] == 2) {
    v <- entry[[xmax_name]] * (entry[['Oa']]-1) * entry[['To']]
  } 
  else {
    v <- entry[[xmax_name]] * entry[[ymax_name]] * (entry[['Oa']]*2-1) * entry[['To']]
  }

  bytes(v, entry)
}

reuse_resource_predicate <- function(resource, entry) {
  # Does the shape allow data reuse?
  rd_size <- bytes(get_reuse_distance(entry, xmax_name='xmax', ymax_name='ymax'), entry)

  if(rd_size < resource[['rcapacity']]) {
    #print_if_2d('entry has ts within capacity', entry[['Dim']], ts_size, resource[['rid']], resource[['rcapacity']], entry=entry)
    return(TRUE)
  }

  FALSE
}


get_total_size <- function(entry, xmax_name='xmax', ymax_name='ymax', zmax_name='zmax') {
  bytes(get_ts_size(entry, xmax_name, ymax_name, zmax_name), entry)
}

get_cpi_slowdown <- function(entry) {
  if(entry[['divsLup']] > 0) {
    return(0.5)
  }
  1
}

pick <- function(entry, df2, predicate, ..., default=df2[4,]) {
#  return(df2[4,])
  for (i in seq_along(1:nrow(df2))) {
    row = df2[i,]
    if(predicate(row,entry,...))
      return(row)
  }
  default
}

set_blocksize <- function(entry) {

  if (entry[['optId']] == "orig") {
    tibble(
           bt = 1,
           bx = entry[['xmax']], 
           by = entry[['ymax']], 
           bz = entry[['zmax']]
           )
  }
  else if(entry[['Dim']] == 1) {
    tibble(
           bt = 128,
           bx = 512,
           by = -1,
           bz = -1
           )
  } 
  else if(entry[['Dim']] == 2) {
    tibble(
           bt = 64,
           bx = 32,
           by = 256,
           bz = -1
           )
  } 
  else if(entry[['Dim']] == 3) {
    tibble(
           bt = 64,
           bx = 32,
           by = 8,
           bz = 256
           )
  } 
  else {
    stop("unrecognized dimension")
  }
}
