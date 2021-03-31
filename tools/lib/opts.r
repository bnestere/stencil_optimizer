stencil_clioptions = list(
#                   make_option(c("-f", "--file"), type="character", default="stencilbench.csv", 
#                               help="dataset file name", metavar="character"),
                   make_option(c("--to"), type="integer", default=NULL, 
                               help="temporal order", metavar="number"),
                   make_option(c("--oa"), type="integer", default=NULL, 
                               help="order of accuracy", metavar="number"),
                   make_option(c("--dim"), type="integer", default=NULL, 
                               help="dimension", metavar="number"),
                   make_option(c("--nt"), type="integer", default=NULL, 
                               help="number of threads", metavar="number"),
                   make_option(c("-e", "--equation"), type="character", default=NULL, 
                               help="equation to isolate results for", metavar="character"),
                   make_option(c("-E", "--no-equation"), type="character", default=NULL, 
                               help="equation to remove from results"),
                   make_option(c("--original"), type="logical", default=FALSE, 
                               help="include original code data in report", metavar="include for true", 
                               action="store_true"),
                   make_option(c("--tileopt"), type="logical", default=FALSE, 
                               help="include blocking optimizations", metavar="include for true", 
                               action="store_true"),
                   make_option(c("--paropt"), type="logical", default=FALSE, 
                               help="include parallelization optimizations", metavar="include for true", 
                               action="store_true"),
                   make_option(c("--comboopt"), type="logical", default=FALSE, 
                               help="include tiling+parallelization optimizations", metavar="include for true", 
                               action="store_true")
#                   make_option(c("-x", "--xaxis"), type="character", default="Dim", 
#                               help="variable for the x axis"),
#                   make_option(c("-y", "--yaxis"), type="character", default="MLUPS", 
#                               help="variable for the y axis"),
#                   make_option(c("-g", "--groupby"), type="character", default="optId", 
#                               help="field to group data by in output graph")
                   )
#opt_parser = OptionParser(option_list=option_list)
#opt = parse_args(opt_parser)
