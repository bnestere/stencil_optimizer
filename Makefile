include gmtt.mk

PCG=pcg
PCGFLAGS=-L../..
DBG=

CC=gcc
CFLAGS=-O3
LDFLAGS:=-lm -lstencil_config -lacc_timer -lbenchtool
LDPARFLAGS:=-fopenmp
REF_TEST_DIR=#!date +%Y_%m_%d:%H_%M
PLUTO_HOME=/home/brandon/opt_devel/pluto/target
PLC=$(PLUTO_HOME)/bin/polycc
PLCFLAGS=--pet
PLCPARFLAGS=--parallel
PLCTILEFLAGS=--noparallel --tile
#PLCTILEFLAGS=--tile
PLCALLFLAGS=--parallel --tile
PLCFULLDIAMONDFLAGS=--noparallel --full-diamond-tile
PLCFULLDIAMONDPARFLAGS=--full-diamond-tile $(PLCPARFLAGS)
PLCNODIAMONDFLAGS=--noparallel --nodiamond-tile --tile 
PLCNODIAMONDPARFLAGS=--nodiamond-tile --tile $(PLCPARFLAGS)
PLCRECTSPATIALPARFLAGS=$(PLCNODIAMONDPARFLAGS)
PLCRECTSPATIALL2PARFLAGS=$(PLCRECTSPATIALPARFLAGS) --second-level-tile
PLCTWOTILEFLAGS=--noparallel --second-level-tile
PLCTWOTILEPARFLAGS=--second-level-tile $(PLCPARFLAGS)
PLCRECTL2PARFLAGS=$(PLCNODIAMONDPARFLAGS) --second-level-tile
PLCFULLDL2PARFLAGS=$(PLCFULLDIAMONDPARFLAGS) --second-level-tile 
POCHOIRFLAGS=-O3 -DNDEBUG -funroll-loops -xHOST -fno-alias -fno-fnalias -fp-model precise -std=c++0x

PLCNONEFLAGS:=--nointratileopt --noparallel --notile --nodiamond-tile

EXT_ORIG=orig
EXT_TILE=tile
EXT_PAR=par
EXT_TWOTILE=twotile
EXT_TWOTILEPAR=twotilepar
EXT_FULLDIAMOND=fulldiamond
EXT_FULLDIAMONDPAR=fulldiamondpar
EXT_NODIAMOND=nodiamond
EXT_NODIAMONDPAR=nodiamondpar
EXT_RECTSPATIALPAR=rectspatialpar
EXT_RECTSPATIALL2PAR=rectspatiall2par
EXT_RECTL2PAR=rectl2par
EXT_FULLDL2PAR=fulldl2par
EXT_ALLOPT=all
EXT_POCHOIR=pochoir

EQ_HEAT:=heat
EQ_WAVE:=wave
EQ_LAPLA:=lapla
EQ_LAPLAINV:=laplainv
EQ_ALLENCAHN:=allencahn
EQ_ISO:=iso
EQ_ANISO:=aniso
EQ_ACINV:=acinv
EQ_WAVEMOD:=wavemod
EQ_GENERIC:=generic
EQ_POISSON:=poisson

TEST_MODE=0
RUN_ID=1

OUTPUT_DIR=target
STAGE_DIR=stage
STAGE_PATH=$(OUTPUT_DIR)/$(STAGE_DIR)
TEST_DIR=test_results
TEST_PATH=$(OUTPUT_DIR)/$(TEST_DIR)/$(CC)/$(EQ)/$(OPT)/$(RUN_ID)

#SRC=$(EQ)$(D)d_$(OA)oa
SRC=out

###################################
#### Compiler Parameters ##########
###################################
ifeq ($(CC),gcc)
CFLAGS:=$(CFLAGS) -march=native -mtune=native
#CFLAGS:=$(CFLAGS) -march=native -mtune=native -fopt-info-vec-optimized -mavx -msse2
endif

ifeq ($(CC),icc)
CFLAGS:=$(CFLAGS) -xHost -ansi-alias -ipo -fp-model precise
LDPARFLAGS=-qopenmp
endif

###################################
#### Dimensionality Parameters ####
###################################
T_MAX=250
D1_X = 3200000
D2_X = 5600
D2_Y = 5600
D3_X = 315
D3_Y = 315
D3_Z = 315

ifeq ($(D), 1)
#DIM_PARAMS=16000 0 0 $(T_MAX)
DIM_PARAMS=$(D1_X) 0 0 $(T_MAX)
#DIM_PARAMS=6400000 0 0 $(T_MAX)
endif

ifeq ($(D), 2)
#DIM_PARAMS=100 100 0 $(T_MAX)
DIM_PARAMS=$(D2_X) $(D2_Y) 0 $(T_MAX)
#DIM_PARAMS=8000 8000 0 $(T_MAX)
endif

ifeq ($(D), 3)
#DIM_PARAMS=20 20 20 $(T_MAX)
DIM_PARAMS=$(D3_X) $(D3_Y) $(D3_Z) $(T_MAX)
#DIM_PARAMS=400 400 400 $(T_MAX)
endif

###################################
#### Tiling Overrides  ############
###################################
#ifeq ($(EQ),$(EQ_WAVE))
#TILE_INPUT:=$(TILE_INPUT).2to
#endif
#ifeq ($(OPT),$(EXT_RECTSPATIALPAR))
#TILE_INPUT:=$(TILE_INPUT).spatial
#endif
#ifeq ($(OPT),$(EXT_RECTSPATIALL2PAR))
#TILE_INPUT:=$(TILE_INPUT).spatial
#endif
#TILE_INPUT:=$(TILE_INPUT).sizes

TILE_INPUT:= tile.$(D)d
TILE_OPT_SUFFIX:=

# partial diamond
ifeq ($(OPT),$(EXT_ALLOPT))
TILE_OPT_SUFFIX:=.pd
endif
ifeq ($(OPT),$(EXT_TILE))
TILE_OPT_SUFFIX:=.pd
endif

# full diamond
ifeq ($(OPT),$(EXT_FULLDIAMOND))
TILE_OPT_SUFFIX:=.fd
endif
ifeq ($(OPT),$(EXT_FULLDIAMONDPAR))
TILE_OPT_SUFFIX:=.fd
endif

# rect / no diamond tiling
ifeq ($(OPT),$(EXT_NODIAMOND))
TILE_OPT_SUFFIX:=.rect
endif
ifeq ($(OPT),$(EXT_NODIAMONDPAR))
TILE_OPT_SUFFIX:=.rect
endif

TILE_INPUT:=$(TILE_INPUT)$(TILE_OPT_SUFFIX).sizes


# rule names and file names cause conflict in target names, so we don't want any implicit rules
.SUFFIXES:

prep_parallel:
	$(eval LDFLAGS :=$(LDFLAGS) $(LDPARFLAGS))

#prep_tile: tile.$(D)d.sizes
# do nothing because tiling should be taken care of by predict opt
prep_tile:
#prep_tile: tile_sizes/$(TILE_INPUT)
#	cp tile_sizes/$(TILE_INPUT) tile.sizes

untile: 
	rm tile.sizes

prepare:
	mkdir -p $(STAGE_PATH)
	cp tile_sizes/$(TILE_INPUT) tile.sizes
	cp unopt/$(SRC).c $(SRC).c 
	cp exe_configs/$(SRC).json exe_config.json
	Rscript predict_opt.R --in $(SRC).c
	mv out.c $(STAGE_PATH)/$(SRC).c

clean:
	rm -rf $(OUTPUT_DIR)

gen_id=$(dim)$(oa)$(to)$(div)$(mi)$(T_MAX)$(D3_X)

$(SRC).$(EXT_ORIG).c: $(SRC).c prepare
	cp $(SRC).c $(STAGE_PATH)/$(SRC).$(OPT).c

$(SRC).$(EXT_TILE).c: $(SRC).c prepare prep_tile
	$(PLC) $(SRC).c $(PLCTILEFLAGS) $(PLCFLAGS) -o $(STAGE_PATH)/$(SRC).$(OPT).c

$(SRC).$(EXT_NODIAMOND).c: $(SRC).c prepare prep_tile
	$(PLC) $(SRC).c $(PLCNODIAMONDFLAGS) $(PLCFLAGS) -o $(STAGE_PATH)/$(SRC).$(OPT).c

$(SRC).$(EXT_NODIAMONDPAR).c: $(SRC).c prepare prep_tile prep_parallel
	$(PLC) $(SRC).c $(PLCNODIAMONDPARFLAGS) $(PLCFLAGS) -o $(STAGE_PATH)/$(SRC).$(OPT).c

$(SRC).$(EXT_FULLDIAMOND).c: $(SRC).c prepare prep_tile
	$(PLC) $(SRC).c $(PLCFULLDIAMONDFLAGS) $(PLCFLAGS) -o $(STAGE_PATH)/$(SRC).$(OPT).c
	
$(SRC).$(EXT_FULLDIAMONDPAR).c: $(SRC).c prepare prep_tile prep_parallel
	$(PLC) $(SRC).c $(PLCFULLDIAMONDPARFLAGS) $(PLCFLAGS) -o $(STAGE_PATH)/$(SRC).$(OPT).c

$(SRC).$(EXT_TWOTILE).c: $(SRC).c prepare prep_tile
	$(PLC) $(SRC).c $(PLCTWOTILEFLAGS) $(PLCFLAGS) -o $(STAGE_PATH)/$(SRC).$(OPT).c

$(SRC).$(EXT_TWOTILEPAR).c: $(SRC).c prepare prep_tile prep_parallel
	$(PLC) $(SRC).c $(PLCTWOTILEPARFLAGS) $(PLCFLAGS) -o $(STAGE_PATH)/$(SRC).$(OPT).c

$(SRC).$(EXT_RECTL2PAR).c: $(SRC).c prepare prep_tile prep_parallel
	$(PLC) $(SRC).c $(PLCRECTL2PARFLAGS) $(PLCFLAGS) -o $(STAGE_PATH)/$(SRC).$(OPT).c

$(SRC).$(EXT_FULLDL2PAR).c: $(SRC).c prepare prep_tile prep_parallel
	$(PLC) $(SRC).c $(PLCFULLDL2PARFLAGS) $(PLCFLAGS) -o $(STAGE_PATH)/$(SRC).$(OPT).c

$(SRC).$(EXT_RECTSPATIALPAR).c: $(SRC).c prepare untile prep_parallel
	$(PLC) $(SRC).c $(PLCRECTSPATIALPARFLAGS) $(PLCFLAGS) -o $(STAGE_PATH)/$(SRC).$(OPT).c

$(SRC).$(EXT_RECTSPATIALL2PAR).c: $(SRC).c prepare prep_tile prep_parallel
	$(PLC) $(SRC).c $(PLCRECTSPATIALL2PARFLAGS) $(PLCFLAGS) -o $(STAGE_PATH)/$(SRC).$(OPT).c

$(SRC).$(EXT_PAR).c: $(SRC).c prepare prep_parallel
	$(PLC) $(SRC).c $(PLCPARFLAGS) $(PLCFLAGS) -o $(STAGE_PATH)/$(SRC).$(OPT).c

$(SRC).$(EXT_ALLOPT).c: $(SRC).c prepare prep_parallel prep_tile
	$(PLC) $(SRC).c $(PLCALLFLAGS) $(PLCFLAGS) -o $(STAGE_PATH)/$(SRC).$(OPT).c

$(SRC).$(EXT_POCHOIR).c: $(SRC).cpp prepare
	pochoir -o $(STAGE_PATH)/$(SRC).$(OPT) $(POCHOIRFLAGS) $(SRC).$(EXT)

EQID=10
EQ=generic
ifeq ($(EQ), $(EQ_HEAT))
EQID=1
endif
ifeq ($(EQ), $(EQ_WAVE))
EQID=2
endif
ifeq ($(EQ), $(EQ_ALLENCAHN))
EQID=3
endif
ifeq ($(EQ), $(EQ_ISO))
EQID=4
endif
ifeq ($(EQ), $(EQ_ACINV))
EQID=5
endif
ifeq ($(EQ), $(EQ_WAVEMOD))
EQID=6
endif
ifeq ($(EQ), $(EQ_LAPLA))
EQID=7
endif
ifeq ($(EQ), $(EQ_LAPLAINV))
EQID=8
endif
ifeq ($(EQ), $(EQ_ANISO))
EQID=9
endif
ifeq ($(EQ), $(EQ_GENERIC))
EQID=10
endif
ifeq ($(EQ), $(EQ_POISSON))
EQID=11
endif

OPTID=0
ifeq ($(OPT), $(EXT_ORIG))
OPTID=1
endif
ifeq ($(OPT), $(EXT_NODIAMOND))
OPTID=2
endif
ifeq ($(OPT), $(EXT_TILE))
OPTID=3
endif
ifeq ($(OPT), $(EXT_FULLDIAMOND))
OPTID=4
endif
ifeq ($(OPT), $(EXT_PAR))
OPTID=5
endif
ifeq ($(OPT), $(EXT_NODIAMONDPAR))
OPTID=8
endif
ifeq ($(OPT), $(EXT_ALLOPT))
OPTID=9
endif
ifeq ($(OPT), $(EXT_FULLDIAMONDPAR))
OPTID=10
endif
ifeq ($(OPT), $(EXT_POCHOIR))
OPTID=11
endif


build:
	$(PCG) $(PCGFLAGS) $(DBG) -pdim=$(dim) -ptime=$(T_MAX) -poa=$(oa) -pto=$(to) -pdiv=$(div) -pmi=$(mi) ../../stencil-gen-pochoir.pt > outcpp

EXT=c

REPORT_PATH=$(OUTPUT_DIR)/casestudies.csv
N_THREADS=1
opt: export BENCHTOOL_EQID=$(EQID)
opt: export BENCHTOOL_RUNID=$(RUN_ID)
opt: export BENCHTOOL_OPTID=$(OPTID)
opt: export BENCHTOOL_FILE=$(REPORT_PATH)
opt: export BENCHTOOL_NTHREADS=$(N_THREADS)
opt: export OMP_NUM_THREADS=$(N_THREADS)
#opt: export CILK_NWORKERS=$(N_THREADS)

#opt: prepare $(SRC).$(OPT).c
opt: prepare 
	@echo "Working {dim:${dim}, oa:${oa}, to:${to}, ct:${coeftype}, div:${div}, mi:${mi}, nt:${N_THREADS}, rid:${RUN_ID}}"
	$(eval EXE_OUT=x$(SRC)_$(CC))
	$(CC) $(CFLAGS) $(STAGE_PATH)/$(SRC).c -o $(OUTPUT_DIR)/$(EXE_OUT) $(LDFLAGS); 
	@if [ "$(TEST_MODE)" = "1" ]; then\
		echo "INFO: Running: ./$(OUTPUT_DIR)/$(EXE_OUT) $(DIM_PARAMS) > $(TEST_PATH)/$(EXE_OUT);"; \
		mkdir -p $(TEST_PATH); \
		./$(OUTPUT_DIR)/$(EXE_OUT) $(DIM_PARAMS) > $(TEST_PATH)/$(EXE_OUT); \
	fi

#	$(CC) $(CFLAGS) $(STAGE_PATH)/$(SRC).$(OPT).c -o $(OUTPUT_DIR)/$(EXE_OUT) $(LDFLAGS)

stencil_id="generic"

# 8 threads for Poisson because it was pre-determinedvia the model from the 16 max
poisson:
	@$(MAKE) dim=2 D=2 oa=4 to=0 div=1 mi=1 CC=gcc N_THREADS=8 D2_X=2863 D2_Y=2863 T_MAX=256 SRC=poisson2d_4oa EQ=poisson OPT=all opt

run-poisson:
	@$(MAKE) TEST_MODE=1 poisson

wave:
	@$(MAKE) dim=2 D=2 oa=4 to=2 div=0 mi=0  CC=gcc N_THREADS=16 D1_X=1600000 D2_X=4096 D2_Y=4096 D3_X=256 D3_Y=256 D3_Z=256 T_MAX=128 SRC=wave2d_4oa EQ=wave OPT=fulldiamondpar opt

run-wave:
	@$(MAKE) TEST_MODE=1 wave

heat:
	@$(MAKE) dim=3 D=3 oa=2 to=1 div=0 mi=0 CC=gcc N_THREADS=16 D1_X=6400000 D2_X=8000 D2_Y=8000 D3_X=216 D3_Y=216 D3_Z=216 T_MAX=64 SRC=heat3d_2oa EQ=heat OPT=all opt

run-heat:
	@$(MAKE) TEST_MODE=1 heat

lapla:
	@$(MAKE) dim=1 D=1 oa=4 to=0 div=0 mi=0 CC=gcc D1_X=16384000 N_THREADS=8 T_MAX=32 SRC=lapla1d_4oa EQ=lapla OPT=nodiamondpar opt

run-lapla:
	@$(MAKE) TEST_MODE=1 lapla

allencahn:
	@$(MAKE) dim=3 D=3 oa=2 to=1 div=7 mi=1 CC=gcc N_THREADS=1 D1_X=64000000 D2_X=8000 D2_Y=8000 D3_X=161 D3_Y=161 D3_Z=161 T_MAX=32 SRC=allencahn3d_2oa EQ=allencahn EXT=cpp OPT=all opt

run-allencahn:
	@$(MAKE) TEST_MODE=1 allencahn

run-norms: 
	@$(MAKE) run-allencahn
	@$(MAKE) run-lapla
	@$(MAKE) run-heat
	@$(MAKE) run-poisson
	@$(MAKE) run-wave

run-norm-opts:
	@$(MAKE) EXT=c OPT=orig run-norms
#	@$(MAKE) EXT=c OPT=all run-norms
#	@$(MAKE) EXT=c OPT=fulldiamondpar run-norms
#	@$(MAKE) EXT=c OPT=nodiamond run-norms
#	@$(MAKE) EXT=c OPT=tile run-norms
#	@$(MAKE) EXT=c OPT=fulldiamond run-norms

run-pochoirs:
	@$(MAKE) OPT=pochoir EXT=cpp run-allencahn
	#@$(MAKE) OPT=pochoir EXT=cpp run-wave

run-all:
	@$(MAKE) run-norm-opts
	#@$(MAKE) run-poisson
#	@$(MAKE) run-pochoirs

run-all-n:
	@$(MAKE) RUN_ID=1 run-all
	@$(MAKE) RUN_ID=2 run-all
	@$(MAKE) RUN_ID=3 run-all
	@$(MAKE) RUN_ID=4 run-all

	#@$(MAKE) dim=2 D=2 oa=4 to=0 div=1 mi=1 OPT=fulldiamondpar CC=gcc N_THREADS=16 T_MAX=32 SRC=poisson2d_4oa EQ=poisson opt
