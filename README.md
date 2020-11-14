# stencil_optimizer
Automated Stencil Optimizer Tool

Supplementary tool for the dissertation **Supporting Performance in Scientific Applications with Energy and Resilience Constraints from Modern Architectures**

Required libraries/tools for use:
 * R > 3.6 
 * PLUTO Polyhedral Compiler: https://github.com/bondhugula/pluto
 * POET: https://github.com/qingyi-yan/POET

 _Note: This tool is very preliminary and has only been tested on our examples in the respective directory_

 Current working examples:
  * The Wave equation stencil in subdirectory examples/wave3d\_4oa
  * The Heat equation stencil in subdirectory examples/heat1d\_2oa

To test an example, repeat the steps for a given {case\_study}, where each step must be performed in the root diretory of this tool

```
cp examples/{case_study}/* .
Rscript predict_opt.R  --in {case_study}.c
gcc out.c -O3 -fopenmp -march -native -o opt_{case_study}.exe
```

_Note: Compilation of the output optimized stencil code is delegated to the user, who may use their preferred C compiler_
 
**Full Example**

```
cp example/heat1d_2oa/* .
Rscript predict_opt.R  --in heat1d_2oa.c
gcc out.c -O3 -fopenmp -march=native -lm -o opt_heat1d_2oa.exe
./opt_heat1d_2oa 16000000 0 0 64
```
