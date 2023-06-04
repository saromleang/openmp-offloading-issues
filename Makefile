SHELL:=/bin/bash

.PHONY : clean

all : dgemm.f03.x dgemm-offload.f03.x

dgemm.f03.x : dgemm.f03
	source ${MODULESHOME}/init/bash; \
	module load cpe-22.12 cce/15.0.0 craype-accel-amd-gfx90a rocm/5.4.0; \
	module use /lustre/orion/world-shared/chm135/hipfort/crusher/cpe-22.12/modulefiles/cce/15.0/rocm5.4; \
	module load hipfort; \
	hipfc -sdefault64 -M878,3086 -lhipblas dgemm.f03 -o dgemm.f03.x

dgemm-offload.f03.x : dgemm-offload.f03
	source ${MODULESHOME}/init/bash; \
	module load cpe-22.12 cce/15.0.0 craype-accel-amd-gfx90a rocm/5.4.0; \
	module use /lustre/orion/world-shared/chm135/hipfort/crusher/cpe-22.12/modulefiles/cce/15.0/rocm5.4; \
	module load hipfort; \
	hipfc -homp -sdefault64 -M878,3086 -lhipblas dgemm-offload.f03 -o dgemm-offload.f03.x

clean : 
	rm -f dgemm.f03.x dgemm-offload.f03.x
