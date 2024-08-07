device=gpu

#configure machine
include ../../config/default.mk
include ../../config/version.mk
include ../../config/*.mk
include ../../config/devices/devices.mk

INCLUDE+= -I../         # other src libraries
INCLUDE+= -I../../inc   # other project libraries

all: thermal thermal_hpc thermal_mpi thermaldiag

thermaldiag: thermaldiag.cpp thermaldiag.h
	$(CC) $(OPT) $(CFLAGS) $< -o $@ $(INCLUDE) $(LIBS) $(JSONLIB) -g $(VERSION_FLAGS)

thermal: thermal.cpp thermaldiag.h thermal.h init.h parameters.h init_from_file.h
	$(CC) $(OPT) $(CFLAGS) $< -o $@ $(INCLUDE) $(GLFLAGS) $(LIBS) $(JSONLIB) -g $(VERSION_FLAGS)

thermal_hpc: thermal.cpp thermaldiag.h thermal.h init.h parameters.h init_from_file.h
	$(CC) -g $(OPT) $(CFLAGS) $< -o $@ $(INCLUDE) $(LIBS) $(JSONLIB) -DWITHOUT_GLFW $(VERSION_FLAGS)

thermal_mpi: thermal.cpp thermaldiag.h thermal.h init.h parameters.h init_from_file.h
	$(MPICC) $(OPT) $(MPICFLAGS) $< -o $@ $(INCLUDE) $(LIBS) $(JSONLIB) -DWITH_MPI -DWITHOUT_GLFW $(VERSION_FLAGS)

.PHONY: clean doc

doc:
	pdflatex -shell-escape thermal.tex
	bibtex thermal.aux;
	pdflatex -shell-escape thermal.tex
	pdflatex -shell-escape thermal.tex

clean:
	rm -f thermal thermal_hpc thermal_mpi thermaldiag thermal.bbl thermal.blg thermal.aux thermal.log thermal.out thermal.toc
