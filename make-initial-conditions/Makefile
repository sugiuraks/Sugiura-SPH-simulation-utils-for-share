PS_PATH = -I ./src/ -I ./my-src -I ./

CC = mpicxx
CFLAGS = -O3 
#CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
#CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

CPPOBJS = $(patsubst %.cpp, %.o, $(wildcard ./my-src/*.cpp))
CPPOBJM = $(patsubst %.cpp, %.o, $(wildcard ./*.cpp))
CPPHDRS = $(wildcard *.h) $(wildcard ./my-src/*.h)
PROGRAM = sph.out

.PHONY:	clean all

all:	$(CPPOBJS) $(CPPOBJM) $(CPPHDRS)
	@echo "Linking object files..."
	@$(CC) $(CFLAGS) $(WARNINGS) $(CPPOBJM) $(CPPOBJS) -o $(PROGRAM) $(LIBS) $(PS_PATH)
	@echo "Link Success! [$(PROGRAM)]"

%.o:	%.cpp $(CPPHDRS) Makefile
	@echo "Building $< ..."	
	@$(CC) -c $< -o $(patsubst %.cpp, %.o, $<) $(CFLAGS) $(WARNINGS) $(PS_PATH)
	@echo "[$< OK]"

clean:
	-rm *.out *.o ../my-src/*.o

convertfromcgstomks:
	@$(CC) convert-from-cgs-to-mks.c -o convert-from-cgs-to-mks.out $(PS_PATH)
	@echo "make convert-from-cgs-to-mks!"

makedifferentresolution:
	@$(CC) make-different-resolution-from-relaxed-condition.c -o make-different-resolution-from-relaxed-condition.out $(PS_PATH)
	@echo "make make-different-resolution-from-relaxed-condition!"

makedifferentscale:
	@$(CC) make-different-scale-from-relaxed-condition.c -o make-different-scale-from-relaxed-condition.out $(PS_PATH)
	@echo "make make-different-scale-from-relaxed-condition!"

makespherefromrelaxedcondition:
	@$(CC) make-sphere-from-relaxed-condition.c -o make-sphere-from-relaxed-condition.out $(PS_PATH)
	@echo "make make-sphere-from-relaxed-condition!"

makespherefromrelaxedconditionANEOS:
	@$(CC) make-sphere-from-relaxed-condition-for-ANEOS.c -o make-sphere-from-relaxed-condition-for-ANEOS.out $(PS_PATH)
	@echo "make make-sphere-from-relaxed-condition-for-ANEOS!"

makeimpactconditionfromtwospheres:
	@$(CC) make-impact-condition-from-two-spheres.c -o make-impact-condition-from-two-spheres.out $(PS_PATH)
	@echo "make make-impact-condition-from-two-spheres!"

makevariousinitials:
	@$(CC) make-various-impact-conditions-from-two-spheres.c -o make-various-initials.out $(PS_PATH)
	@echo "make make-various-initials!"
