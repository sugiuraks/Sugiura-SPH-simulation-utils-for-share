PS_PATH = -I ./src/ -I ./my-src-for-flaw-number-change -I ./

CC = mpicc
CFLAGS = -O3
#CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
#CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

CPPOBJS = $(patsubst %.cpp, %.o, $(wildcard ../my-src/*.cpp))
CPPOBJM = $(patsubst %.cpp, %.o, $(wildcard ./*.cpp))
CPPHDRS = $(wildcard *.h) $(wildcard ../my-src/*.h)
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

makeinitial:
	@$(CC) make-initial-from-relaxed-condition.c -o make-initial.out $(PS_PATH)
	@echo "make make-initial!"

makeinitialfromsphere:
	@$(CC) make-initial-from-static-sphere.c -o make-initial-from-static-sphere.out $(PS_PATH)
	@echo "make make-initial-from-static-sphere!"

makedifferentresolution:
	@$(CC) make-different-resolution-from-relaxed-condition.c -o make-different-resolution-from-relaxed-condition.out $(PS_PATH)
	@echo "make make-different-resolution-from-relaxed-condition!"

makedifferentscale:
	@$(CC) make-different-scale-from-relaxed-condition.c -o make-different-scale-from-relaxed-condition.out $(PS_PATH)
	@echo "make make-different-scale-from-relaxed-condition!"

putrotation:
	@$(CC) put-rotation.c -o put-rotation.out $(PS_PATH)
	@echo "make put-rotation!"

makevariousinitialsfromsphere:
	@$(CC) make-various-initials-from-static-sphere.c -o make-various-initials-from-static-sphere.out $(PS_PATH)
	@echo "make make-various-initials-from-static-sphere!"

changeflawnumber:
	@$(CC) change-flaw-number.c -o change-flaw-number.out $(PS_PATH)
	@echo "make change-flaw-number!"
