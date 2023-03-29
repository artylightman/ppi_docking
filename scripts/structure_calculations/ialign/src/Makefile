##############################
#### Fortran compiler
##############################

### Intel fortran
FC     = ifort
FFLAGS = -extend_source 132 -real-size 64 -fast
### debug purpose
#FFLAGS = -extend_source 132 -real-size 64 -g -check all

### PGF fortran
#FC	= pgf90
#FFLAGS	= -Mextend -Bstatic -fastsse -O3 -Mr8

###############################
#### C compiler
###############################

CC       = cc
CFLAGS   = -ansi -O6 -fomit-frame-pointer -Wall -W -Wcast-qual -Wpointer-arith -Wcast-align -pedantic -fno-schedule-insns -fschedule-insns2 -fstrict-aliasing -m64
CLINKERS = -lm -static 


#############
PROGRAMS = IS-score extint
IS_OBJS  = IS-score.o
EXT_OBJS = extint.o structures.o

all:	$(PROGRAMS)
 
#############

IS-score:	$(IS_OBJS)
	$(FC) $(FFLAGS) -o $@ $(IS_OBJS)

extint:		$(EXT_OBJS) structures.h
	$(CC) $(CFLAGS) -o $@ $(EXT_OBJS) $(CLINKERS)


#############
clean:
	rm -f *.o $(PROGRAMS)
#############


