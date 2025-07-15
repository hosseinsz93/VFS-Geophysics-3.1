CC      = mpicxx

EXEC = VFS-Geophysics-3.1-mhk


LIBFLAG =       -lpthread -lrt -ldl -lstdc++  \
                -lpetsc -lHYPRE  -lgfortran -lflapack -lfblas

SOURCEC =       bcs.c bmv.c canopy.c compgeom.c ibm.c ibm_io.c init.c \
                main.c metrics.c poisson.c rhs.c timeadvancing.c \
                timeadvancing1.c variables.c fsi.c implicitsolver.c\
                fsi_move.c solvers.c rhs2.c wallfunction.c \
                les.c k-omega.c distance.c level.c momentum.c \
                poisson_hypre.c rotor_model.c \
                convection_diffusion.c sediment_transport.c turbinecontrol.c turbinestructure.c wave.c\

OBJSC   = $(SOURCEC:%.c=%.o)

CPPFLAGS =	-DNDEBUG        -DTECIO=1 -O3      

test:   $(OBJSC)
	$(CC) -o $(EXEC) $(OBJSC) $(CPPFLAGS)  $(LIBFLAG)


VFS-Geophysics-3.1-mhk:   $(OBJSC)
	$(CC) -o VFS-Geophysics-3.1-mhk $(OBJSC) $(CPPFLAGS) $(LIBDIR) $(LIBFLAG)

data-3.1: data.o
	$(CC) -o data-3.1 data.o $(CPPFLAGS) $(LIBFLAG) -ltecio64

data-small-3.1: data-small.o
	$(CC) -o data-small-3.1 data-small.o $(CPPFLAGS) $(LIBDIR) $(LIBFLAG) tecio64.a

data05: data05.o
	$(CC) -o data05 data05.o $(CPPFLAGS) $(LIBDIR) $(LIBFLAG)

itfcsearch: itfcsearch.o
	$(CC) -o itfcsearch itfcsearch.o $(CPPFLAGS) $(LIBDIR) $(LIBFLAG)

shear: shear.o

##data1: data1.o
##$(CC) -o data1 data1.o $(CPPFLAGS) $(LIBDIR) $(LIBFLAG) tecio64.a

hill: hill.o
	$(CC) -o hill hill.o $(CPPFLAGS) $(LIBDIR) $(LIBFLAG)

#xyz: xyz2Plot3d.o
#      $(CC) -o xyz xyz2Plot3d.o $(CPPFLAGS) $(LIBDIR) $(LIBFLAG)

clean:
	-rm -f *.o
	-rm -f data-3.1
	-rm -f VFS-Geophysics-3.1-mhk*

turbine: TurbineLoc.o
	$(CC) -o TurbineLoc TurbineLoc.o $(CPPFLAGS) $(LIBDIR) $(LIBFLAG)
terrain: terrain.o
	$(CC) -o terrain terrain.o $(CPPFLAGS) $(LIBDIR) $(LIBFLAG)
Tec2UCD: OSLTec2UCD.o
	$(CC) -o Tec2UCD OSLTec2UCD.o $(CPPFLAGS) $(LIBDIR) $(LIBFLAG)

#GenerateInflow: GenerateInflow__.o
#       $(CC) -o GenerateInflow GenerateInflow__.o $(CPPFLAGS) $(LIBDIR) $(LIBFLAG) tecio64.a

#ReadKevin_3dhilldat: ReadKevin_3dhilldat.o
#       $(CC) -o ReadKevin_3dhilldat ReadKevin_3dhilldat.o $(CPPFLAGS) $(LIBDIR) $(LIBFLAG) tecio64.a

#UCD2Plot3d: ucd2Plot3d_OSL.o
#       $(CC) -o UCD2Plot3d ucd2Plot3d_OSL.o $(CPPFLAGS) $(LIBDIR) $(LIBFLAG) tecio64.a

#GenInflow2: GenInflow2.o
#       $(CC) -o GenInflow2 GenInflow2.o $(CPPFLAGS) $(LIBDIR) $(LIBFLAG) tecio64.a


#SPAnal: SPAnal.o
#       $(CC) -o SPAnal SPAnal.o $(CPPFLAGS) $(LIBDIR) $(LIBFLAG) tecio64.a

#GenInflow3: GenInflow3.o
#       $(CC) -o GenInflow3 GenInflow3.o $(CPPFLAGS) $(LIBDIR) $(LIBFLAG) tecio64.a
