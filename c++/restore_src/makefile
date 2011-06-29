COMPILE=icpc
FLAGS=-g -O3
LINKS=-L/opt/netCDF/lib -lnetcdf -lm
LIBS=-L/opt/netCDF/lib
INCS=-I/opt/netCDF/include

OBJ= timers.o iounits.o $(GRID) $(CONSERV) remap_conserv.o $(REMAP) remap_write.o $(IO) $(UTILS) $(NEIGHBOR)
GRID=grids.o
IO=io.o namelist.o
UTILS= nc_error.o iounits.o utils.o 
REMAP=remap.o remap_vars.o 
CONSERV= intersection.o line_integral.o store_link_cnsrv.o
NEIGHBOR= neighborhood.o gradient.o

all: sc2rip
sc2rip: $(OBJ)
	$(COMPILE) $(FLAGS) sc2rip.cxx $(OBJ) -o sc2rip $(INCS) $(LINKS)
neighborhood.o: remap_conserv.o
	$(COMPILE) $(FLAGS) $(INCS) neighborhood.cxx -c
gradient.o:
	$(COMPILE) $(FLAGS) gradient.cxx -c
remap_conserv.o: $(GRID) $(REMAP) $(CONSERV)
	$(COMPILE) $(FLAGS) remap_conserv.cxx $(INCS) $(LIBS) -c
remap_write.o: $(GRID) $(UTILS) $(REMAP)
	$(COMPILE) $(FLAGS) remap_write.cxx $(INCS) $(LIBS) -c
remap.o:
	$(COMPILE) $(FLAGS) remap.cxx $(INCS) $(LIBS) -c
remap_vars.o: $(GRID) $(UTILS)
	$(COMPILE) $(FLAGS) remap_vars.cxx $(INCS) $(LIBS) -c
io.o:
	$(COMPILE) $(FLAGS) io.cxx -c
intersection.o: utils.o
	$(COMPILE) $(FLAGS) intersection.cxx -lm -c
line_integral.o: utils.o
	$(COMPILE) $(FLAGS) line_integral.cxx $(INCS) $(LIBS) -c
store_link_cnsrv.o: remap_vars.o
	$(COMPILE) $(FLAGS) store_link_cnsrv.cxx $(INCS) $(LIBS) -c
namelist.o:
	$(COMPILE) $(FLAGS) namelist.cxx -c
utils.o:
	$(COMPILE) $(FLAGS) utils.cxx -c
grids.o:
	$(COMPILE) $(FLAGS) grids.cxx $(INCS) $(LIBS) -c
nc_error.o:
	$(COMPILE) $(FLAGS) nc_error.cxx $(INCS) $(LIBS) -c
iounits.o: 
	$(COMPILE) $(FLAGS) iounits.cxx -c
timers.o:
	$(COMPILE) $(FLAGS) timers.cxx -c
clean:
	rm -f *.o *.out interp_remap.nc
