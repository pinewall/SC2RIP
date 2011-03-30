COMPILE=icpc
FLAGS=-g -O3
LINKS=-L/opt/netCDF/lib -lnetcdf -lm
LIBS=-L/opt/netCDF/lib
INCS=-I/opt/netCDF/include

OBJ= timers.o iounits.o nc_error.o grids.o utils.o namelist.o remap_vars.o remap.o
GRID= nc_error.o iounits.o utils.o namelist.o
all: sc2rip
sc2rip: $(OBJ)
	$(COMPILE) $(FLAGS) sc2rip.cxx $(OBJ) -o sc2rip $(INCS) $(LINKS)
remap.o:
	$(COMPILE) $(FLAGS) remap.cxx $(INCS) $(LIBS) -c
remap_vars.o: grids.o utils.o
	$(COMPILE) $(FLAGS) remap_vars.cxx $(INCS) $(LIBS) -c
namelist.o:
	$(COMPILE) $(FLAGS) namelist.cxx -c
utils.o:
	$(COMPILE) $(FLAGS) utils.cxx -c
grids.o: $(GRID)
	$(COMPILE) $(FLAGS) grids.cxx $(INCS) $(LIBS) -c
nc_error.o:
	$(COMPILE) $(FLAGS) nc_error.cxx $(INCS) $(LIBS) -c
iounits.o: 
	$(COMPILE) $(FLAGS) iounits.cxx -c
timers.o:
	$(COMPILE) $(FLAGS) timers.cxx -c
clean:
	rm *.o sc2rip
