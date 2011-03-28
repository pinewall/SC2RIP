COMPILE=icpc
FLAGS=-g -O3
LINKS=-L/opt/netCDF/lib -lnetcdf -lm
LIBS=-L/opt/netCDF/lib
INCS=-I/opt/netCDF/include

all: test
test: test.cxx timers.o iounits.o nc_error.o grids.o scrip.o
	$(COMPILE) $(FLAGS) test.cxx timers.o iounits.o nc_error.o grids.o scrip.o utils.o -o test $(INCS) $(LINKS)
utils.o:
	$(COMPILE) $(FLAGS) utils.cxx -c
scrip.o:
	$(COMPILE) $(FLAGS) scrip.cxx $(INCS) -c
grids.o: nc_error.o iounits.o utils.o
	$(COMPILE) $(FLAGS) grids.cxx $(INCS) $(LIBS) -c
nc_error.o:
	$(COMPILE) $(FLAGS) nc_error.cxx $(INCS) $(LIBS) -c
iounits.o: 
	$(COMPILE) $(FLAGS) iounits.cxx -c
timers.o:
	$(COMPILE) $(FLAGS) timers.cxx -c
clean:
	rm *.o test
