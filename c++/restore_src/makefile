COMPILE=icpc
FLAGS=-g -O0
NETCDF_LIBS=-L/opt/netCDF/lib
NETCDF_INCS=-I/opt/netCDF/include
CSPARSE_INCS=-I/home/pine/Works/WangXiaoGe/SparseMatrix/CSparse/Include
CSPARSE_LIBS=-L/home/pine/Works/WangXiaoGe/SparseMatrix/CSparse/Lib
#LINKS=-lnetcdf -lm
#LINKS=-lnetcdf -lm -lcsparse
LINKS=-lnetcdf -lm /home/pine/Works/WangXiaoGe/SparseMatrix/CSparse/Lib/libcsparse.a

OBJ= timers.o iounits.o $(GRID) $(CONSERV) remap_conserv.o $(REMAP) remap_write.o $(IO) $(UTILS) $(NEIGHBOR)
GRID=grids.o
IO=io.o namelist.o
UTILS= nc_error.o iounits.o utils.o 
REMAP=remap.o remap_vars.o 
CONSERV= intersection.o line_integral.o store_link_cnsrv.o
NEIGHBOR= neighborhood.o gradient.o

all: sc2rip sc2rip_test
sc2rip: $(OBJ)
	$(COMPILE) $(FLAGS) $(NETCDF_INCS) $(CSPARSE_INCS) -o sc2rip sc2rip.cxx $(OBJ) $(NETCDF_LIBS) $(CSPARSE_LIBS) $(LINKS)
sc2rip_test: $(OBJ)
	$(COMPILE) $(FLAGS) $(NETCDF_INCS) -o sc2rip_test  sc2rip_test.cxx $(OBJ) $(NETCDF_LIBS) $(LINKS)
neighborhood.o: remap_conserv.o
	$(COMPILE) $(FLAGS) $(NETCDF_INCS) neighborhood.cxx -c
gradient.o:
	$(COMPILE) $(FLAGS) gradient.cxx -c
remap_conserv.o: $(GRID) $(REMAP) $(CONSERV)
	$(COMPILE) $(FLAGS) remap_conserv.cxx $(NETCDF_INCS) -c
remap_write.o: $(GRID) $(UTILS) $(REMAP)
	$(COMPILE) $(FLAGS) remap_write.cxx $(NETCDF_INCS) -c
remap.o:
	$(COMPILE) $(FLAGS) remap.cxx $(NETCDF_INCS) -c
remap_vars.o: $(GRID) $(UTILS)
	$(COMPILE) $(FLAGS) remap_vars.cxx $(NETCDF_INCS) -c
io.o:
	$(COMPILE) $(FLAGS) $(NETCDF_INCS) io.cxx -c
intersection.o: utils.o
	$(COMPILE) $(FLAGS) intersection.cxx $(NETCDF_INCS) -c
line_integral.o: utils.o
	$(COMPILE) $(FLAGS) line_integral.cxx $(NETCDF_INCS) -c
store_link_cnsrv.o: remap_vars.o
	$(COMPILE) $(FLAGS) store_link_cnsrv.cxx $(NETCDF_INCS) -c
namelist.o:
	$(COMPILE) $(FLAGS) namelist.cxx -c
utils.o:
	$(COMPILE) $(FLAGS) utils.cxx -c
grids.o:
	$(COMPILE) $(FLAGS) grids.cxx $(NETCDF_INCS) -c
nc_error.o:
	$(COMPILE) $(FLAGS) nc_error.cxx $(NETCDF_INCS) -c
iounits.o: 
	$(COMPILE) $(FLAGS) iounits.cxx -c
timers.o:
	$(COMPILE) $(FLAGS) timers.cxx -c
clean:
	rm -f *.o
distclean:
	rm -f *.o *.out interp_remap.nc sc2rip sc2rip_test
