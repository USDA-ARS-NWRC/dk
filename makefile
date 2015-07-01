ADDL_OPTIONS=-Wall -fopenmp
NETCDF_INC=-I/opt/local/include -DNDEBUG 
NETCDF_LIBS=-L/opt/local/lib -lnetcdf

dk : dk.o arcout.o array.o caldate.o dist.o getln.o\
     grassout.o index.o interp.o ipwout.o isleap.o krige.o lusolv.o\
     medfit.o netcdfout.o period1.o period2.o readcnfg.o readcsv.o readdata.o\
     readgrid.o sca_grid.o sreg.o storm1.o storm2.o\
     swe1.o swe2.o wyjdate.o zoneout.o
	gcc  -o dk $(ADDL_OPTIONS) $(NETCDF_INC) $(NETCDF_LIBS) dk.o arcout.o array.o caldate.o \
	dist.o getln.o grassout.o index.o interp.o ipwout.o \
	isleap.o krige.o lusolv.o medfit.o netcdfout.o period1.o period2.o readcnfg.o \
	readcsv.o readdata.o readgrid.o sca_grid.o sreg.o storm1.o \
	storm2.o swe1.o swe2.o wyjdate.o zoneout.o  -lm

dk.o : dk.c dk_m.h
	gcc $(ADDL_OPTIONS) -c dk.c 

arcout.o : arcout.c dk_x.h
	gcc -c $(ADDL_OPTIONS) arcout.c

array.o : array.c
	gcc -c $(ADDL_OPTIONS) array.c 

caldate.o : caldate.c dk_x.h
	gcc -c $(ADDL_OPTIONS) caldate.c

dist.o : dist.c
	gcc -c $(ADDL_OPTIONS) dist.c 

getln.o : getln.c
	gcc -c $(ADDL_OPTIONS) getln.c

grassout.o : grassout.c dk_x.h
	gcc -c $(ADDL_OPTIONS) grassout.c

index.o : index.c
	gcc -c $(ADDL_OPTIONS) index.c

interp.o : interp.c dk_x.h
	gcc -c $(ADDL_OPTIONS) interp.c

ipwout.o : ipwout.c dk_x.h
	gcc -c $(ADDL_OPTIONS) ipwout.c

isleap.o : isleap.c dk_x.h
	gcc -c $(ADDL_OPTIONS) isleap.c

krige.o : krige.c dk_x.h
	gcc -c $(ADDL_OPTIONS) krige.c

lusolv.o : lusolv.c
	gcc -c $(ADDL_OPTIONS) lusolv.c 

medfit.o : medfit.c
	gcc -c $(ADDL_OPTIONS) medfit.c
	
netcdfout.o : netcdfout.c
	gcc -c $(ADDL_OPTIONS) $(NETCDF_INC) $(NETCDF_LIBS) netcdfout.c
#	gcc -c $(ADDL_OPTIONS) `nc-config --cflags` netcdfout.c `nc-config --libs`

period1.o : period1.c dk_m.h dk_x.h
	gcc -c $(ADDL_OPTIONS) period1.c

period2.o : period2.c dk_x.h
	gcc -c $(ADDL_OPTIONS) period2.c

readcnfg.o : readcnfg.c dk_x.h dk_m.h
	gcc -c $(ADDL_OPTIONS) readcnfg.c

readcsv.o : readcsv.c dk_x.h
	gcc -c $(ADDL_OPTIONS) readcsv.c

readdata.o : readdata.c dk_x.h
	gcc -c $(ADDL_OPTIONS) readdata.c

readgrid.o : readgrid.c dk_x.h
	gcc -c $(ADDL_OPTIONS) readgrid.c

sca_grid.o : sca_grid.c dk_x.h
	gcc -c $(ADDL_OPTIONS) sca_grid.c

sreg.o : sreg.c
	gcc -c $(ADDL_OPTIONS) sreg.c 

storm1.o : storm1.c dk_x.h
	gcc -c $(ADDL_OPTIONS) storm1.c 

storm2.o : storm2.c dk_x.h
	gcc -c $(ADDL_OPTIONS) storm2.c 

swe1.o : swe1.c dk_x.h
	gcc -c $(ADDL_OPTIONS) swe1.c 

swe2.o : swe2.c dk_x.h
	gcc -c $(ADDL_OPTIONS) swe2.c

wyjdate.o : wyjdate.c dk_x.h
	gcc -c $(ADDL_OPTIONS) wyjdate.c

zoneout.o : zoneout.c dk_x.h
	gcc -c $(ADDL_OPTIONS) zoneout.c
