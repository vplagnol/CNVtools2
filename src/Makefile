RDIST=/software/R-2.6.0/

INCLUDES= -I. -I$(RDIST)/include -I/usr/lib64/R/include -I/usr/share/R/include/
LIBS = -lg2c

PROF= #-pg
DEB=  -g
OPT = -O 

CPPFLAGS= -pthread -fPIC -Wall $(OPT) $(PROF) $(DEB) $(INCLUDES)

fitmodel.so: fitmodel.h 
	R CMD SHLIB fitmodel.cpp routine.cpp glm_test.cpp mla.cpp zeroin.c -o CNVtools.so

fitmodel.o : fitmodel.cpp fitmodel.h 
routine.o : routine.cpp fitmodel.h
glm_test.o : glm_test.cpp fitmodel.h
mla.o : mla.cpp fitmodel.h
zeroin.o : zeroin.c fitmodel.h

clean:
	rm -f CNVtools.so *.o test

listing:	
	nm -g CNVtools.so
