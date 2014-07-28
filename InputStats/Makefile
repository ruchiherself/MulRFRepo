OUTEXEC=InputStats

#For Windows: uncomment the next two lines 
#cpp=c++ -g -O -funroll-loops -Wno-long-long
#cc=gcc -O3 -fomit-frame-pointer -funroll-loops

#For Mac: uncomment the next two lines 
#cpp=c++ -g -O3 -funroll-loops -Wno-long-long -mmacosx-version-min=10.0
#cc=gcc -O3 -funroll-loops -mmacosx-version-min=10.0

#For Linux: uncomment the next two lines
#cpp=c++ -g -O3 -static
#cc=gcc -O3 

INCLUDE=-I./include

all: InputStats

InputStats: main.o 
	${cpp} main.o ${INCLUDE} -o ${OUTEXEC}

main.o: main.cpp Makefile 
	${cpp} ${INCLUDE} -c $<

clean:
	rm -f *.o *~ core ${OUTEXEC}



