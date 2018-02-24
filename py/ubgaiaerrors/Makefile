# Variables
EXEC=gaia_errors
OBJS= main.o gaia_errors_DR2.o
FORT = gfortran
FFLAGS = -O 

${EXEC} : ${OBJS}
	${FORT} -o ${EXEC} ${OBJS} ${FFLAGS}

.F.o :
	${FORT} ${FFLAGS} -c $<

clean:  
	rm -f *.o
