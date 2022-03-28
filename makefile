# Makefile: http://www.eng.hawaii.edu/Tutor/Make/

# Compiling flags
CFLAGS = -O3 -Wno-unused-result

# Linking flags
LFLAGS = -lm

# Compilers
CC = gcc

satpredict: satpredict.o sgdp4.o satutl.o deep.o ferror.o
	$(CC) -o satpredict satpredict.o sgdp4.o satutl.o deep.o ferror.o -lm

clean:
	rm -f *.o
	rm -f *~
