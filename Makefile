


all: sembDip_3d.exe

sembDip_3d.exe: sembDip_3d.o segyIO_class.o
	gcc -O2 -o  sembDip_3d.exe sembDip_3d.o segyIO_class.o -lm -fopenmp

sembDip_3d.o: sembDip_3d.c
	gcc -O2 -c sembDip_3d.c -lm -fopenmp

segyIO_class.o: segyIO_class.c
	gcc -O2 -c segyIO_class.c 


