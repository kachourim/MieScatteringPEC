APPNAME := MiePEC
FLAGS   := -Wall -O3 -lpthread -fopenmp -lgfortran -m64

$(APPNAME):
	g++ -o $(APPNAME) main.cpp Mie_Functions.cpp libslatec.a $(FLAGS)
	strip -s $(APPNAME)
	mv $(APPNAME) ..

