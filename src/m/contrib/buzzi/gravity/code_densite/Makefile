all:
		g++ density_rock_uneven.cpp  \
			-I$(ISSM_DIR)/externalpackages/gsl/install/include \
			-L$(ISSM_DIR)/externalpackages/gsl/install/lib -lgsl -lgslcblas \
		   -I$(ISSM_DIR)/externalpackages/mpich/install/include/ \
		   -L$(ISSM_DIR)/externalpackages/mpich/install/lib/ -lmpich -lpmpich -lmpl -o forward.exe 

forward:
	   mpirun -np 4 ./forward.exe
run -np 4 ./forward.execlean:
		rm ./forward.exe
