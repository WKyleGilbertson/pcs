CC=CL
CURRENT_HASH := '"$(shell git rev-parse HEAD)"'
CURRENT_DATE := '"$(shell date /t)"'
CURRENT_NAME := '"pcs"'

pcs: pcs.cpp NCO.cpp G2INIT.cpp FftComplex.cpp
	$(CC) pcs.cpp NCO.cpp G2INIT.cpp PCSEngine.cpp FftComplex.cpp /EHsc \
	 /DCURRENT_HASH=$(CURRENT_HASH) /DCURRENT_DATE=$(CURRENT_DATE) \
	 /DCURRENT_NAME=$(CURRENT_NAME)
	del *.obj

clean:
	del *.obj
	del *.exe
	del *.pdb