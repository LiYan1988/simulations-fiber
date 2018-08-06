# Windows MinGW (gcc):
EXT = .dll
MEX = C:/MATLAB7/bin/win32/mex.bat
INCLUDE = include
MEXOPTS = mexopts.bat
FFTWLIB = lib/libfftw3.a
FFTWLIBF = lib/libfftw3f.a

# Linux (gcc)
# EXT = .mexglx
# MEX = mex
# INCLUDE = include
# MEXOPTS = mexopts.sh
# FFTWLIB = libfftw3.a
# FFTWLIBF = libfftw3f.a

all: sspropc$(EXT) sspropcf$(EXT) sspropvc$(EXT) sspropvcf$(EXT)

double: sspropc$(EXT) sspropvc$(EXT)

single: sspropcf$(EXT) sspropvcf$(EXT)

clean:
	rm sspropc$(EXT) sspropcf$(EXT) sspropvc$(EXT) sspropvcf$(EXT)


sspropc$(EXT): sspropc.c
	$(MEX) -f $(MEXOPTS) -I$(INCLUDE) $? $(FFTWLIB) -output sspropc

sspropcf$(EXT): sspropc.c
	$(MEX) -f $(MEXOPTS) -I$(INCLUDE) -DSINGLEPREC $? $(FFTWLIBF) -output sspropcf

sspropvc$(EXT): sspropvc.c
	$(MEX) -f $(MEXOPTS) -I$(INCLUDE) $? $(FFTWLIB) -output sspropvc

sspropvcf$(EXT): sspropvc.c
	$(MEX) -f $(MEXOPTS) -I$(INCLUDE) -DSINGLEPREC $? $(FFTWLIBF) -output sspropvcf
