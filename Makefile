# GNU Makefile for PCS (Using MSVC Compiler)
CC = cl
CXX = cl

CFLAGS = /EHsc /O2 /I. /Dkiss_fft_scalar=int16_t /nologo
LDFLAGS = /nologo

# Versioning
# -- Metadata --
# Note: Using your specific quoting style as requested
CURRENT_HASH := '"$(shell git rev-parse HEAD)"'
CURRENT_DATE := '"$(shell date /t)"'
CURRENT_NAME := '"pcs"'

# Target and Objects
TARGET = pcs.exe
# Ensure AcqUtils.obj is included in the link stage
OBJS = pcs.obj NCO.obj g2init.obj PCSEngine.obj kiss_fft.obj AcqUtils.obj

# Simplified escaping for MSVC
DEFS = /DGIT_HASH=$(CURRENT_HASH) /DBUILD_DATE=$(CURRENT_DATE) \
	/DAPP_NAME=$(CURRENT_NAME)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) /Fe:$(TARGET) /link $(LDFLAGS) /DEBUG

# Rule for AcqUtils specifically (includes DEFS)
AcqUtils.obj: AcqUtils.cpp AcqUtils.hpp
	$(CXX) $(CFLAGS) $(DEFS) /c AcqUtils.cpp /Fo:AcqUtils.obj

# Generic rule for all other CPP files (added DEFS so pcs.cpp can see version info too)
%.obj: %.cpp
	$(CXX) $(CFLAGS) $(DEFS) /c $< /Fo:$@

# Generic rule for C files (kiss_fft.c)
%.obj: %.c
	$(CC) $(CFLAGS) /c $< /Fo:$@

clean:
	@del /f /q $(OBJS) $(TARGET) *.pdb *.ilk 2>nul