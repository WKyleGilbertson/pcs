# GNU Makefile for PCS (Using MSVC Compiler)
CC = cl
CXX = cl

# Flags
# Note: /nologo keeps the output clean
CFLAGS = /EHsc /O2 /I. /Dkiss_fft_scalar=float /nologo
LDFLAGS = /nologo

# Versioning Information
# We use $(shell ...) for GNU Make to run external commands
CURRENT_HASH := "$(shell git rev-parse --short HEAD)"
CURRENT_DATE := "$(shell echo %date%)"
CURRENT_NAME := "pcs"

# Target and Objects
TARGET = pcs.exe
OBJS = pcs.obj NCO.obj g2init.obj PCSEngine.obj kiss_fft.obj

# Definitions for versioning
DEFS = /DCURRENT_HASH='$(CURRENT_HASH)' \
       /DCURRENT_DATE='$(CURRENT_DATE)' \
       /DCURRENT_NAME='$(CURRENT_NAME)'

# Default Rule
all: $(TARGET)

# Linker Step
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) /Fe:$(TARGET) /link $(LDFLAGS)

# Pattern Rule for C++ files
# $< is the source file, $@ is the target object
%.obj: %.cpp
	$(CXX) $(CFLAGS) $(DEFS) /c $< /Fo:$@

# Special Rule for C files (kiss_fft.c)
%.obj: %.c
	$(CC) $(CFLAGS) /c $< /Fo:$@

# Cleanup
clean:
	@del /f /q $(OBJS) $(TARGET) *.pdb *.ilk 2>nul