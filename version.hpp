#pragma once
#include <cstdio>
#include <cstring>

#ifndef MAJOR_VERSION
#define MAJOR_VERSION 0
#endif
#ifndef MINOR_VERSION
#define MINOR_VERSION 1
#endif
#ifndef PATCH_VERSION
#define PATCH_VERSION 0
#endif
#ifndef CURRENT_HASH
#define CURRENT_HASH "unknown"
#endif
#ifndef BUILD_DATE
#define BUILD_DATE "unknown"
#endif
#ifndef APP_NAME
#define APP_NAME "unknown"
#endif

typedef struct {
    unsigned short int Major, Minor, Patch;
    unsigned int GitTag;
    char GitCI[41];
    char BuildDate[17];
    char Name[13];
} SWV;

class Version {
private:
    SWV v; // Store the structure as member data

public:
    // 1. The Default Constructor: Initializes the data immediately upon creation
    Version() {
        v.Major = MAJOR_VERSION;
        v.Minor = MINOR_VERSION;
        v.Patch = PATCH_VERSION;
        v.GitTag = 0;
        if (strcmp(CURRENT_HASH, "unknown") != 0) {
            sscanf(CURRENT_HASH, "%8x", &v.GitTag);
        }
        snprintf(v.GitCI, sizeof(v.GitCI), "%s", CURRENT_HASH);
        snprintf(v.BuildDate, sizeof(v.BuildDate), "%s", BUILD_DATE);
        snprintf(v.Name, sizeof(v.Name), "%s", APP_NAME);
    }

    // No destructor needed!

    // Getter to fetch the underlying structure
    SWV data() const { 
        return v; 
    }

    // Prints the instance data
    void print(FILE* stream = stdout) const {
        fprintf(stream, "%s v%d.%d.%d Built:%s GitCI:%s\n",
                v.Name, v.Major, v.Minor, v.Patch, v.BuildDate, v.GitCI);
    }
};