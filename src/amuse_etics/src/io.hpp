#pragma once

struct ParametersStruct {
    int N;
    std::string Filename;
    int Seed;
    Real dT1, dT2, Tcrit, ConstantStep;
    std::string Prefix;
    std::string OutputFormat;
    int DeviceID;
};

void ReadICsASCII(std::string Filename, int N, Particle **P_h, int *FileSnapshotNum, Real *FileTime);
void WriteSnapshotASCII(std::string Prefix, int SnapNumber, Particle *hostP, int N, Real T);
void ParseInput(int argc, char *argv[], ParametersStruct *Params);
#ifdef ETICS_HDF5
void ReadICsHDF5(std::string Filename, int N, Particle **P_h, int *FileSnapshotNum, Real *FileTime);
void WriteSnapshotHDF5(std::string Prefix, int SnapNumber, Particle *hostP, int N, Real T);
#endif
