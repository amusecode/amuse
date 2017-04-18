#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <climits>
#include <fstream>
#include <sstream>
#include "common.hpp"
#include "io.hpp"

#ifndef NOBOOST
    #include <boost/property_tree/ptree.hpp>
    #include <boost/property_tree/ini_parser.hpp>
#endif

using namespace std;

void ReadICsASCII(string Filename, int N, Particle **P_h, int *FileSnapshotNum, Real *FileTime) {
    *P_h = new Particle[N];
    string Line;
    ifstream InputFile(Filename.c_str());
    int LineNum;

    int i = 0;
    if (InputFile.is_open()) {
        getline(InputFile, Line); // First line is the snapshot number.
        int Count = sscanf(Line.c_str(), "%d", FileSnapshotNum);
        if (Count != 1) {
            cerr << "Problem in input file, line 1." << endl;
            exit(1);
        }
        getline(InputFile, Line); // Second line is the expected number of particles.
        int ExpectedN;
        Count = sscanf(Line.c_str(), "%d", &ExpectedN); // We ignore it actually.
        if (Count != 1) {
            cerr << "Problem in input file, line 2." << endl;
            exit(1);
        }
        getline(InputFile, Line);
        double T0;
        Count = sscanf(Line.c_str(), "%lf", &T0);
        if (Count != 1) {
            cerr << "Problem in input file, line 3." << endl;
            exit(1);
        }
        *FileTime = (Real)T0;
        LineNum = 3;
        while (getline(InputFile, Line))
        {
            double m, x, y, z, vx, vy, vz;
            int id;
            Particle p;
            Count = sscanf(Line.c_str(), "%d %lf %lf %lf %lf %lf %lf %lf", &id, &m, &x, &y, &z, &vx, &vy, &vz);
            if (Count != 8) {
                cerr << "Problem in input file, line " << LineNum + 1 << "." << endl;
                exit(1);
            }
            p.ID = id;
            p.m = (Real)m;
            p.pos = vec3((Real)x, (Real)y, (Real)z);    // This is to ensure proper casting in the case of single prec.
            p.vel = vec3((Real)vx, (Real)vy, (Real)vz);
            p.Status = 0;
            p.CalculateR2();
            (*P_h)[i] = p;
            i++;
            LineNum++;
            if (i == N) break;
        }
        InputFile.close();
    } else {
        cerr << "Can't open file" << Filename << "." << endl;
        exit(1);
    }
    if (i < N) {
        cerr << "Was only able to read " << i << " particles while " << N << " were requested." << endl;
        exit(1);
    }
}

void WriteSnapshotASCII(string Prefix, int SnapNumber, Particle *P_h, int N, Real T) {
    char S[512];
    sprintf(S, "%s%04d.dat", Prefix.c_str(), SnapNumber);
    ofstream SnapshotFile;
    SnapshotFile.open(S);
    sprintf(S, "%06d\n", SnapNumber); SnapshotFile << S;
    sprintf(S, "%06d\n", N); SnapshotFile << S;
    sprintf(S, "%.16E\n", T); SnapshotFile << S;
    for (int i = 0; i < N; i++) {
        sprintf(S, "%06d%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e\n", P_h[i].ID, P_h[i].m, P_h[i].pos.x, P_h[i].pos.y, P_h[i].pos.z, P_h[i].vel.x, P_h[i].vel.y, P_h[i].vel.z); SnapshotFile << S;
    }
    SnapshotFile.close();
}

#define THROW_EXCEPTION(str, stat) {cerr << str << endl; exit((stat));}
void ParseInput(int argc, char *argv[], ParametersStruct *Params) {
    ParametersStruct P;
#ifndef NOBOOST
    if (argc < 2) THROW_EXCEPTION("Usage: etics [FILE]...", 1)
    ifstream TestFileObject(argv[1]);
    if (!TestFileObject) THROW_EXCEPTION("Problem reading file...", 1)
    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini(argv[1], pt); // what if file doesn't exist?

    P.N = pt.get<int>("N", 0);
    if (P.N <= 0) THROW_EXCEPTION("Could not read number of particles (N) from ini file.", 1)

    P.Tcrit = pt.get<Real>("Tcrit", -1);
    if (P.Tcrit < 0) THROW_EXCEPTION("Could not read finish time (Tcrit) from ini file.", 1)

    P.dT1 = pt.get<Real>("dT1", 0);
    if (P.dT1 <= 0) THROW_EXCEPTION("Could output time interval (dT1) from ini file.", 1)

    P.dT2 = pt.get<Real>("dT2", 0);
    if (P.dT2 <= 0) THROW_EXCEPTION("Could snapshot time interval (dT2) from ini file.", 1)

    P.ConstantStep = pt.get<Real>("StepSize", -1);
    if (P.ConstantStep < 0) THROW_EXCEPTION("Could not read step size (StepSize) from ini file.", 1)

    P.Filename = pt.get<string>("Filename", "\n");
    if (P.Filename == "\n") THROW_EXCEPTION("Could not read initial condition file name (Filename) from ini file.", 1)

    P.Prefix = pt.get<string>("Prefix", "");
    P.OutputFormat = pt.get<string>("OutputFormat", "ascii");
    P.DeviceID = pt.get<int>("device", -1);
    if (argc >= 3) { // If there is a argument after the file, it must be either --device or -d
        int Version1 = strncmp(argv[2], "--device", 8);
        int Version2 = strncmp(argv[2], "-d", 2);
        if ((Version1!=0) && (Version2!=0)) THROW_EXCEPTION("Commandline argument after input file must be device number.", 1)
        char *DeviceIDStr = strchr(argv[2], '=');
        if (DeviceIDStr!=NULL) DeviceIDStr = DeviceIDStr + 1;
        else {
            if (argc <= 3) THROW_EXCEPTION("Device number must follow.", 1)
            DeviceIDStr = argv[3];
        }
        P.DeviceID = atoi(DeviceIDStr);
        if ((P.DeviceID==0) && (strcmp(DeviceIDStr, "0")!=0)) THROW_EXCEPTION("Error understanding device number.", 1)
    }
#else
    // If no BOOST, we include a file with the parameters and compile it.
    int N, DeviceID = -1;
    Real dT1, dT2, Tcrit, StepSize;
    std::string Filename, Prefix = "";
    #include "noboost.inc"
    P.N = N;
    P.Filename = Filename;
    P.dT1 = dT1;
    P.dT2 = dT2;
    P.Tcrit = Tcrit;
    P.ConstantStep = StepSize;
    P.DeviceID = DeviceID;
#endif
    P.Seed = INT_MIN;
    if ((P.Filename[0]=='_') && (P.Filename.rfind("_") > 0)) {
        int Break = P.Filename.rfind("_") + 1;
        string SeedStr = P.Filename.substr(Break, P.Filename.length() - Break);
        int SeedInt;
        istringstream iss(SeedStr);
        iss >> ws >> P.Seed >> ws;
        if(!iss.eof()) THROW_EXCEPTION("Could not understand random seed (in Filename).", 1)
        P.Filename = P.Filename.substr(0, Break);
    }
    if (P.Seed == INT_MIN) P.Seed = (int)time(NULL);
    *Params = P;
}


// make safety: if single precision, fail to compile
#ifdef ETICS_HDF5
#include "H5Cpp.h"

int H5IterNum;
int *H5IterSnapNumbers;
string *H5IterGroupNames;

herr_t file_info(hid_t loc_id, const char *name, void *opdata) {
    int res = sscanf(name, "Step#%d", H5IterSnapNumbers + H5IterNum);
    if (res != 1) {
        cerr << "Problem understanding group \"" << name << "\" in HDF5 file" << endl;
        exit(1);
    }
    H5IterGroupNames[H5IterNum] = name;
    H5IterNum++;
    return 0;
}


void ReadICsHDF5(string Filename, int N, Particle **P_h, int *FileSnapshotNum, Real *FileTime) {
    double *Mass, *X, *Y, *Z, *VX, *VY, *VZ;
    int *ID;
    Mass  = new double[N];
    X  = new double[N];
    Y  = new double[N];
    Z  = new double[N];
    VX = new double[N];
    VY = new double[N];
    VZ = new double[N];
    ID = new int[N];

    H5::H5File file;
    file = H5::H5File(Filename, H5F_ACC_RDONLY);

    H5::Group group = file.openGroup("/");
    int NumberOfGroups = group.getNumObjs();
    H5IterSnapNumbers = new int[NumberOfGroups];
    H5IterGroupNames = new string[NumberOfGroups];
    H5IterNum = 0;
    int ret = file.iterateElems("/", NULL, file_info, NULL);
    int MaxSnapIndex = 0;
    for (int i=1; i < NumberOfGroups; i++) MaxSnapIndex = (H5IterSnapNumbers[i]>H5IterSnapNumbers[MaxSnapIndex])?i:MaxSnapIndex;
    string GroupName = H5IterGroupNames[MaxSnapIndex];


    group = file.openGroup("/" + GroupName);
    H5::Attribute attribute = group.openAttribute("Time");
    attribute.read(H5::PredType::NATIVE_DOUBLE, FileTime);
    *FileSnapshotNum = H5IterSnapNumbers[MaxSnapIndex];

    H5::DataSet dataset = file.openDataSet("/" + GroupName + "/ID");
    int NfromFile = dataset.getStorageSize()/sizeof(int);
    if (NfromFile < N) {
        cerr << "Was only able to read " << NfromFile << " particles while " << N << " were requested." << endl;
        exit(1);
    }
    dataset.read(ID, H5::PredType::NATIVE_UINT32);

    dataset = file.openDataSet("/" + GroupName + "/Mass");  dataset.read(Mass,  H5::PredType::NATIVE_DOUBLE);
    dataset = file.openDataSet("/" + GroupName + "/X");     dataset.read(X,  H5::PredType::NATIVE_DOUBLE);
    dataset = file.openDataSet("/" + GroupName + "/Y");     dataset.read(Y,  H5::PredType::NATIVE_DOUBLE);
    dataset = file.openDataSet("/" + GroupName + "/Z");     dataset.read(Z,  H5::PredType::NATIVE_DOUBLE);
    dataset = file.openDataSet("/" + GroupName + "/VX");    dataset.read(VX, H5::PredType::NATIVE_DOUBLE);
    dataset = file.openDataSet("/" + GroupName + "/VY");    dataset.read(VY, H5::PredType::NATIVE_DOUBLE);
    dataset = file.openDataSet("/" + GroupName + "/VZ");    dataset.read(VZ, H5::PredType::NATIVE_DOUBLE);

    dataset.close();
    group.close();
    file.close();

    *P_h = new Particle[N];
    for (int i = 0; i < N; i++) {
        Particle p;
        p.ID = ID[i];
        p.m = (Real)Mass[i];
        p.pos = vec3((Real)X[i], (Real)Y[i], (Real)Z[i]);    // This is to ensure proper casting in the case of single prec.
        p.vel = vec3((Real)VX[i], (Real)VY[i], (Real)VZ[i]);
        p.Status = 0;
        p.CalculateR2();
        (*P_h)[i] = p;
    }
}

bool H5FirstSnapshot = true;

void WriteSnapshotHDF5(string Prefix, int SnapNumber, Particle *P_h, int N, Real T) {
    string Filename = Prefix + ".h5part";

    H5::H5File file;

    if (H5FirstSnapshot && std::ifstream(Filename.c_str())) {
        cerr << "File \"" << Filename << "\" already exist!" << endl;
        exit(1);
    }
    H5FirstSnapshot = false;

    if (std::ifstream(Filename.c_str())) {
        file = H5::H5File(Filename, H5F_ACC_RDWR);
    }
    else file = H5::H5File(Filename, H5F_ACC_TRUNC);

    char GroupName[64];
    sprintf(GroupName, "/Step#%d", SnapNumber);

    H5::Group group;
    try {
        H5::Exception::dontPrint();
        group = H5::Group(file.createGroup(GroupName));
    }
    catch (  H5::FileIException error ) {
        cerr << "Couldn't create group \"" << GroupName << "\" in HDF5 file, maybe it already exists?"  << endl;
        exit(1);
    }

    double *Mass, *X, *Y, *Z, *VX, *VY, *VZ, *AX, *AY, *AZ;
    Mass  = new double[N];
    X  = new double[N];
    Y  = new double[N];
    Z  = new double[N];
    VX = new double[N];
    VY = new double[N];
    VZ = new double[N];
    AX = new double[N];
    AY = new double[N];
    AZ = new double[N];
    int *ID;
    ID = new int[N];

    for (int i = 0; i < N; i++) {
        Particle p = P_h[i];
        ID[i] = p.ID;
        Mass[i]  = p.m;
        X[i]  = p.pos.x;
        Y[i]  = p.pos.y;
        Z[i]  = p.pos.z;
        VX[i] = p.vel.x;
        VY[i] = p.vel.y;
        VZ[i] = p.vel.z;
        AX[i] = p.acc.x;
        AY[i] = p.acc.y;
        AZ[i] = p.acc.z;
    }

    H5::DataSpace attr_dataspace = H5::DataSpace(H5S_SCALAR);

    H5::Attribute attribute = group.createAttribute("Time", H5::PredType::NATIVE_DOUBLE, attr_dataspace);
    attribute.write( H5::PredType::NATIVE_DOUBLE, &T);

    attribute = group.createAttribute("TotalN", H5::PredType::NATIVE_UINT32, attr_dataspace);
    attribute.write( H5::PredType::NATIVE_UINT32, &N);

    hsize_t dimsxxx = N;
    H5::DataSpace dataspacexxx(1, &dimsxxx);
    H5::DataSet dataset3;
    dataset3 = H5::DataSet(group.createDataSet("ID", H5::PredType::NATIVE_UINT32, dataspacexxx));
    dataset3.write(ID, H5::PredType::NATIVE_UINT32);
    dataset3 = H5::DataSet(group.createDataSet("Mass",  H5::PredType::NATIVE_DOUBLE, dataspacexxx));
    dataset3.write(Mass, H5::PredType::NATIVE_DOUBLE);
    dataset3 = H5::DataSet(group.createDataSet("X",  H5::PredType::NATIVE_DOUBLE, dataspacexxx));
    dataset3.write(X, H5::PredType::NATIVE_DOUBLE);
    dataset3 = H5::DataSet(group.createDataSet("Y",  H5::PredType::NATIVE_DOUBLE, dataspacexxx));
    dataset3.write(Y, H5::PredType::NATIVE_DOUBLE);
    dataset3 = H5::DataSet(group.createDataSet("Z",  H5::PredType::NATIVE_DOUBLE, dataspacexxx));
    dataset3.write(Z, H5::PredType::NATIVE_DOUBLE);
    dataset3 = H5::DataSet(group.createDataSet("VX", H5::PredType::NATIVE_DOUBLE, dataspacexxx));
    dataset3.write(VX, H5::PredType::NATIVE_DOUBLE);
    dataset3 = H5::DataSet(group.createDataSet("VY", H5::PredType::NATIVE_DOUBLE, dataspacexxx));
    dataset3.write(VY, H5::PredType::NATIVE_DOUBLE);
    dataset3 = H5::DataSet(group.createDataSet("VZ", H5::PredType::NATIVE_DOUBLE, dataspacexxx));
    dataset3.write(VZ, H5::PredType::NATIVE_DOUBLE);
    dataset3 = H5::DataSet(group.createDataSet("AX", H5::PredType::NATIVE_DOUBLE, dataspacexxx));
    dataset3.write(AX, H5::PredType::NATIVE_DOUBLE);
    dataset3 = H5::DataSet(group.createDataSet("AY", H5::PredType::NATIVE_DOUBLE, dataspacexxx));
    dataset3.write(AY, H5::PredType::NATIVE_DOUBLE);
    dataset3 = H5::DataSet(group.createDataSet("AZ", H5::PredType::NATIVE_DOUBLE, dataspacexxx));
    dataset3.write(AZ, H5::PredType::NATIVE_DOUBLE);

    dataset3.close();
    dataspacexxx.close();

    attr_dataspace.close();
    attribute.close();
    group.close();
    file.close();

    delete [] Mass;
    delete [] X;
    delete [] Y;
    delete [] Z;
    delete [] VX;
    delete [] VY;
    delete [] VZ;
    delete [] AX;
    delete [] AY;
    delete [] AZ;
    delete [] ID;
}

#endif
