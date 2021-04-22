#!/usr/bin/env python
# This script converts phiGRAPE snapshot files to HDF5.
# Author: Yohai Meiron

import argparse, os, sys, numpy, re

if __name__ == "__main__":
    # Set commandline parameters
    parser = argparse.ArgumentParser(description='Convert phiGRPAE snapshot files into HDF5 (h5part) format.')
    parser.add_argument('InputFile', metavar='InputFile', type=str, help='first input file (numbers in the file name must represet numerical order)')
    parser.add_argument('--num', type=int, help='number of snapshot to convert (by default, all are converted)')
    parser.add_argument('-o', '--output', default='output', help='output file name (or pattern)')
    parser.add_argument('--particles', type=int, help='number of particles (if not the whole thing)')
    parser.add_argument('--maxfs', type=float, help='maximum size in MB of a single HDF5 output')
    parser.add_argument('--single', action='store_true', help='use single precision instead of double')
    parser.add_argument('--bh', type=int, nargs='?', const=1, default=0, help='number of black holes to grab from the end of the file')

    # Print message and parse input
    print('Welcome to Yohai\'s format converter! (phiGRAPE ---> HDF5)\n')
    args = parser.parse_args()

    # Parameter conflict tests
    if (args.bh > 0) and (args.particles==None):
        parser.print_usage()
        print('%s: error: cannot --bh without --particles' % sys.argv[0])
        sys.exit(2)
    Nbh = args.bh

    # Verify h5py
    try:
        import h5py
    except:
        print('%s: error: booooo, can\'t find h5py.' % sys.argv[0])
        sys.exit(1)

    # Check if input file exists, separate name from path and get file name pattern
    if not os.path.exists(args.InputFile):
        print('File `%s` not found!' % args.InputFile)
        sys.exit(1)
    InputFileName = os.path.basename(args.InputFile)
    InputFilePath = os.path.dirname(os.path.realpath(args.InputFile))
    Pattern = re.sub(r'\d', r'\d', InputFileName)
    Pattern = '^' + Pattern + "$"
    r = re.compile(Pattern)

    # Look for the other files
    ListDir = os.listdir(InputFilePath)
    FileList = []
    for Filename in ListDir:
        m = r.search(Filename)
        if not m is None:
            FileList.append(Filename)
    FileList.sort()
    FileList = FileList[FileList.index(InputFileName):]

    # Verify that enough files were found; shorten list if needed
    if (not args.num==None) and args.num > len(FileList):
        print('Requested to convert %d files, but only %d found!' % (args.num, len(FileList)))
        sys.exit(1)
    if (not args.num==None): num = args.num
    else: num = len(FileList)
    FileList = FileList[:num]

    print('Will convert %d files in the directory `%s`.' % (num, InputFilePath))

    NumperOfStepsPerFile = 2147483647
    OutputFileNum = 2
    StepCounter = 0
    FirstTime = True
    OutputFile = False

    for Filename in FileList:
        if StepCounter > NumperOfStepsPerFile:
            OutputFile.close()
            OutputFile = h5py.File('output(%d).h5part' % OutputFileNum, 'w')
            StepCounter = 0
            OutputFileNum += 1
        FullFilename = os.path.join(InputFilePath, Filename)
        with open(FullFilename, 'r') as File:
            Line = File.readline()
            try:
                StepNumer = int(Line)
            except:
                print('Error to read snapshot number (line 1 in file `%s`)' % Filename)
                sys.exit(1)
            File.readline()
            Line = File.readline()
            try:
                Time = numpy.double(Line)
            except:
                print('Error to read snapshot time (line 3 in file `%s`)' % Filename)
                sys.exit(1)
        Data = numpy.loadtxt(FullFilename, skiprows=3)
        if (not args.particles==None) and (args.particles > Data.shape[0]):
            print('Requested to use %d paticles, but file %s has only %d!' % (args.particles, Filename, Data.shape[0]))
            sys.exit(1)
        if (not args.particles==None): N = args.particles
        else: N = Data.shape[0]
        if FirstTime:
            if Nbh == 0: print('Reading %d particles.' % (N))
            else: print('Reading %d normal particles + %d black hole(s).' % (N, Nbh))
            if args.single: EstimatedDataSize = 32*N/1024.0**2
            else: EstimatedDataSize = 60*N/1024.0**2
            if not args.maxfs==None: NumperOfStepsPerFile = (int)(args.maxfs/EstimatedDataSize)
            FirstTime = False
            if NumperOfStepsPerFile > num: OutputFile = h5py.File('output.h5part', 'w')
            else: OutputFile = h5py.File('output(1).h5part', 'w')

        DataBH = Data[-Nbh-1:-1]
        Data = Data[:N]
        Data = numpy.vstack((Data, DataBH))
        Group = OutputFile.create_group('Step#%d' % (StepNumer))
        Group.attrs['Time'] = Time
        Group.attrs['TotalN'] = N + Nbh
        if args.single: T='f'
        else: T='d'
        DataSet = Group.create_dataset('ID',   (N+Nbh,), dtype='i'); DataSet[...] = Data[:,0]
        DataSet = Group.create_dataset('Mass', (N+Nbh,), dtype=T);   DataSet[...] = Data[:,1]
        DataSet = Group.create_dataset('X',    (N+Nbh,), dtype=T);   DataSet[...] = Data[:,2]
        DataSet = Group.create_dataset('Y',    (N+Nbh,), dtype=T);   DataSet[...] = Data[:,3]
        DataSet = Group.create_dataset('Z',    (N+Nbh,), dtype=T);   DataSet[...] = Data[:,4]
        DataSet = Group.create_dataset('VX',   (N+Nbh,), dtype=T);   DataSet[...] = Data[:,5]
        DataSet = Group.create_dataset('VY',   (N+Nbh,), dtype=T);   DataSet[...] = Data[:,6]
        DataSet = Group.create_dataset('VZ',   (N+Nbh,), dtype=T);   DataSet[...] = Data[:,7]
        OutputFile.flush()
        print('Completed /Step#%d' % (StepNumer))
        StepCounter += 1
    OutputFile.close()
