import numpy as np


def dust_temperature(path, filename=None, ext=None, binary=False):

    if (filename == None):
        if (ext == None):
            if binary:
                filename = path + "dust_temperature.bdat"
            else:
                filename = path + "dust_temperature.dat"
        else:
            if binary:
                filename = path + "dust_temperature_"+str(ext)+".bdat"
            else:
                filename = path + "dust_temperature_"+str(ext)+".dat"

    if binary:
        f = open(filename, "rb")
        data = np.fromfile(filename)
    else:
        f = open(filename,"r")

    if binary:
        int.from_bytes(f.read(8), byteorder="little")
        int.from_bytes(f.read(8), byteorder="little")
        ncells = int.from_bytes(f.read(8), byteorder="little")
        nspecies = int.from_bytes(f.read(8), byteorder="little")
    else:
        f.readline()
        ncells = int(f.readline())
        nspecies = int(f.readline())

    temperature = []
    index = 3
    for i in range(nspecies):
        temp = np.empty((ncells,))

        for j in range(ncells):
            if binary:
                temp[j] = data[index+j]
            else:
                temp[j] = float(f.readline())

        temperature.append(temp)

    f.close()

    return temperature


def localfield(path, filename=None):

    # if (filename == None):
    #     filename = path + "mean_intensity.out"


    # f = open(filename,"r")

    # f.readline()
    # ncells = int(f.readline())
    # nfreq = int(f.readline())
    # f.readline()

    # localfield = []
    # for i in range(nfreq*ncells):
    #     localfield.append(float(f.readline()))


    # # for i in range(nfreq):
    # #     field = np.empty((ncells,))

    # #     for j in range(ncells):
    # #         field[j] = float(f.readline())

    # #     localfield.append(field)

    # f.close()

    # return localfield


    if (filename == None):
        filename = path + "mean_intensity.out"
  
    localfield = np.loadtxt(filename, skiprows=4)

    return localfield