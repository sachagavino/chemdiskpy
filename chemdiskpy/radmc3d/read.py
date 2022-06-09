import numpy as np


def grid(path, filename=None):

    if (filename == None):
        filename = path + "amr_grid.inp"
 


    f = open(filename,"r")


    f.readline()
    grid_type = int(f.readline())

    if grid_type == 0:
        coordsystem = int(f.readline())
        f.readline()
        incl_x, incl_y, incl_z = [int(x) for x in f.readline().split()]
        nx, ny, nz = [int(x) for x in f.readline().split()]

        grid = []

        x, y, z = [], [], []

        for ix in range(nx+1):
            x.append(float(f.readline()))
        for iy in range(ny+1):
            y.append(float(f.readline()))
        for iz in range(nz+1):
            z.append(float(f.readline()))

    else:
        print('grid type must be regular (1) for the moment.')
        exit(1)

    f.close()

    return nx, ny, nz, x, y, z



def dust_density(path, filename=None, ext=None, binary=False):

    if (filename == None):
        if (ext == None):
            if binary:
                filename = path + "dust_density.binp"
            else:
                filename = path + "dust_density.inp"
        else:
            if binary:
                filename = path + "dust_density_"+str(ext)+".binp"
            else:
                filename = path + "dust_density_"+str(ext)+".inp"

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

    density = []
    index = 3
    #for i in range(nspecies):
    dens = np.empty((nspecies*ncells,))

    for j in range(nspecies*ncells):
        if binary:
            dens[j] = data[index+j]
        else:
            dens[j] = float(f.readline())

    density.append(dens)

    f.close()

    return density

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
    #for i in range(nspecies):
    temp = np.empty((nspecies*ncells,))

    for j in range(nspecies*ncells):
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
    # # for i in range(nfreq*ncells):
    # #     localfield.append(float(f.readline()))

    # for i in range(nfreq):
    #     field = np.empty((ncells,))

    #     for j in range(ncells):
    #         field[j] = float(f.readline())

    #     localfield.append(field)

    # f.close()

    # return localfield


    if (filename == None):
        filename = path + "mean_intensity.out"
  
    localfield = np.loadtxt(filename, skiprows=4)

    return localfield