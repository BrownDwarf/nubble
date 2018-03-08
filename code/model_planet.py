from Starfish.model import *
from Starfish.model import ThetaParam as ThetaParamBase


class ThetaParamPlanet(ThetaParamBase):
    '''
    Inherit the base model and extend it to include an extra parameter: scattering

    :param grid: parameters corresponding to the dimensions of the grid.
    :type grid: 1D np.array
    '''
    def __init__(self, grid, grid2, vz=0.0, vsini=0.0, logOmega=0.0, Av=0.0, teff2=6000.0, logOmega2=-0.5):
        self.grid = grid
        self.grid2 = grid2
        self.vz = vz
        self.vsini = vsini
        self.logOmega = logOmega #log10Omega
        self.Av = Av
        self.teff2 = teff2
        self.logOmega2 = logOmega2


    def __repr__(self):
        return "grid:{} grid2:{} vz:{} vsini:{} logOmega:{} Av:{} teff2:{} ff:{}".format(self.grid, self.grid2, self.vz, self.vsini, self.logOmega, self.Av, self.teff2, self.logOmega2)
