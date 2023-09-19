import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from autodeer.dataset import Dataset
from autodeer.utils import sop
from autodeer.classes import Parameter
from scipy import signal
from scipy.linalg import eig
from scipy.sparse import bsr_array
import deerlab as dl


def create_Nmodel(mwFreq):


    def model_func(B, Boffset, gy,gz,axy,az,GB):
        B = B.astype(np.float64)
        gx  =-0.0025 * az + 2.0175
        system = SpinSystem([1/2],[1],[gx, gy, gz], [axy*28.0328,axy*28.0328,az*28.0328])
        system.gn = np.array([0.4038])

        _,y =build_spectrum(system, mwFreq, B,Guass_broadening=GB);
        y_new = np.interp(B,B+Boffset,y)
        return y_new
    

    mymodel = dl.Model(model_func,constants='B')
    # mymodel.mwFreq
    
    # mymodel.gx.par0 = 2.007
    # mymodel.gx.lb = mymodel.gx.par0 - 5e-3
    # mymodel.gx.ub = mymodel.gx.par0 + 5e-3

    mymodel.Boffset.par0 = 0.7
    mymodel.Boffset.lb=-2
    mymodel.Boffset.ub=2
    mymodel.Boffset.unit = 'mT'

    mymodel.gy.par0 = 2.006
    mymodel.gy.lb = mymodel.gy.par0 - 5e-3
    mymodel.gy.ub = mymodel.gy.par0 + 5e-3
    mymodel.gy.freeze(2.0061)

    mymodel.gz.par0 = 2.003
    mymodel.gz.lb = mymodel.gz.par0 - 5e-3
    mymodel.gz.ub = mymodel.gz.par0 + 5e-3
    mymodel.gz.freeze(2.0021)

    # mymodel.axy.par0 = 15
    # mymodel.axy.lb = mymodel.ax.par0 - 10
    # mymodel.axy.ub = mymodel.ax.par0 + 10
    # mymodel.axy.freeze(13.7)
    mymodel.axy.par0 = 0.488
    mymodel.axy.lb = mymodel.axy.par0 - 0.2
    mymodel.axy.ub = mymodel.axy.par0 + 0.2
    mymodel.axy.freeze(0.488)
    mymodel.axy.unit = 'mT'


    # mymodel.az.par0 = 100
    # mymodel.az.lb = mymodel.az.par0 - 10
    # mymodel.az.ub = mymodel.az.par0 + 10
    mymodel.az.par0 = 3.66
    mymodel.az.lb = mymodel.az.par0 - 0.5
    mymodel.az.ub = mymodel.az.par0 + 0.5
    mymodel.az.unit = 'mT'

    mymodel.GB.par0=0.45
    mymodel.GB.lb = 0.15
    mymodel.GB.ub = 0.65

    mymodel.addlinear('scale',lb=0)

    return mymodel
        
class FieldSweepAnalysis():

    def __init__(self, dataset: Dataset) -> None:
        """Analysis and calculation of FieldSweep Experiment. 

        Parameters
        ----------
        dataset : Dataset
            _description_
        """
        self.axis = dataset.axes[0]
        self.data = dataset.data
        self.dataset = dataset
        if hasattr(self.dataset,"LO"):
            self.LO = self.dataset.LO
        pass

    def find_max(self) -> float:
        """Calculates the maximum field

        Returns
        -------
        float
            Max field
        """
        max_id = np.argmax(np.abs(self.data))
        self.max_field = self.axis[max_id]

        return self.max_field

    def calc_gyro(self, LO: float=None) -> float:
        """Calculates the gyromagnetic ratio for a given frequency

        Parameters
        ----------
        det_frq : float
            The detection frequency for the field sweep.

        Returns
        -------
        float
            The gyromagnetic ratio in G/GHz.
        """

        if not hasattr(self, "max_field"):
            self.find_max()

        if LO is None:
            if hasattr(self,"LO"):
                LO = self.LO.value
            else:
                raise ValueError("A LO frequency must eithe be in the dataset or specified as an argument")
            
        self.LO = LO
        self.gyro = LO/self.max_field
        hf_x = LO - self.gyro*self.axis
        self.fs_x = LO + hf_x
        self.fs_x = LO - self.gyro*self.axis
        return self.gyro
    
    def fit(self, spintype='N', **kwargs):

        if spintype != 'N':
            raise ValueError("Currently the fit function only supports Nitroxide spins")

        if isinstance(self.LO,Parameter):
            mymodel  = create_Nmodel(self.LO.value*1e3)

        else:
            mymodel  = create_Nmodel(self.LO*1e3)
        B = np.linspace(self.axis.min(), self.axis.max(), self.data.shape[0])*0.1
        result = dl.fit(mymodel,self.data.real,B,verbose=2,reg=False,  **kwargs)
        self.results = result
        self.model = mymodel
        return result

    def plot(self, norm: bool = True, axis: str = "time") -> Figure:
        """Generate a field sweep plot

        Parameters
        ----------
        norm : bool, optional
            Nomarlisation of the plot to a maximum of 1, by default True
        axis : str, optional
            plot field sweep on either the "time" axis or "freq" axis

        Returns
        -------
        Matplotlib.Figure
            matplotlib figure
        """

        if norm is True:
            data = self.data
            data /= np.max(np.abs(data))

        if axis.lower() == 'time':
            fig, ax = plt.subplots()
            ax.plot(self.axis, np.abs(data), label='abs')
            ax.plot(self.axis, np.real(data), label='real')
            ax.plot(self.axis, np.imag(data), label='imag')
            ax.legend()
            ax.set_xlabel('Field G')
            ax.set_ylabel('Normalised Amplitude')
            if hasattr(self, "max_field"):
                min_value = np.min(np.hstack([np.real(data), np.imag(data)]))
                max_value = np.min(np.hstack([np.real(data), np.imag(data)]))
                ax.vlines(
                    self.max_field, min_value, max_value, label="Maximum")
        elif axis.lower() == 'freq':
            if not hasattr(self, "fs_x"):
                raise RuntimeError("Please run fieldsweep.calc_gyro() first")
            fig, ax = plt.subplots()
            ax.plot(self.fs_x, np.abs(data), label='abs')
            ax.set_xlabel('Frequency GHz')
            ax.set_ylabel('Normalised Amplitude')

        return fig

class SpinSystem:

    def __init__(self,espins, nspin, g,A) -> None:
        self.g = np.array(g)
        self.A = np.array(A)
        self.I = np.array(nspin)
        self.S = np.array(espins)
        self.nElectrons = len(espins)
        self.nNuclei = len(nspin)
        self.Spins = np.concatenate([espins, nspin])

        self.nStates = np.prod(2*self.Spins +1)


        # Defaults
        self.gnscale = 1
        
def erot(*args):
    if len(args) == 0:
        print("erot(varargin)")
        return

    option = ""
    if len(args) == 1 or len(args) == 2:
        angles = np.asarray(args[0])
        if angles.size != 3:
            raise ValueError("Three angles (either separately or in a 3-element array) expected.\n"
                             "You gave an {}x{} array with {} elements.".format(angles.shape[0], angles.shape[1],
                                                                                   angles.size))
        gamma = angles[2]
        beta = angles[1]
        alpha = angles[0]
        if len(args) == 2:
            option = args[1]
    elif len(args) == 3 or len(args) == 4:
        alpha = args[0]
        beta = args[1]
        gamma = args[2]
        if len(args) == 4:
            option = args[3]
    else:
        raise ValueError("Wrong number of input arguments!")

    if not isinstance(option, str):
        raise ValueError("Last argument must be a string, either 'rows' or 'cols'.")

    if option == "":
        return_rows = False
    elif option == "rows":
        return_rows = True
    elif option == "cols":
        return_rows = False
    else:
        raise ValueError("Last argument must be a string, either 'rows' or 'cols'.")

    if (len(args) == 3 or len(args) == 4) and len(args) != 3:
        raise ValueError("3 outputs required if you specify 'rows' or 'cols'.")

    if (len(args) == 3 or len(args) == 4) and len(args) != 3:
        raise ValueError("Please specify whether erot() should return the 3 columns or the 3 rows of the matrix.")

    # Check angles
    if np.isnan(alpha) or np.isnan(beta) or np.isnan(gamma):
        raise ValueError("At least one of the angles is NaN. Angles must be numbers.")

    # Precalculate trigonometric functions of angles
    sa = np.sin(alpha)
    ca = np.cos(alpha)
    sb = np.sin(beta)
    cb = np.cos(beta)
    sg = np.sin(gamma)
    cg = np.cos(gamma)

    # Compute passive rotation matrix
    R = np.array([[cg*cb*ca - sg*sa, cg*cb*sa + sg*ca, -cg*sb],
                  [-sg*cb*ca - cg*sa, -sg*cb*sa + cg*ca, sg*sb],
                  [sb*ca, sb*sa, cb]])

    if len(args) == 3:
        if return_rows:
            return R[0, :], R[1, :], R[2, :]
        else:
            return R[:, 0], R[:, 1], R[:, 2]
    else:
        return R



def eyekron(M):
    size = np.shape(M)[0]
    return np.kron(np.identity(size),M)

def kroneye(M):
    size = np.shape(M)[0]
    return np.kron(M,np.identity(size))

def ham(SpinSystem, elspins=None, nucspins=None):
    # Only using the hyperfine part at the moment as only section present in Nitroxides
    SpinVec = SpinSystem.Spins
    nStates = int(SpinSystem.nStates)
    nElectrons = SpinSystem.nElectrons
    nNuclei = SpinSystem.nNuclei

    Hhf = bsr_array((nStates,nStates),dtype=np.complex128);
    
    if nNuclei == 0: # If there are no Nuclei there this no hyperfine component
        return Hhf
    
    if elspins is None:
        elspins = np.arange(0,nElectrons)

    if nucspins is None:
        nucspins = np.arange(0,nNuclei)

    AMatrix = np.atleast_2d(SpinSystem.A)
    fullAMatrix =  np.size(AMatrix,axis=0) > nNuclei

    # Generate Hamiltonian for hyperfine interaction
    for eSp in elspins:
        eidx = np.arange((eSp - 1) * 3, eSp * 3)
        for nsp in nucspins:

            if SpinSystem.I[nsp] == 0:
                continue

            if fullAMatrix:
                A = AMatrix[np.arange((nsp - 1) * 3, nsp * 3),eidx]
            else:
                A = np.diag(AMatrix[nsp, eidx])


            # TODO: Transform matrix into molecular frame representation 

            for c1,s1 in enumerate(['x','y','z']):
                for c2, s2 in enumerate(['x','y','z']):
                    comps = ['e'] * len(SpinVec)
                    comps[eSp] = s1
                    comps[nElectrons+nsp] = s2
                    comps = ''.join(comps)
                    Hhf += A[c1,c2]*sop(SpinVec,comps)
    
    return (Hhf + Hhf.conj())/2


def ham_ez(SpinSystem, B=None, espins=None):
    bmagn = 9.274010078300000e-24
    planck = 6.626070150000000e-34
    spins = SpinSystem.Spins;
    nElectrons = SpinSystem.nElectrons;
    nStates = int(SpinSystem.nStates);


    if espins is None:
        espins = np.arange(0,nElectrons)

    muxM = bsr_array((nStates,nStates),dtype=np.complex128);
    muyM = bsr_array((nStates,nStates),dtype=np.complex128);
    muzM = bsr_array((nStates,nStates),dtype=np.complex128);
    pre = -bmagn/planck*SpinSystem.g # Hz/T
    pre = pre/1e9 #GHz/T = MHz/mT
    g= np.diag(pre)
    for i in espins:


        for k in range(3):
            comps = ['e'] * len(spins)
            comps[i] = ['x','y','z'][k]
            comps = ''.join(comps)
            Sk = sop(spins,comps)
            muxM += g[k,0]*Sk
            muyM += g[k,1]*Sk
            muzM += g[k,2]*Sk
        

    if B is None:
        return muxM, muyM, muzM
    else:
        return -(muxM*B[0] + muyM*B[1] + muzM*B[2])
    

def ham_nz(SpinSystem, B=None, nspins=None):
    bmagn = 9.274010078300000e-24
    nmagn = 5.050783746100000e-27

    planck = 6.626070150000000e-34
    spins = SpinSystem.Spins;
    nElectrons = SpinSystem.nNuclei;
    nStates = int(SpinSystem.nStates);


    if nspins is None:
        nspins = np.arange(0,nElectrons)

    muxM = bsr_array((nStates,nStates),dtype=np.complex128);
    muyM = bsr_array((nStates,nStates),dtype=np.complex128);
    muzM = bsr_array((nStates,nStates),dtype=np.complex128);
    pre = +nmagn/planck * SpinSystem.gn * SpinSystem.gnscale # Hz/T
    pre = pre/1e9 #GHz/T = MHz/mT
    g= np.diag(pre)


    for i in nspins:

        #TODO: add sigma

        sigma = np.identity(3)
        
        for k in range(3):
            comps = ['e'] * len(spins)
            comps[nElectrons+i] = ['x','y','z'][k]
            comps = ''.join(comps)
            Sk = sop(spins,comps)
            muxM += pre[i]*Sk*sigma[0,k]
            muyM += pre[i]*Sk*sigma[1,k]
            muzM += pre[i]*Sk*sigma[2,k]
        

    if B is None:
        return muxM, muyM, muzM
    else:
        return -(muxM*B[0] + muyM*B[1] + muzM*B[2])

def sphgrid(nOctants,maxPhi,GridSize,closedPhi):

    dtheta = (np.pi/2)/(GridSize - 1)

    if nOctants > 0:

        if nOctants == 8:
            nOct = 4
        else:
            nOct = nOctants

        nOrientations = int(GridSize + nOct * GridSize * (GridSize - 1)/2)

        phi = np.zeros(nOrientations)
        theta = np.zeros(nOrientations)
        Weights = np.zeros(nOrientations)

        sindth2 = np.sin(dtheta/2)
        w1 = 0.5

        if not closedPhi:
            w1 = w1 + 1/2

        phi[0] = 0.
        theta[0] = 0.
        Weights[0] = maxPhi * (1 -np.cos(dtheta/2))

        # Non Equatorial Points
        Start = 1
        for iSlice in np.arange(2, GridSize):
            nPhi = nOct * (iSlice - 1) + 1
            dphi = maxPhi / (nPhi - 1)
            idx = Start + np.arange(0,stop = nPhi)
            Weights[idx] = 2 * np.sin((iSlice - 1) * dtheta) * sindth2 * dphi * np.concatenate([[w1], np.ones(nPhi - 2), [0.5]])
            phi[idx] = np.linspace(0,maxPhi,nPhi)
            theta[idx] = (iSlice - 1) * dtheta
            Start = Start + nPhi



        # Equatorial Points
        nPhi = nOct * (GridSize - 1) + 1;
        dPhi = maxPhi / (nPhi - 1);
        idx = Start + np.arange(0,stop = nPhi)
        phi[idx] = np.linspace(0,maxPhi,nPhi)
        theta[idx] = np.pi / 2
        Weights[idx] = sindth2 * dPhi * np.concatenate([[w1], np.ones(nPhi - 2), [0.5]])

        # Border Removal

        if not closedPhi:
            rmv = np.cumsum(nOct * (np.arange(1,stop = GridSize - 1)) + 1) + 1
            phi =np.delete(phi,rmv)
            theta = np.delete(theta,rmv)
            Weights = np.delete(Weights,rmv)

        if nOctants == 8: # C1 Symmetry
            idx = np.size(theta) - np.arange(nPhi+1,stop = -1,step = 1) 
            phi = np.concatenate([phi, phi[idx]])
            theta = np.concatenate([theta , np.pi - theta[idx]])
            Weights[idx] = Weights[idx] / 2
            Weights = np.concatenate([Weights , Weights[idx]])

        Weights = 2 * (2 * np.pi / maxPhi) * Weights # Sum = 4* pi

    elif nOctants == 0: # Dinfh symmetry (quarter of meridian in xz plane)
        phi = np.zeros(1,GridSize)
        theta = np.linspace(0,phi/2,GridSize)
        Weights = -2 * (2 * np.pi) * np.diff(np.cos([0, np.arange(dtheta/2,stop = dtheta,step = np.pi/2), np.pi/2]))

    elif nOctants == -1: # O3 Symmetry

        phi = 0;
        theta = 0;
        Weights =0;

    else:
        raise ValueError("Unsupported value ",nOctants," for nOctants")


    return phi,theta,Weights


def resfields(system, Orientations, mwFreq, computeIntensities = True,
              RejectionRatio = 1e-8, Range = (0,1e8),Threshold = 0, computeFreq2Field = True):


    # Generate orientations
    nOrientations = Orientations.shape[0]
    averageOverChi = True

    H0 = ham(system)
    [muxe,muye,muze] = ham_ez(system)
    [muxn,muyn,muzn] = ham_nz(system)
    [mux,muy,muz] = [muxe + muxn,muye + muyn,muze + muzn]


    A = eyekron(H0) - kroneye(H0.conj()) + mwFreq*np.eye(H0.shape[0]**2);
    E = np.diag(eig(A, right=False))

    if computeIntensities:

        mux_vec = mux.flatten()
        muy_vec = muy.flatten()
        muz_vec = muz.flatten()

    EigenFields = []
    Intensities = []

    for iOri,Ori in enumerate(Orientations):
        
        [xLab,yLab,zLab] = erot(Ori,'rows')
        muzL = zLab[0]*mux + zLab[1]*muy + zLab[2]*muz

        B = - kroneye(muzL.conj()) + eyekron(muzL)
        
        if computeIntensities:
            [Fields,Vecs] = eig(A,B);
            idx  = np.argsort(Fields);
            Fields = Fields[idx]
            Vecs = Vecs[:, idx]

        mask = np.abs(Fields.imag) < (np.abs(Fields.real)* RejectionRatio) 
        mask &= np.greater(Fields, 0)
        mask &= np.isfinite(Fields)
        mask &= np.greater(Fields, Range[0])
        mask &= np.less(Fields, Range[1])

        if np.equal(mask, False).all():
            EigenFields.append([])
            Intensities.append([])
        else:
            EigenFields.append(Fields[mask].real)
            Vecs = Vecs[:,mask]

        # Normalize eigenvectors to unity
        Norms = np.sqrt(np.sum(np.abs(Vecs)**2,axis=0))
        Vecs /= Norms[None,:]

        # Assuming never parallel mode
        muxL_vec = xLab[0]*mux_vec + xLab[1]*muy_vec + xLab[2]*muz_vec
        if averageOverChi:
            muyL_vec = yLab[0]*mux_vec + yLab[1]*muy_vec + yLab[2]*muz_vec
            TransitionRate = (np.abs((muxL_vec[:,None]*Vecs).sum(axis=0))**2 + np.abs((muyL_vec[:,None]*Vecs).sum(axis=0)**2))/2
        else:
            TransitionRate = np.abs((muxL_vec[:,None]*Vecs).sum(axis=0))**2


        Polarization = 1;
        Polarization = Polarization/np.prod(2*system.I+1);

        if computeFreq2Field:

            n = H0.shape[0]
            Vecs = np.reshape(Vecs,(n,n, int(Vecs.size/n**2)),order='F')
            dBdE = np.zeros(Vecs.shape[2])
            for iVec in range(Vecs.shape[2]):
                V = Vecs[:,:,iVec]
                dBdE[iVec] = 1/np.abs(np.trace(-muzL@(V@V.conj().T - V.conj().T@V)))
        else:
            dBdE = np.ones(TransitionRate.shape)
            
        # Combine factors
        Intensities.append(Polarization * np.real(TransitionRate*dBdE).T)
        mask = Intensities[iOri] >= Threshold*Intensities[iOri].max();
        EigenFields[iOri] = EigenFields[iOri][mask];
        Intensities[iOri] = Intensities[iOri][mask];


    return EigenFields, Intensities

def build_spectrum(system, mwFreq, Range, knots=19,npoints = 1000, Guass_broadening=0.25):
    """Build a field sweep spectrum

    Parameters
    ----------
    system : SpinSystem
        The spin system it must include: I & S spins, g, A, gn
    mwFreq : float
        The microwave frequency in MHz
    Range : float
        The field range in mT
    knots : int, optional
        The number of knots of orientation averaging, by default 19
    npoints : int, optional
        The number of points in the spectrum, by default 1000

    Returns
    -------
    xAxis: np.ndarray
        The xAxis in mT
    y: np.ndarray
        The spectrum intensities normalised to 1  
    """

    phi,theta,Weights = sphgrid(1,np.pi/2,knots,True)
    Orientations = np.vstack([phi,theta, np.zeros(phi.shape)]).T
    nOrientations = Orientations.shape[0]
    EigenFields, Intensities = resfields(system,Orientations,mwFreq)

    nReson = 0
    for k in EigenFields:
        nReson += k.size

    nSites = 1

    if isinstance(Range, np.ndarray):
        xAxis = Range
        xmin = Range.min()
        xmax = Range.max()
        npoints = Range.shape[0]
        prefactor = (npoints - 1)/(Range.max()-Range.min())
    else:
        xAxis = np.linspace(*Range,npoints)
        xmin = Range[0]
        xmax = Range[1]

        prefactor = (npoints - 1)/(Range[1]-Range[0])
    dx = xAxis[1]-xAxis[0]
    spec = np.zeros(npoints)
    

    for iOri in range(nOrientations):
        thisP = EigenFields[iOri]
        Amplitudes = Intensities[iOri]
        idxPos = np.around(1+prefactor*(thisP-xmin))
        outofRange = np.less(idxPos,1) | np.greater(idxPos, npoints)
        spec[idxPos[~outofRange].astype(int)] += Amplitudes[~outofRange] * Weights[iOri]

    # Convolution broadening


    win = signal.windows.gaussian(npoints, Guass_broadening/dx)
    filtered = signal.convolve(spec, win, mode='same')
    filtered /= filtered.max()

    return xAxis, filtered