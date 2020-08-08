import numpy as np
import copy
from Control import *

N_REACT = 5  # D, Siga, nuSigf, Xi, sigr
DIF    = 0
SIGA   = 1
NUSIGF = 2
XI     = 3
SIGR   = 4

class Material:
    def __init__(self, val=None, ng=ng):
        """
        constructor
        """
        self.x  = np.zeros((ng, N_REACT))
        self.sm = np.zeros((ng, ng))
        if not (val is None):
            self.set(val)

    def read_file(self,path,file):
        # print("Reading XS file from",file)
        data=np.loadtxt(path+"/"+file)
        # #self.__init__(data[0])
        ig=int(data[0])
        if(ig != self.ng):
            print(colored("Error:","red"),"read XS file, group number doesn't match!")
        # self.enband=np.zeros(ig+1)
        nusigf=np.zeros(ig)
        siga=np.zeros(ig)
        sign2n=np.zeros(ig)
        sigt = np.zeros(ig)
        d=np.zeros(ig)
        chi=np.zeros(ig)
        sigs=np.zeros((ig,ig))

        ipl=int(data[1])
        pl=ipl
        b1=2
        b2=2+ig+1
        for i in range(b1,b2):
            self.enband[i-b1]=data[i]
        b1=b2
        b2=b1+ig
        for i in range(b1,b2):
            nusigf[i-b1]=data[i]
        b1=b2
        b2=b1+ig
        for i in range(b1,b2):
            siga[i-b1]=data[i]
        b1=b2
        b2=b1+ig
        for i in range(b1,b2):
            chi[i-b1]=data[i]
        b1=b2
        b2=b1+ig
        for i in range(b1,b2):
            sign2n[i-b1]=data[i]
        b1=b2
        b2=b1+ig
        for i in range(b1,b2):
            sigt[i-b1]=data[i]           
        # 1D macroscopic total cross section 
        # weighted by high-order angular moment of neutron flux
        sigt1=np.zeros(ig)
        b1=b2
        b2=b1+ig
        for i in range(b1,b2):
            sigt1[i-b1]=data[i]
        #Average and anisotropic diffusion coefficient
        # d dr dz
        b1=b2
        b2=b1+ig
        for i in range(b1,b2):
            d[i-b1]=data[i]
        b1=b2
        b2=b1+ig
        # for i in range(b1,b2):
        #     dr[i-b1]=data[i]
        b1=b2
        b2=b1+ig
        # for i in range(b1,b2):
        #     dz[i-b1]=data[i]
        if(pl!=0):
            sigs_pl=np.zeros((ig,ig,pl))
            for l in range(pl):
                for ii in range(ig):
                    b1=b2
                    tmp1=int(data[b1])
                    tmp2=int(data[b1+1])
                    b1=b1+2
                    b2=b1+tmp2-tmp1+1
                    for jj in range(b1,b2):
                        sigs_pl[ii,jj-b1+tmp1,l]=data[jj]
            sigs=sigs_pl[:,:,0]
        else:
            sigs_pl=np.zeros((ig,ig))
            for ii in range(ig):
                b1=b2
                tmp1=int(data[b1])
                tmp2=int(data[b1+1])
                b1=b1+2
                b2=b1+tmp2-tmp1+1
                for jj in range(b1,b2):
                    sigs_pl[ii,jj-b1+tmp1]=data[jj]
            sigs=sigs_pl
        for kg in range(ig):
            self.x[kg, DIF]   = d[kg].tolist()
            self.x[kg, SIGA]  = siga[kg].tolist()
            self.x[kg, NUSIGF]= nusigf[kg].tolist()
            self.x[kg, XI]    = chi[kg].tolist()
            for kkg in range(ig):
                self.sm[kg, kkg]=sigs[kg,kkg].tolist()

            self.x[kg, SIGR] = self.x[kg, SIGA]
            for kkg in range(self.ng):
                if kg != kkg:
                    self.x[kg, SIGR] += self.sm[kg,kkg]

    def set(self, val):
        if type(val) == Material:
            self.x = copy.copy(val.x)
            self.sm = copy.copy(val.sm)
        else:
            for kg in range(self.ng()):
                for i in range(len(val[kg])):
                    self.x[kg, i] = val[kg][i]

    def set_d(self, kg, val):
        self.x[kg, DIF] = val

    def set_siga(self, kg, val):
        self.x[kg, SIGA] = val

    def set_nusigf(self, kg, val):
        self.x[kg, NUSIGF] = val

    def set_xi(self, kg, val):
        self.x[kg, XI] = val

    def set_smat(self, mat):
        for kg in range(self.ng()):
            for i in range(self.ng()):
                self.sm[kg, i] = mat[kg][i]

    def calc_sigr(self):
        for kg in range(self.ng()):
            self.x[kg, SIGR] = self.x[kg, SIGA]
            for kkg in range(self.ng()):
                if kg != kkg:
                    self.x[kg, SIGR] += self.sm[kg,kkg]

    def calc_matrix_a(self, k):
        if(cplx):
            k=k.real
        else:
            k=k
        for kg in range(self.ng()):
            self.x[kg, SIGR] = self.x[kg, SIGA]
            for kkg in range(self.ng()):
                if kg != kkg:
                    self.x[kg, SIGR] += self.sm[kg,kkg]

        a1 = np.zeros((self.ng(), self.ng()))
        a2 = np.zeros((self.ng(), self.ng()))
        for i in range(self.ng()):
            for j in range(self.ng()):
                a1[i,j] = self.x[i,XI]*self.x[j,NUSIGF] / (self.x[i,DIF]*k)
                if i == j:
                    a2[i,j] = - self.x[i,SIGR]/self.x[i,DIF]
                else:
                    a2[i,j] = self.sm[j,i]/self.x[i,DIF]
        return (a1+a2)

    def calc_matrix_a_sp3(self, k):
        if(cplx):
            k=k.real
        else:
            k=k
        if (self.ng()>2):
            print("error")
        a = np.zeros((2,2))
        d1 = self.x[0,DIF]
        st = self.x[0,SIGA]
        st0 = 0.013 
        d3 = 9./(35*st)
        nsf = self.x[0,NUSIGF]
        
        a[0,0]=-st0/d1+1./d1/keff*nsf
        a[0,1]=-2./d1/keff*nsf+2./d1*st0
        a[1,0]=-35./27./d3*(-0.4*st0+0.4/keff*nsf)
        a[1,1]=-35./27./d3*(st+0.8*st0-0.8/keff*nsf)
        
        return a
                
    def ng(self):
        return self.x.shape[0]

    def dif(self, kg):
        return self.x[kg, DIF]

    def siga(self, kg):
        return self.x[kg, SIGA]

    def sigr(self, kg):
        return self.x[kg, SIGR]

    def nusigf(self, kg):
        return self.x[kg, NUSIGF]

    def xi(self, kg):
        return self.x[kg, XI]

    def sigs(self, kg, kkg):
        return self.sm[kg, kkg]



    def show_xs(self):
        teams_list=['nusigf','siga','chi','diff']
        f1 ="{:^9}"
        print( f1.format( "nusigf"),f1.format("siga"),f1.format("chi"),f1.format("diff"))
        for i in range(self.ng()):
            print("g = ",str(i+1),"&","%7.5e"%self.x[i,NUSIGF],"&", "%7.5e"%self.x[i,SIGA],"&", "%7.5e"%self.x[i,XI], "&","%7.5e"%self.x[i,DIF],"\\\\")

    def show_scatter_xs(self):
        for i in range(self.ng()):
            print("g = ",str(i+1),end="")
            for j in range(self.ng()):
                print("&", "%7.5e"%self.sm[i,j],end='')
            print("\\\\")
        # print(self.sm)


    def debug(self):
        print("-" * 9 + " XS " + "-" * 29)
        print("kg\tD\tSiga\tSigr\tNuSigf\tXi")
        for kg in range(self.x.shape[0]):
            print(kg, self.x[kg, DIF], self.x[kg, SIGA], self.x[kg, SIGR], self.x[kg, NUSIGF], self.x[kg, XI])
        print( "smat")
        print( self.sm )
        print("-"*42)
