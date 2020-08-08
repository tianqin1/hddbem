import numpy as np
import copy
from Control import *
from termcolor import colored
import pathlib

N_REACT = 5  # D, Siga, nuSigf, Xi, sigr
DIF    = 0
SIGA   = 1
NUSIGF = 2
XI     = 3
SIGR   = 4

class Material:

    def __init__(self, val):
        """
        constructor
        """
        self.ng=val
        self.x  = np.zeros((self.ng, N_REACT))
        self.sm = np.zeros((self.ng, self.ng))
        self.enband = np.zeros(self.ng+1)

        # if not (val is None):
        #     self.set(val)
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
    def show_xs(self):
        teams_list=['higher','lower','nusigf','siga','chi','diff']
        f1 ="{:^9}"
        print( f1.format("higher"),f1.format("lower"),f1.format( "nusigf"),f1.format("siga"),f1.format("chi"),f1.format("diff"))
        for i in range(self.ng):
            print("%7.5e"%self.enband[i],"%7.5e"%self.enband[i+1],"%7.5e"%self.x[i,NUSIGF], "%7.5e"%self.x[i,SIGA], "%7.5e"%self.x[i,XI], "%7.5e"%self.x[i,DIF])

    def show_scatter_xs(self):
        print(self.sm)

    def set(self, data):
        if type(data) == Material:
            self.x = copy.copy(data.x)
            self.sm = copy.copy(data.sm)
        else:
            for kg in range(self.ng):
                for i in range(len(data[kg])):
                    self.x[kg, i] = data[kg][i]

    def set_d(self, kg, val):
        self.x[kg, DIF] = val

    def set_siga(self, kg, val):
        self.x[kg, SIGA] = val

    def set_nusigf(self, kg, val):
        self.x[kg, NUSIGF] = val

    def set_xi(self, kg, val):
        self.x[kg, XI] = val

    def set_smat(self, mat):
        for kg in range(self.ng):
            for i in range(self.ng):
                self.sm[kg, i] = mat[kg][i]

    def calc_sigr(self):
        for kg in range(self.ng):
            self.x[kg, SIGR] = self.x[kg, SIGA]
            for kkg in range(self.ng):
                if kg != kkg:
                    self.x[kg, SIGR] += self.sm[kg,kkg]

    def calc_matrix_a(self, k):
        if(cplx):
            k=k.real
        else:
            k=k
        self.calc_sigr()
        a1 = np.zeros((self.ng, self.ng))
        a2 = np.zeros((self.ng, self.ng))
        for i in range(self.ng):
            for j in range(self.ng):
                a1[i,j] = self.x[i,XI]*self.x[j,NUSIGF] / (self.x[i,DIF]*k)
                if i == j:
                    a2[i,j] = - self.x[i,SIGR]/self.x[i,DIF]
                else:
                    a2[i,j] = self.sm[j,i]/self.x[i,DIF]
        return (a1+a2)

    def calc_matrix_bc(self,k):
        amat=self.calc_matrix_a(k)
        bmat=np.zeros(self.ng)
        cmat=np.zeros((self.ng,self.ng))
        bmat,cmat=np.linalg.eig(amat)
        for i in range(self.ng-1):
            cmat[i,:]=cmat[i,:]/cmat[-1,:]
        cmat[-1,:]=np.ones(self.ng)
        return bmat,cmat
    def get_up_scatter_ngroup(self):
        for i in range(1,self.ng):
            for j in range(i):
                if(self.sm[i,j]!=0):
                    return i
        return i
    def cal_flx_infinite_homogeneous(self):
        #print(self.ng)
        self.calc_sigr()
        f_s = 1.0
        s = np.zeros(self.ng)
        phi = np.zeros(self.ng)

        group_boundary = self.get_up_scatter_ngroup()

        #self.show_scatter_xs()
        for i in range(group_boundary):
            s[i]=0
            for j in range(i):
                s[i]+=self.sm[j,i]*phi[j]
            phi[i]=(s[i]+self.x[i,XI]*f_s)/self.x[i, SIGR]
        error=0
        max_error=0
        k_inf_old=0
        k_inf=0
        for iter in range(100):
            max_error=0
            for i in range(group_boundary, self.ng):
                s[i]=0
                for j in range(self.ng):
                    if(i!=j):
                        s[i]+=self.sm[j,i]*phi[j]
                    else:
                        s[i]+=0
                phi[i]=(s[i]+self.x[i,XI]*f_s)/self.x[i,SIGR]
            f_s_new=0
            k_inf_old=k_inf
            for i in range(self.ng):
                f_s_new+=self.x[i,NUSIGF]*phi[i]
            k_inf=f_s_new/f_s
            if(k_inf_old >0):
                error=abs(k_inf-k_inf_old)/k_inf_old
            if(error>max_error):
                max_error=error
                if(max_error<1e-6):
                    break
        #norm_phi=phi/sum(phi)
        return k_inf, phi#,norm_phi


    # def ng(self):
    #     return self.x.shape[0]

    def dif(self, kg):
        return self.x[kg, DIF]

    def sig_a(self, kg):
        return self.x[kg, SIGA]

    def sig_r(self, kg):
        return self.x[kg, SIGR]

    def nusigf(self, kg):
        return self.x[kg, NUSIGF]

    def chi(self, kg):
        return self.x[kg, XI]

    def sig_s(self, kg, kkg):
        return self.sm[kg, kkg]


    def debug(self):
        print("-" * 9 + " XS " + "-" * 29)
        print("kg\tD\tSiga\tSigr\tNuSigf\tXi")
        for kg in range(self.x.shape[0]):
            print(kg, self.x[kg, DIF], self.x[kg, SIGA], self.x[kg, SIGR], self.x[kg, NUSIGF], self.x[kg, XI])
        print( "smat")
        print( self.sm )
        print("-"*42)
# mat = Material(ng)
# path="/home/tian/Desktop/hddbem/bench_final/XS/CBGXS"
# mat.read_file(path,"SOL1")
# k,phi=mat.cal_flx_infinite_homogeneous()
# import matplotlib.pyplot as plt
# cbz_inflx=np.loadtxt("./bench_final/cbz_inflx.txt")
# plt.subplot(2,1,1)
# plt.plot(phi,'-sr',label="homogeneous slowing-down equation")
# plt.legend()
# plt.grid(True)
# plt.subplot(2,1,2)
# plt.plot(cbz_inflx[0:107,1],'-sk',label="homogeneous B1 equation from CBZ")
# plt.legend()
# plt.grid(True)
# plt.savefig("b1_equaiton.png",bbox_inches='tight',dpi=600)
# plt.show()
# #print(cbz_inflx[0:107,0])
