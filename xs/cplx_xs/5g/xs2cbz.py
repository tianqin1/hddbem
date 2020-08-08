import numpy as np

file1 = 'cplx_1d_xs.txt'
file2 = 'cplx_2d_xs.txt'
xs1d=np.loadtxt(file1)
ng=int(len(xs1d)/12)
readxs1d = np.reshape(xs1d,(2,ng,6))
xs1d=np.zeros((2,ng,4))
xs1d[:,:,0]=readxs1d[:,:,3]
xs1d[:,:,1]=readxs1d[:,:,2]
xs1d[:,:,2]=readxs1d[:,:,1]
xs1d[:,:,3]=readxs1d[:,:,-1]
xs2d = np.reshape(np.loadtxt(file2),(2,ng,ng))
for nr in range(2):
    print("//### med[",nr,"]")
    for nk in range(4):
        for g in range(ng):
            print("xs[",g,"]=",xs1d[nr,g,nk],";")
        if(nk==0):
            print("med[",nr,"].PutData(d,xs);")
        elif(nk==1):
            print("med[",nr,"].PutData(siga,xs);")
        elif(nk==2):
            print("med[",nr,"].PutData(nusigf,xs);")
        elif(nk==3):
            print("med[",nr,"].PutData(chi,xs);")


    for ig in range(ng):
        for igg in range(ng):
            print("xx[",ng*ig+igg,"]=",xs2d[nr,ig,igg],";")
    print("med[",nr,"].PutDataSigs(0,xx);")
    for ig in range(ng):
        print("xs[",ig,"]=",xs1d[nr,ig,1]+sum(xs2d[nr,ig,:]),";")
    print("med[",nr,"].PutData(sigt,xs);")
