import os
import numpy as np
import matplotlib.pyplot as plt
from Material import Material
import Green
from Control import *
import Control as ctrol

if(CBGXS):
    scf = Material(ng)
    rbl = Material(ng)
    rdr = Material(ng)

    fuel=Material(ng)
    water=Material(ng)
    path = "../xs/_xs_high_up_scattering/"
    path = path+str(ng)+"g"
    # scf.read_file(path, "r1.txt")
    # rbl.read_file(path, "r2.txt")
    # rdr.read_file(path, "r3.txt")
    fuel.read_file(path,"fuel.txt")
    water.read_file(path,"water.txt")

#med = np.array([fuel, water, fuel])
med = np.array([fuel, water,water ])

def bem(xinp):

    if(cplx):
        k_ini = xinp[0].real
    else:
        k_ini = xinp[0]
    ini_mode = xinp[1:]
    ini_mode = np.reshape(ini_mode, (nmed-1, ng))
    # region 1
    modex0 = np.zeros(ng, dtype=dataType)
    modex1_left = ini_mode[0]
    dmodex1_left = np.zeros(ng, dtype=dataType)
    # region 2
    modex1_right = np.zeros((nmed-2,ng), dtype=dataType)
    dmodex1_right = np.zeros((nmed-2,ng), dtype=dataType)
    modex2_left = np.zeros((nmed-2,ng), dtype=dataType)
    for i in range(nmed-2):
        modex2_left[i]=ini_mode[i+1]
    dmodex2_left = np.zeros((nmed-2,ng), dtype=dataType)
    # region 3
    modex2_right = np.zeros(ng, dtype=dataType)
    dmodex2_right = np.zeros(ng, dtype=dataType)

    amat = np.zeros((nmed, ng, ng))
    bmat = np.zeros((nmed, ng), dtype=dataType)
    cmat = np.zeros((nmed, ng, ng), dtype=dataType)
    for i in range(nmed):
        bmat[i], cmat[i] = med[i].calc_matrix_bc(k_ini)
    dmat = np.zeros((nmed, ng, ng))
    for i in range(nmed):
        for g in range(ng):
            dmat[i, g, g] = med[i].dif(g)

    # region 1
    for j in range(ng):
        b2 = bmat[0][j]
        tmp1 = modex1_left[j]
        tmp2 = 0
        tmp = 0
        if (b2 > 0):
            b = np.sqrt(b2)
            tmp2 = tmp1/np.cos(b*width0)
            tmp = (-tmp2+np.cos(b*width0)*tmp1)*b/(np.sin(b*width0))
        else:
            b = np.sqrt(-b2)
            tmp2 = tmp1*2.0/(np.exp(b*width0)+np.exp(-b*width0))
            tmp = b*tmp1-b*tmp2*np.exp(-b*width0)
        modex0[j] = tmp2
        dmodex1_left[j] = tmp
    # region 2

    for m in range(nmed-2):
        modex1_right[m] = np.dot(np.linalg.inv(cmat[m+1]), np.dot(cmat[m], modex1_left))
        for j in range(ng):
            b2 = bmat[m+1][j]
            # print(b2)
            tmp1 = modex2_left[m,j]
            tmp2 = modex1_right[m,j]
            tmp = 0
            tmp0 = 0
            if (b2 > 0):
            	b = np.sqrt(b2)
            	tmp0 = (tmp1-tmp2*np.cos(b*width1))/(np.sin(b*width1)/(2*b))
            	tmp = (tmp2-tmp1*np.cos(b*width1))/(np.sin(b*width1)/(2*b))
            else:
                b = np.sqrt(-b2)
                array_in = np.array([tmp1,tmp2])
                array_oup = np.dot(
                    np.linalg.inv(np.array([[np.exp(-b*width1)/b, -1/b], [1/b, -np.exp(-b*width1)]])),
                    np.dot(np.array([[-np.exp(-b*width1), 1], [1, -np.exp(-b*width1)]]),array_in))
                tmp=array_oup[0]
                tmp0=array_oup[1]
            dmodex1_right[m,j] = tmp0
            dmodex2_left[m,j] = tmp
    # for reflector region with infinite length
    modex2_right = np.dot(np.linalg.inv(cmat[nmed-1]), np.dot(cmat[nmed-2], modex2_left[nmed-3]))

    for j in range(ng):
        b2 = bmat[2][j]
        # print("here")
        b = np.sqrt(-b2)
        dmodex2_right[j] = -b*modex2_right[j]
    curr_left = np.zeros((nmed-1, ng), dtype=dataType)
    curr_right = np.zeros((nmed-1, ng), dtype=dataType)
    
    curr_left[0] = np.dot(dmat[0], np.dot(cmat[0], dmodex1_left))
    curr_right[0] = np.dot(dmat[1], np.dot(cmat[1], dmodex1_right[0]))

    curr_left[1] = np.dot(dmat[1], np.dot(cmat[1], dmodex2_left[0]))
    curr_right[1] = np.dot(dmat[2], np.dot(cmat[2], dmodex2_right))
    
    curr = curr_left - curr_right

    curr=np.reshape(curr,((nmed-1)*ng))
    oup = np.zeros((nmed-1)*ng+1, dtype=dataType)
    norm = 1-np.sum(xinp[1:])

    for i in range((nmed-1)*ng+1):
        if(i==0):
            oup[i] = norm
        else:
            oup[i] = curr[i-1]
    if(ctrol.converged):
	    print("==>> Saving converged data ...")
	    np.savetxt('.modex0.txt', modex0)
	    np.savetxt('.modex1_left.txt', modex1_left)
	    np.savetxt('.xinp.txt', xinp)
	    np.savetxt('.dmodex1_left.txt', dmodex1_left)
	    np.savetxt('.modex1_right.txt', modex1_right)
	    np.savetxt('.dmodex1_right.txt', dmodex1_right)
	    np.savetxt('.modex2_left.txt', modex2_left)
	    np.savetxt('.dmodex2_left.txt', dmodex2_left)
	    np.savetxt('.modex2_right.txt', modex2_right)
	    np.savetxt('.dmodex2_right.txt', dmodex2_right)
	    np.savetxt('.bmat.txt', bmat)
	    np.savetxt('.cmat.txt', np.reshape(cmat, (nmed*ng*ng)))
    return oup


def save_converged(inp):
    bem(inp)


def get_flx(inp):
    # if(not pplot):
    # 	print("Please turn on the pplot flag on ...")
    # else:
    save_converged(inp)
    # bem(inp)
    print("==>> Reconstructing the neutron flux ...")

    modex0 = np.loadtxt('.modex0.txt', dtype=dataType)
    modex1_left = np.loadtxt('.modex1_left.txt', dtype=dataType)
    dmodex1_left = np.loadtxt('.dmodex1_left.txt', dtype=dataType)
    # region 2
    modex1_right = np.loadtxt('.modex1_right.txt', dtype=dataType)
    dmodex1_right = np.loadtxt('.dmodex1_right.txt', dtype=dataType)
    modex2_left = np.loadtxt('.modex2_left.txt', dtype=dataType)
    dmodex2_left = np.loadtxt('.dmodex2_left.txt', dtype=dataType)
    modex2_right = np.loadtxt('.modex2_right.txt', dtype=dataType)
    dmodex2_right = np.loadtxt('.dmodex2_right.txt', dtype=dataType)
    bmat = np.loadtxt('.bmat.txt', dtype=dataType)
    cmat = np.reshape(np.loadtxt('.cmat.txt', dtype=dataType), (nmed, ng, ng))
    # print(abs(np.dot(cmat[1],modex2_left)-np.dot(cmat[2],modex2_right)))
    # dx = width/div
    dx0 = width0/div0
    dx1 = width1/div1
    dx2 = width2/div2 

    x = dx0*0.5
    mode_tmp_f = np.zeros(ng, dtype=dataType)
    xflx_f = np.zeros((div0, ng+1), dtype=dataType)
    xmodeflx_f = np.zeros((div0, ng+1), dtype=dataType)

    for i in range(div0):
        for j in range(ng):
            b2 = bmat[0, j]
            if(b2 > 0):
                b = np.sqrt(b2)
                mode_tmp_f[j] = -0.5/b*np.sin(b*(width0-x))*dmodex1_left[j]\
                    + 0.5*np.cos(b*(width0-x))*modex1_left[j]\
                    + 0.5*np.cos(-b*x)*modex0[j]
            else:
                b = np.sqrt(-b2)
                mode_tmp_f[j] = 0.5/b*np.exp(-b*(width0-x))*dmodex1_left[j]\
                    + 0.5*np.exp(-b*(width0-x))*modex1_left[j]\
                    + 0.5*np.exp(-b*x)*modex0[j]
        xmodeflx_f[i, 0] = x
        xmodeflx_f[i, 1:ng+1] = mode_tmp_f
        xflx_f[i, 0] = x
        xflx_f[i, 1:ng+1] = np.dot(cmat[0], mode_tmp_f)
        x += dx0
    # region 2
    x = dx1*0.5
    mode_tmp_m = np.zeros(ng, dtype=dataType)
    xflx_m = np.zeros((div1, ng+1), dtype=dataType)
    xmodeflx_m = np.zeros((div1, ng+1), dtype=dataType)
    for i in range(div1):
        for j in range(ng):
            b2 = bmat[1, j]
            if(b2 > 0):
                b = np.sqrt(b2)
                mode_tmp_m[j] = -0.5/b*np.sin(b*(width1-x))*dmodex2_left[j]\
                    + 0.5/b*np.sin(b*x)*dmodex1_right[j]\
                    + 0.5*np.cos(b*(width1-x))*modex2_left[j]\
                    + 0.5*np.cos(-b*x)*modex1_right[j]
            else:
                b = np.sqrt(-b2)
                mode_tmp_m[j] = 0.5/b*np.exp(-b*(width1-x))*dmodex2_left[j]\
                    - 0.5/b*np.exp(-b*x)*dmodex1_right[j]\
                    + 0.5*np.exp(-b*(width1-x))*modex2_left[j]\
                    + 0.5*np.exp(-b*x)*modex1_right[j]

        xmodeflx_m[i, 0] = x+width0
        xmodeflx_m[i, 1:ng+1] = mode_tmp_m
        xflx_m[i, 0] = x+width0
        xflx_m[i, 1:ng+1] = np.dot(cmat[1], mode_tmp_m)
        x += dx1
    # region 3
    x = dx2*0.5
    mode_tmp_r = np.zeros(ng, dtype=dataType)
    xflx_r = np.zeros((div2, ng+1), dtype=dataType)
    xmodeflx_r = np.zeros((div2, ng+1), dtype=dataType)
    tmp=0
    # reflector retion
    for i in range(div2):
        for j in range(ng):
            if(ng == 1):
                b2 = bmat[2]
            else:
                b2 = bmat[2, j]
            b = np.sqrt(-b2)
            mode_tmp_r[j] = -0.5/b*np.exp(-b*x)*dmodex2_right[j]\
            			+0.5*np.exp(-b*x)*modex2_right[j]
        xmodeflx_r[i, 0] = x+width0+width1
        xmodeflx_r[i, 1:ng+1] = mode_tmp_r
        xflx_r[i, 0] = x+width0+width1
        xflx_r[i, 1:ng+1] = np.dot(cmat[2], mode_tmp_r)
        x += dx2
    flx = np.array([xflx_f, xflx_m, xflx_r])
    mflx = np.array([xmodeflx_f, xmodeflx_m, xmodeflx_r])
    # print(cmat[1]==cmat[2])
    # exit()
    plt.xlabel("Distance from the core center ${x} [m]$")

    flux = np.concatenate((flx[0], flx[1], flx[2]), axis=0)
    mflux = np.concatenate((mflx[0], mflx[1], mflx[2]), axis=0)
    # print(flux.shape)
    # print(mflux.shape)
    # exit()

    print("==>> Saving the neutron flux and mode neutron flux ...")
    np.savetxt(".flux.txt", flux)
    np.savetxt(".mode.flux.txt", mflux)
    # print(mflx.shape)
    print("==>> Plotting the neutron flux, [Figure 1]... ")
    plt.figure(1)
    plt.title("Neutron flux distribution of 1-D half core")

    for i in range(0, ng):
        plt.plot(flux[:, 0].real, flux[:, i+1].real,label="Group " + str(i+1))
    plt.xlabel("Distance from the core center ${x} [m]$")

    plt.grid(which='both', axis="both")
    plt.ylabel("Neutron Flux ${\phi}$[x]")

    plt.legend()

    # print("==>> Plotting the mode neutron flux, [Figure 2] ...")
    # plt.figure(2)

    # for g in range(ng):
    #     plt.subplot(ng, nmed, g*nmed+1)
    #     if(g == 0):
    #         plt.title("Core region")
    #     for mode in range(ng):
    #         plt.plot(mflx[0, :, 0].real, np.real(cmat[0, g, mode]*mflx[0, :, mode+1]),
    #                  ".-", linewidth=0.7, markersize=1, label="Mode "+str(mode+1))
    #     plt.plot(mflx[0, :, 0].real, np.real(np.dot(np.transpose(cmat[0, g, :]),
    #                                                 np.transpose(mflx[0, :, 1:ng+1]))), "k-", markersize=1, label="Combination")
    #     plt.grid(which='both', axis="both")
    #     if(ng < 5):
    #         plt.legend(loc='best', prop={'size': 7})
    #     plt.xlabel(" ${x} [m]$")
    #     plt.ylabel("Group "+str(g+1))

    #     plt.subplot(ng, nmed, g*nmed+2)
    #     if(g == 0):
    #         plt.title("Reflector region")
    #     for mode in range(ng):
    #         plt.plot(mflx[1, :, 0].real, np.real(cmat[1, g, mode]*mflx[1, :, mode+1]),
    #                  ".-", linewidth=0.7, markersize=1, label="Mode "+str(mode+1))
    #     plt.plot(mflx[1, :, 0].real, np.real(np.dot(np.transpose(cmat[1, g, :]),
    #                                                 np.transpose(mflx[1, :, 1:ng+1]))), "k-", markersize=1, label="Combination")
    #     plt.xlabel(" ${x} [m]$")
    #     plt.grid(which='both', axis="both")
    #     if(ng < 5):
    #         plt.legend(loc='best', prop={'size': 7})
    fig_file = "mode_flux_"+str(ng)+"g"
    dirr = "../xs/"+str(ng)+"g/"+fig_file
    plt.savefig(dirr, bbox_inches='tight', dpi=600)
    plt.show()
    print("===================================================")

    # print(cmat)

    return flux


def get_initial():
    # ffile = "./cplx_xs/"+str(ng)+"g/interface.txt"
    # interface_flx=np.loadtxt(ffile)

    ffile = "./cplx_xs/"+str(ng)+"g/flux_cbz.txt"
    flux = np.loadtxt(ffile)
    nmesh = int(flux.shape[0]/ng)
    flx = np.zeros((nmesh, ng))
    x = np.zeros(nmesh)
    x = flux[int(nmesh/2)-1:int(nmesh/2)+1, 0]
    for i in range(nmesh):
        for g in range(ng):
            flx[i, g] = flux[nmesh*g+i, 1]
    fac = (flx[int(nmesh/2)-1, -1]+flx[int(nmesh/2), -1])/2
    interface_flx = (flx[int(nmesh/2)-1, :]+flx[int(nmesh/2), :])/(2*fac)

    keff = "./cplx_xs/"+str(ng)+"g/keff.txt"
    keff = np.loadtxt(keff)[0]
    bmat, cmat = fuel.calc_matrix_bc(keff)
    interface_mode_flx = np.dot(np.linalg.inv(cmat), interface_flx)
    iinput = np.zeros(ng+1, dtype=dataType)
    iinput[0] = complex(keff, 0.0)
    iinput[1:] = interface_mode_flx
    # print(interface_mode_flx)
    return iinput


def get_initial_mode():
    ffile = "./cplx_xs/"+str(ng)+"g/flux_bem.txt"
    keff = 1.0
    flux = np.loadtxt(ffile)

    div = 500
    dx = width/div
    x = dx*0.5
    mode_tmp_f = np.zeros(ng, dtype=dataType)
    tmp = 0

    xflx_f = np.zeros((div, ng+1), dtype=dataType)
    xflx = np.zeros((2*div, ng+1), dtype=dataType)
    for i in range(2*div):
        for g in range(ng):
            xflx[i, 0] = x
            xflx[i, g+1] = flux[g*2*div+i, 1]
        x += dx

    xflx_f = xflx[0:div, :]

    keff = np.loadtxt(keff)[0]
    bmat, cmat = fuel.calc_matrix_bc(keff)
    xmodeflx_f = np.zeros((div, ng+1), dtype=dataType)
    for i in range(div):
        xmodeflx_f[i, 1:ng+1] = np.dot(np.linalg.inv(cmat), xflx_f[i, 1:ng+1])
    # print(flux.shape)

    modex2_right = np.zeros(ng, dtype=dataType)
    modex1_left = np.zeros(ng, dtype=dataType)

    modex0 = np.zeros(ng, dtype=dataType)

    dmodex1_left = np.zeros(ng, dtype=dataType)
    dmodex2_right = np.zeros(ng, dtype=dataType)

    # fuel region
    # for i in range(div):
    # 	for j in range(ng):
    # 		b2 = bmat[j]
    # 		if(b2>0):
    # 			b=np.sqrt(b2)
    # 			tmp=-0.5/b*np.sin(b*(width-x))*dmodex1_left[j]
    # 			tmp=tmp+0.5*np.cos(b*(width-x))*modex1_left[j]+0.5*np.cos(-b*x)*modex0[j]
    # 			mode_tmp_f[j] =tmp
    # 		else:
    # 			b=np.sqrt(-b2)
    # 			tmp=0.5/b*np.exp(-b*(width-x))*dmodex1_left[j]
    # 			tmp=tmp+0.5*np.exp(-b*(width-x))*modex1_left[j]+0.5*np.exp(-b*x)*modex0[j]
    # 			mode_tmp_f[j]=tmp
    # 	xmodeflx_f[i,0]=x
    # 	xmodeflx_f[i,1:ng+1] = mode_tmp_f
    # 	xflx_f[i,0]=x
    # 	xflx_f[i,1:ng+1] = np.dot(cmat,mode_tmp_f)
    # 	x+=dx


# for debug
if __name__ == '__main__':
    x0 = np.loadtxt(".xinp.txt",dtype=complex)
    output = bem(x0)

    print(output)
