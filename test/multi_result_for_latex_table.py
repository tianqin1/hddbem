# 	  "fvm":[1.21810	,1.21562,	1.21666,	1.21767,	1.21652,	1.21813,	1.21801,	1.21689]}
# niter={"hddbem":[4.2,	4.8,	4.6,	4.5,	5.1,	4.9,	4.7,	6],
# 		"fvm":[6,	11,	12,	30,	36,	16,	22,	46]}
# error=[]
# size=len(group_case)
# for i in range(size):
# 	error.append(abs(keff["hddbem"][i]-keff["fvm"][i])/keff["fvm"][i])
# formatter1="{:.5f}".format
# formatter2 = "{:.2e}".format
# for i in range(size):
# 	print(str(group_case[i])+'-group'+" & "+ formatter1(keff["hddbem"][i])+" & "+str(niter["hddbem"][i])+" & "+formatter1(keff["fvm"][i])+" & "+str(niter["fvm"][i])+" & "+formatter2(error[i])+"\\\\")
# # print(keff["hddbem"][0])

import numpy as np
group_case = [3, 8, 11, 15, 20, 25, 26, 27]
size = len(group_case)
location="../src/"
keff = {"hddbem": [1.60088,	1.60453,	1.60459,	1.60278,	1.60133,	1.59864,	1.59833,	1.59904],
        "fvm": [1.60088,	1.60454,	1.60460,	1.60278,	1.60134,	1.59865,	1.59833,	1.59903]}
niter = {"hddbem": [3.8,	3.6,	4.4,	4.4,	5.5,	5.3,	5.1,	17],
         "fvm": [10,	10,	29,	28,	18,	31,	38,	29]}
for i in range(size):
	g=group_case[i]
	tmp=np.loadtxt(location+str(g)+"g.info")
	# print(tmp.shape[0])
	niter["hddbem"][i]=sum(tmp[:,0])/tmp.shape[0]

error = []
for i in range(size):
    error.append(abs(keff["hddbem"][i]-keff["fvm"][i])/keff["fvm"][i])
formatter1 = "{:.5f}".format
formatter2 = "{:.2e}".format
formatter3 = "{:.1f}".format
for i in range(size):
    if(i != 6):

        print(str(group_case[i])+'-group'+" & " + formatter1(keff["hddbem"][i])+" & "+formatter3(niter["hddbem"][i]) +
              " & "+formatter1(keff["fvm"][i])+" & "+str(niter["fvm"][i])+" & "+formatter2(error[i])+"\\\\")
# print(keff["hddbem"][0])
