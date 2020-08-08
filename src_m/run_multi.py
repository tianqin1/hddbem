import numpy as np
import os
seeds=np.arange(10)
# print(seeds)
for i in range(len(seeds)):
	sseed=seeds[i]
	tmp=np.array([sseed,0])
	np.savetxt(".seed",tmp)
	cmd='python main.py'
	os.system(cmd)
