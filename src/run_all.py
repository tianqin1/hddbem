import numpy as np
import os
import sys

for i in range(41,108):
	os.system("touch .ng_var")
	f = open('.ng_var', 'w')
	f.write(str(i)+' 0')  # python will convert \n to os.linesep
	f.close()
	os.system("python Main.py")
	os.system("rm .ng_var")
