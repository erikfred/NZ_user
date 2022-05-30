from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

x = np.arange(15)
print(x)

xmin=4
xmax=10

ix=np.where((x>xmin) & (x<xmax))

x2=x[ix]
print(x2)
