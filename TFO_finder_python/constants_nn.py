''' The data neede to compute my RNA-DNA2 nearest neighbor model
'''
import numpy as np
# assume a series of physiological conditions at the moment
VAL_PH = float(7.0)
CONC_TFO = float(12.0)
DUPLEX_CONC = 0.05
SPEC = "RNA"
# Read the coefficients of the model
NN_DH = np.array([-10.95, -5.73, -6.44])
NN_DG = np.array([-1.891, -0.758, -0.331, 2.646])
COEF_PH = np.array([0.893, -0.005])
COND = [VAL_PH, CONC_TFO, DUPLEX_CONC]
