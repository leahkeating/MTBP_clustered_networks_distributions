import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

fft_data = pd.read_csv("joint_dist_fft.csv")

P_transformed = fft_data['P_re'] + fft_data['P_im']*1j

P_transformed = np.array(P_transformed)
P_transformed = P_transformed.reshape(50,50)

pmf = np.fft.ifft2(P_transformed)
pmf = pmf.real

np.savetxt('simple_contagion_probabilities.csv',pmf)
