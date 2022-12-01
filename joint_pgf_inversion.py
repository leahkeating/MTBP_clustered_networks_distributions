import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# import the approximation of the transformed data
# which was saved to csv in the R file
fft_data = pd.read_csv("joint_dist_fft.csv")

# we saved the real and imaginary parts in seperate cols 
# so now we need to recombine them
P_transformed = fft_data['P_re'] + fft_data['P_im']*1j

# put this in a numpy array
P_transformed = np.array(P_transformed)
# make it into a matrix
P_transformed = P_transformed.reshape(50,50)

# to recover the pmf, we take the 2D ifft
pmf = np.fft.ifft2(P_transformed)
# keep only the real values
pmf = pmf.real

# save the file so that it can be loaded back into R
np.savetxt('net_science_probabilities.csv',pmf)
