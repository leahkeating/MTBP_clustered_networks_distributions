# This GitHub contains the code used in the following preprint

## Below we give a description of the files:

- **functions.R:** This file contains many of the functions required to run the other scripts.
- **joint_distributions_EECC_and_correlation.R:** This file contains the code for generating the joint distribution of cascade size and cumulative depth. Here we calculate the Pearson's correlation coefficient and the edge-disjoint edge clique cover (EECC).
- **joint_pgf_inversion.py:** This file contains the inversion function for the joint pgf, used with ``joint_distributions_EECC_and_correlation.R``. 
- **EECC_on_empirical_networks.R:** This file contains the code for finding the cascade-size distributions for real-world networks, this uses the EECC.
