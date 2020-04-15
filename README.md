# mb_bioequivalence
This repository contains code corresponding to the manuscript "Efficient model-based bioequivalence Testing".

The R-file "be_nlmem_parallel.R" contains code to apply the two model-based approaches, that is the new optimal procedure, and the model-based TOST, to two examplatory datasets, namely dataset_1.txt and dataset_2.txt. These datasets have been generated usind a parallel design and assuming H_0 and H_1, respectively. For a more detailed description see the comments in the code.
As the (model-based) version for crossover designs requires additionally to R the software Monolix, it is not uploaded here, but can be sent upon request (as well as the NCA-based methods for crossover designs). 
The R-file "be_nca_parallel.R" contains code to apply the two NCA-based approaches to the two datasets described above. Estimates for log AUC and log Cmax are obtained by NCA and the tests are conducted in the same manner as for the estimates obtained by fitting NLMEM.
Further datasets can also be provided upon request.
