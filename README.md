# Continous_SIF
This contains the code for generating the CSIF dataset and the NN parameters
## NN parameters
We tested a variety of NN architecture (1-3 layers and 2-9 neurons for each layer), the best formed NN was selected (1 layer 5 neurons). The NN_parameter folder contains the weights and biases for each layer.
This folder also contains the two files required for the normalization of reflectances (.nc files).
## preprocessing of MCD43C4
We fellowed Zhang et al (2017) and generate the 4-day reflectance. Three steps corresponds to the three python files:
1. aggreated the daily reflectance files to 4day using their median values.
2. get a mean seasonal cycle as a reference.
3. reconstruct the SIF for each year.
## generate the CSIF product
Several products are needed to generate the CSIF (clear-day and all-daily)
>4day refletance dataset from the previous step.
>mean and standard deviation of reflectance for each band (in NN parameter folder).
>solar zenith angle calcualted based on the solar time and latitude.
>BESS daily PAR product (only required if all-daily SIF is calculated).
>DEM product to calcualt the clear-sky PAR.
The clear-inst and clear-daily SIF data is first calculated using the neuron network. The all-daily SIF is then calcualted using the clear-inst SIF and the daily PAR from BESS.

