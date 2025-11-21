

This repository contains code and download scripts, that will allow you to make calculations for 1 year (2001) like in article "__LAGRANGIAN MODELLING OF THE BLACK SEA TURBOT EARLY LIFE STAGES: TRANSPORT AND DISTRIBUTION"__ by Krasilnikov D.S., Mizyuk A.I., Khanaychenko A.N., Klimova T.N., Bagaev A.V.

__Disclaimer__: We have used and modified the OceanParcels model code (http://parcels-code.org).

## Plan:

1. Download Input files using bash script download_input.sh: `bash download_input.sh`. As a result in Input./ you will have 14 files
2. Install conda environment from environment.yml using one of the ways below

* `conda env create -f environment.yml` by default it will be named bst_parcles
* `conda env create -f environment.yml -n my_env_name` to specify desired name

3. Activate your environment `conda activate my_env_name `

4. Launch bst_calc.py: `python bst_calc.py`

5. After execution You will have 3 files in Output/. These 3 files are:

* "Parcels_2001.zarr" contains the trajectories and temperatures of the particles that started according to those indicated at the beginning of the file. bst_calc.py conditions. This is prettty raw file and it shouldn't be used to any analysis.

* '"Trajectories_2001.nc"" contains trajectories for those parcels that started in good temperature spawning conditions and were twice checked in order to exclude travelling on land as raw OceanParcels simulations sometimes do :). Suitable for analysis.

* "Start_and_Finish_2001.nc " contains a summary table that stores information about the starting and final coordinates, time and temperature of each particle simulating an bst egg. It also contains information about the stage to which the egg lived, whether it had gastrulation, and for what reasons it died (if the stage is < 0). Suitable for analysis.
