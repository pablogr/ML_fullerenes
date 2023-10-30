# ML_fullerenes
This contains the code necessary for the forecast of eigenvalue renormalisations of fullerenes, as well as the code used to perform the calculations of the corresponding article.

The repository consists of three parts:
* CODE_FOR_CALCULATIONS_OF_THE_PAPER
* CODE_FOR_FORECASTING
* FULLERENE_XYZ_FILES

CODE_FOR_CALCULATIONS_OF_THE_PAPER
This includes the code which was used to perform the regression-based calculations whose results are presented in the article 'Electron-vibrational renormalization in fullerenes through ab initio and machine learning methods', by Pablo Garcıa-Risueno, Eva Armengol, Angel Garcia-Cerdana, David Carrasco-Busturia and Juan M. Garcia-Lastra. 

CODE_FOR_FORECASTING
This includes the code which can be used to make a forecast of the renormalization of the electronic eigenvalues of an arbitrary fullerene. One just need to enter a file with the results of a given DFT ground state calculation and the code outputs an estimate of the electron-vibrational renormalization.

Both Code_for_forescasting and Code_for_calculations_of_the_paper require downloading Python a setting a virtual environment ("venv"). In addition, they require installing some popular libraries, like pandas, numpy and scikit-learn. Both are executed from their respective main.py files, and the user can set the input parameters in the respective module_io.py files.

FULLERENE_XYZ_FILES
This contains the files with the atomic coordinates (in .xyz format) of the fullerenes whose calculations are presented in the mentioned article. In addition, it contains the input files and pseudopotentials to be used with the ab initio DFT simulation code Quantum Espresso. This lies in subfolder called 'info_for_calculations'. This subfolder also includes Fortran code which can be used to perform actual frozen-phonon calculations of the renormalizations (this is no forecasts, but actual accurate, numerically heavy ab initio calculations). How to do it is better explained in the repository found at https://github.com/pablogr/Correction_of_anticrossings_frozen_phonon; the latter repository contains a detailed tutorial which explains how to perform the Quantum Espresso calculations, and eventual corrections of the renormalizations as explained in the article 'Frozen-phonon method for state anticrossing situations and its application to zero-point motion effects in diamondoids' [https://journals.aps.org/prb/abstract/10.1103/PhysRevB.108.125403]. Note that such corrections do not need to be calculated in most of the fullerenes. Nevertheless, the mentioned detailed tutorial includes information on how the Quantum Espresso calculations of actual renormalizations must be done (even if no correction is performed).
The functions of the file info_for_calculations/Fortran_code/extract_data_fullerenes.f90 were used to do the post-processing of the output files of Quantum Espresso which provided the numbers which were later used as input for the ML algorithms used in the code of CODE_FOR_CALCULATIONS_OF_THE_PAPER and CODE_FOR_FORECASTING.




