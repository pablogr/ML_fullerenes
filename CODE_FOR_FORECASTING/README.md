# ML_fullerenes_renormalization_code_for_forecasting

===============================================
                 PROGRAM TO FORECAST THE RENORMALIZATION OF ELECTRONIC EIGENVALUES OF FULLERENES
      See 'ELECTRON-VIBRATIONAL RENORMALIZATION IN FULLERENES THROUGH AB INITION AND MACHINE LEARNING METHODS'
                                       BY PABLO GARCIA-RISUENO ET AL.
===============================================

This program forecasts the renormalizations of electronic eigenvalues of fullerenes using ab initio data calculated with
Quantum Espresso. If you want to have an estimate of the renormalization of electronic eigenvalues in a fullerene,
run a ground state calculation, preferably using Quantum Espresso with the provided pseudopotential files. Example input files
are scf.in (for an ordinary PBE calculation), and scf-B3LYP.in (for the B3LYP run). Note that a geometry optimization
of the atomic position of your fullerenes should be performed first). 
After having run the ground state calculation using Quantum Espresso, a file called eigenval.xml will be produced. Copy it
to the directory of the Python code shared here. Then run the code of module_writing_input_data.py. It will generate the file
called input_regressors_molecule_to_forecast.csv, which contains information of the regressors (note that if you ran the B3LYP 
calculation then you will have to manually modify input_regressors_molecule_to_forecast.csv to include the B3LYP HOMO-LUMO gap).
After having input_regressors_molecule_to_forecast.csv you can run the main program, which will provide you the forecast of the
renormalizations.

Code fully written by Pablo Garcia Risueno [garciaDOTrisuenoATgmailDOTcom], January 2023. Contact me if you have any doubt.
