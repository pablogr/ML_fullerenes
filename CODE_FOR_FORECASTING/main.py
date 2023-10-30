'''
========================================================================================================================
                 PROGRAM USED TO FORECAST THE RENORMALIZATION OF ELECTRONIC EIGENVALUES OF FULLERENES
      See 'ELECTRON-VIBRATIONAL RENORMALIZATION IN FULLERENES THROUGH AB INITION AND MACHINE LEARNING METHODS'
                                       BY PABLO GARCIA-RISUENO ET AL.
========================================================================================================================

This program forecasts the renormalizations of electronic eigenvalues of fullerenes using ab initio data calculated with
Quantum Espresso. If you want to have an estimate of the renormalization of electronic eigenvalues in a fullerene,
run a ground state calculation (preferably using Quantum Espresso with the provided pseudopotential files), then
extract the information of the regressors from its output file, write them to a file (called input_regressors_molecule_to_forecast.csv,
see the example file with that name) and run this program.

Code fully written by Pablo Garcia Risueno [garciaDOTrisuenoATgmailDOTcom], January 2023. Contact me if you have any doubt.
'''

from sys import modules; modules[__name__].__dict__.clear() #clear all variables in memory
import module_io
from module_initialization import initialization
from module_training import ols_multiregression, ml_regression
import numpy as np
import gc


forecast_of_residuals = True
filepath_training_data = f"{module_io.execution_directory}training_data-with_residuals.csv"


for method in [ "RF","NN", "KNN" ]: # ["RF","NN", "KNN"]
    for quantity in [ "HOMO", "LUMO","gap" ]: # ["HOMO","LUMO","gap"]

            array_ml_residual_forecasts = np.zeros(module_io.N_essais)
            my_configuration = module_io.Configuration(quantity,method,forecast_of_residuals)

            # Initialization
            df_in, field_out_to_read, regressor_fields_regression, regressor_fields_ml, different_regressors, ml_method, N_essais, \
            N_data_training, N_data_test, verbosity = initialization( my_configuration )

            for param_sw in my_configuration.param_sweep:

               my_forecast_lr = ols_multiregression(quantity, field_out_to_read, regressor_fields_regression,
                                                    regressor_fields_ml, module_io.filepath_raw_input_data,
                                                    filepath_training_data, verbosity, different_regressors)

               for i in range(N_essais):
                      array_ml_residual_forecasts[i] = ml_regression( my_configuration, param_sw, filepath_training_data, module_io.filepath_regressors_molecule_to_forecast, my_forecast_lr )
                      if (((i+1)%100)==0): gc.collect()

               print("   * The", my_configuration.quantity_to_forecast + " renormalization forecasted through machine learning (" + module_io.methname[my_configuration.ml_method] + ") on top of linear regression is", "{:.3}".format( my_forecast_lr + np.mean(array_ml_residual_forecasts)), "+/-", "{:.3}".format( 1.96*np.std(array_ml_residual_forecasts)) )
               #print(array_ml_residual_forecasts)

print("\n====================== All calculations finished satisfactorily ===============================================")


