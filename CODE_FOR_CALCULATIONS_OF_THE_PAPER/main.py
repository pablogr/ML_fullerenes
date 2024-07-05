'''
========================================================================================================================
                 PROGRAM USED TO PERFORM THE REGRESSION CALCULATIONS PRESENTED IN THE ARTICLE
      'ELECTRON-VIBRATIONAL RENORMALIZATION IN FULLERENES THROUGH AB INITION AND MACHINE LEARNING METHODS'
                                       BY PABLO GARCIA-RISUENO ET AL.
========================================================================================================================

This program checks the effectiveness of several Machine Learning methods to forecast the renormalizations of electronic
eigenvalues of fullerenes due to electron-vibrational interaction (which are obtained from frozen phonon calculations).
It reads the corresponding datafile, splits it into training and test datasets and calculates the average absolute error of
the forecast in the test dataset.
All the input parameters for this calculation must be specified in the module_io.py file. There you can specify arbitrary lists of
regressors (among those available in the .csv files of the Data_ab_inition directory). You can also specify different
machine-learning methods and their hyperparameters in the module_training.py.

Code fully written by Pablo Garcia Risueno [garciaDOTrisuenoATgmailDOTcom] between 2021 and January 2023.
'''


from sys import modules; modules[__name__].__dict__.clear() #clear all variables in memory
import module_io
from module_initialization import initialization, initialize_error_files, middle_message, final_message
from module_training import select_training_test_and_BT_datasets, linear_regression, ml_regression, update_df_errors
import gc
from module_plots import plot_feature_importance, plot_feature_importance_alldata_MDI


'''
# Block for plots. Uncomment it and specify the appropriate paths below if you want to do plots similar to those appearing in the article. 

from module_plots import box_plot_quantiles, integrated_plot
for quantity in ["HOMO"]:
    box_plot_quantiles(  "ml_test_errors_RF_"+quantity+"__fromresid_250"  )
    box_plot_quantiles("ml_test_errors_NN_" + quantity + "__fromresid_100")
    box_plot_quantiles("ml_test_errors_KNN_" + quantity + "__fromresid_22")
    integrated_plot( f"{module_io.execution_directory}Results/ml_training_errors_NN_"+quantity+"__fromresid_100.csv", f"{module_io.execution_directory}Results/ml_test_errors_NN_"+quantity+"__fromresid_100.csv", quantity, "NN", True )
    integrated_plot( f"{module_io.execution_directory}Results/ml_training_errors_RF_"+quantity+"__fromresid_250.csv", f"{module_io.execution_directory}Results/ml_test_errors_RF_"+quantity+"__fromresid_250.csv", quantity, "RF", True )
    integrated_plot( f"{module_io.execution_directory}Results/ml_training_errors_KNN_"+quantity+"__fromresid_22.csv", f"{module_io.execution_directory}Results/ml_test_errors_KNN_"+quantity+"__fromresid_22.csv", quantity, "KNN", True )
exit(0)
'''


for quantity in ["gap"]: #["HOMO","LUMO","gap"]
    for forecast_of_residuals in [True]:  #[True,False]
        for method in ["RF"]: #["RF","NN", "KNN","GB"]

            my_configuration = module_io.Configuration(quantity,method,forecast_of_residuals)

            # Initialization
            df_in, field_out_to_read, regressor_fields_lr, regressor_fields_ml, all_cols, different_regressors, ml_method, N_essais, \
            BT_data, N_data_training, N_data_test, verbosity = initialization( my_configuration, quantity )

            # Sweep in different random selections of rows of the dataset
            for param_sw in my_configuration.param_sweep:

               if (method=="RF"): plot_feature_importance_alldata_MDI(my_configuration, field_out_to_read, regressor_fields_lr, regressor_fields_ml, param_sw, N_essais, forecast_of_residuals)
               exit(0)# xxx

               sum_error_meV = 0.0; sum_error_perc = 0.0; sum_avg_abs_error_val_ols=0.0
               df_errors_train =  module_io.errors_df0.copy()
               df_errors_test  =  module_io.errors_df0.copy()
               features_weights = dict(zip(my_configuration.regressor_fields_ml,[0.0]*(len(my_configuration.regressor_fields_ml) )))

               for i in range(1,N_essais+1):

                      if ( (N_essais > 99) and ((i%50)==0) ): print("Essai",i,"out of ",N_essais)

                      # Random selection of datasets for test, training and (eventually) backtesting:
                      filepath_training_data, filepath_test_data, filepath_BT_data = select_training_test_and_BT_datasets( df_in, N_data_training, N_data_test, i, BT_data  )


                      # Calculation of residuals of the selected datasets
                      if (my_configuration.forecast_of_residuals):
                         linear_regression( field_out_to_read, regressor_fields_lr, regressor_fields_ml, filepath_training_data, filepath_test_data, filepath_BT_data, verbosity, different_regressors, BT_data  )
                         filepath_training_data_ml = filepath_training_data.replace(".csv","-with_residuals.csv")
                         filepath_test_data_ml = filepath_test_data.replace(".csv","-with_residuals.csv")
                      else:
                         filepath_training_data_ml = filepath_training_data
                         filepath_test_data_ml     = filepath_test_data

                      error_meV, error_perc, avg_abs_error_val_ols, df_err_train_1, df_err_test_1, features_weights = ml_regression( my_configuration, param_sw, filepath_training_data_ml, filepath_test_data_ml, features_weights ) #( ml_method, field_out_to_read, regressor_fields_ml, filepath_training_data_ml, filepath_test_data_ml, param_sw, forecast_of_residuals, verbosity )

                      if (i==1): fp_train, fp_test = initialize_error_files(my_configuration, method, quantity, forecast_of_residuals, param_sw, all_cols)
                      df_errors_train = update_df_errors( i, df_errors_train, df_err_train_1, fp_train, my_configuration.forecast_of_residuals, all_cols, 'train' )
                      df_errors_test  = update_df_errors( i, df_errors_test,  df_err_test_1, fp_test, my_configuration.forecast_of_residuals, all_cols, 'test' )
                      if ((i%100)==0): gc.collect()

                      sum_error_meV += error_meV; sum_error_perc += error_perc;
                      if (my_configuration.forecast_of_residuals): sum_avg_abs_error_val_ols += avg_abs_error_val_ols
                      if ( (verbosity != 0) and ((i%200)==0)): print( i, " M.e.%:",sum_error_perc/i, "(",sum_error_meV/i,"meV)")

               if ((method=="RF") and (module_io.plot_FI_MDI_separating_train)):
                  for regressor_name in my_configuration.regressor_fields_ml:
                      features_weights[regressor_name]  /= N_essais
                  plot_feature_importance(features_weights, my_configuration.quantity_to_forecast, my_configuration.ml_method, forecast_of_residuals)

               sum_error_meV /= N_essais; sum_error_perc /= N_essais;
               if (my_configuration.forecast_of_residuals):  sum_avg_abs_error_val_ols/= N_essais # Average of the errors

               middle_message( my_configuration, param_sw, sum_error_meV, sum_error_perc, sum_avg_abs_error_val_ols)

               #if (module_io.store_ml_errors): integrated_plot( fp_train, fp_test, quantity, method, forecast_of_residuals )

            final_message(my_configuration)
            if ((forecast_of_residuals) and (module_io.store_ml_errors) ):
                if (N_essais<100):
                    df_errors_train.to_csv(fp_train, index=False)
                    df_errors_test.to_csv(fp_test,   index=False)
                else:
                    df_errors_train.to_csv(fp_train, mode='a', header=False, index=False)
                    df_errors_test.to_csv( fp_test, mode='a', header=False, index=False)
            del df_errors_train; del df_errors_test


print("====================== All calculations finished satisfactorily ===============================================")


