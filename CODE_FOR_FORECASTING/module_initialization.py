import numpy as np

import module_io
import pandas as pd
from os import makedirs, path, remove

#--------------------------------------------------------------------------------------------------------------

print("\n\n\n=================================================================================================")
print("   PROGRAM FOR FORECASTING OF FULLERENES FROZEN-PHONON RENORMALIZATIONS USING MACHINE LEARNING ")
print("=================================================================================================")

#--------------------------------------------------------------------------------------------------------------

def standardize_regressors(df_in0):
    '''This function standardizes the regressors of the input dataframe.'''

    li_fields_not_to_standardize = ['Name', 'Symmetry (generic)', 'Symmetry (specific)', 'HOMO ren (output)', 'LUMO ren (output)', 'Gap ren (output)']
    li_fields_to_standardize = list(df_in0)
    list_reuse = []
    for field in li_fields_not_to_standardize:
        if (field in li_fields_to_standardize):
            li_fields_to_standardize.remove(field)
            list_reuse.append(field)

    df_in_standardized = df_in0[ list_reuse ]

    for field in li_fields_to_standardize:
        my_array = df_in0[field].to_numpy()
        my_avg = np.mean(my_array); my_stdev = np.std(my_array)
        my_array = ( my_array - my_avg ) / my_stdev
        df_in_standardized[field] = my_array
        #df_in_standardized[field] = pd.Series( my_array )
        #df_in_standardized.loc[:,field] = pd.Series(my_array)
        #for i in range(len(my_array)): df_in_standardized.loc[i, field] = my_array[i]
        #df_in_standardized[field] = (df_in0[field] - my_avg )/my_stdev
        #df_in_standardized.loc[:,field] = df_in0.loc[:,field]
        #df_in_standardized[field] = ( df_in_standardized[field] - my_avg) / my_stdev


    del my_array; del li_fields_to_standardize; del li_fields_not_to_standardize; del list_reuse

    return df_in_standardized

#--------------------------------------------------------------------------------------------------------------

def initialization(  config  ):
    '''
    :return: field_out_to_read (string): Either "HOMO ren (output)", "LUMO ren (output)" or  "Gap ren (output)"
    :return: regressor_fields (list of strings): The list of regressors to use in each case.
    '''


    print("\n\n\n -------------- Now calculating the "+ config.quantity_to_forecast+" RENORMALIZATION  ----------------------------\n")

    if ( config.regressor_fields_regression == config.regressor_fields_ml):
        print("   The regressors are",config.regressor_fields_regression)
    else:
        print("   The",len( config.regressor_fields_regression),"regressors for OLS regression are", config.regressor_fields_regression)
        print("   The", len(config.regressor_fields_ml), "regressors for the ML algorithm (",config.ml_method,") are", config.regressor_fields_ml)
    print("   The chosen calculation method is "+config.ml_method+".")
    print("   The presented results will be the average of "+str((module_io.N_essais))+" random selections of the dataset.")
    if ((module_io.standardize_regressors)): txt=""
    else: txt="not "
    print("   The regressors (input columns) were " + txt  + "standardized.")

    print("   The input and output data are read from",module_io.filepath_raw_input_data)
    df_in = pd.read_csv(module_io.filepath_raw_input_data, header=0, engine='c')

    if module_io.standardize_regressors: df_in = standardize_regressors(df_in)

    N_data_test = 1
    N_data_training = len(df_in)


    print("   Number of data for training:    ",N_data_training,"\n\n" )
    #print("   Number of data for test:        ",N_data_test," (","{:.3}".format( 100*N_data_test/(len(df_in) )),"%)\n\n")


    return  df_in, config.field_out_to_read, config.regressor_fields_regression,  config.regressor_fields_ml, config.different_regressors, config.ml_method, module_io.N_essais,  N_data_training, N_data_test,  module_io.verbosity

#--------------------------------------------------------------------------------------------------------------------------------------

