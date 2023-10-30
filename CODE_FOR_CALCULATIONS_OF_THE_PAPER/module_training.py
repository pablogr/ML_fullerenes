
import module_io
import pandas as pd
import numpy as np
pd.set_option('display.max_columns', 20);# pd.set_option('display.max_rows', 99999)

#-------------------------------------------------------------------------------------------------------------------------------------

def merge_lists(li1,li2):
    '''This function merges two lists removing the duplicates'''
    li3 = li1+li2
    li3=pd.DataFrame(li3)
    li3=li3.drop_duplicates()
    li3=li3[0].values.tolist()
    return li3

#-------------------------------------------------------------------------------------------------------------------------------------

def select_training_test_and_BT_datasets( df_in, N_data_training, N_data_test, random_seed_seed, BTdata  ):
   '''This function randomly selects parts of the input dataset as the datasets for training, test and backtesting.
   INPUT:
   :param filepathin: (string) Path of the file whose data will be split into training, test and BT parts.
   :param random_seed_seed: (integer) Number used to define the random seed to select the data.
   :param BTdata: (boolean) True if part of the whole dataset must be kept for backtesting; false otherwise.
   OUTPUT: Paths of the 3 files which contain the data on training, test and BT.'''

   from sklearn.utils import shuffle

   np.random.seed( int( random_seed_seed * 1234)  )  # To avoid random variation of the result on the same input dataset.

   df_aux = shuffle(df_in)
   df_training = df_aux[0:N_data_training]
   df_test     = df_aux[N_data_training:N_data_training+N_data_test]
   if BTdata: df_BT       = df_aux[N_data_training+N_data_test:len(df_in)]

   filepath_training_data = f"{module_io.execution_directory}training_data.csv"
   filepath_test_data = f"{module_io.execution_directory}test_data.csv"
   filepath_BT_data = f"{module_io.execution_directory}backtesting_data.csv"
   df_training.to_csv(  filepath_training_data ,index=False )
   df_test.to_csv( filepath_test_data, index=False)
   if BTdata:
       df_BT.to_csv( filepath_BT_data, index=False)
       del df_BT

   del df_in; del df_aux; del df_training; del df_test

   return filepath_training_data, filepath_test_data, filepath_BT_data

#-------------------------------------------------------------------------------------------------------------------------------------

def average_files(field_out_to_read, filepathin ):
    '''This function takes the files with the data of the residuals for training and test which correspond to
    different sets of regressors (for an ensemble of two methods, each with its own list of regressors) and writes
    the corresponding file by averaging both.'''

    df1 = pd.read_csv(filepathin.replace(".csv","_1.csv"),header=0)
    df1.set_index("Name",inplace=True)
    df1 = df1.rename(columns={field_out_to_read + "_residual": field_out_to_read + "_residual_1"})
    if (field_out_to_read + "_forecasted" in list(df1)):
        df1 = df1.rename(columns={field_out_to_read + "_forecasted": field_out_to_read + "_forecasted_1"})

    df2 = pd.read_csv(filepathin.replace(".csv", "_2.csv"), header=0)
    df2.set_index("Name", inplace=True)
    if (field_out_to_read + "_forecasted" in list(df2)):
        df2 = pd.DataFrame(df2[ [field_out_to_read + "_residual",field_out_to_read + "_forecasted"] ])
        df2 = df2.rename(columns={field_out_to_read + "_residual": field_out_to_read + "_residual_2"})
        df2 = df2.rename(columns={field_out_to_read + "_forecasted": field_out_to_read + "_forecasted_2"})
    else:
        df2 = pd.DataFrame(df2[field_out_to_read + "_residual"])
        df2 = df2.rename(columns={field_out_to_read + "_residual": field_out_to_read + "_residual_2"})

    df1 = pd.concat([df1,df2],axis=1)
    df1[field_out_to_read+ "_residual"] = ( df1[field_out_to_read+ "_residual_1"] + df2[field_out_to_read+ "_residual_2"] )/2
    if (field_out_to_read + "_forecasted_1" in list(df1)):
        df1[field_out_to_read + "_forecasted"] = (df1[field_out_to_read + "_forecasted_1"] + df2[field_out_to_read + "_forecasted_2"]) / 2

    df1.to_csv(filepathin,index=True)

    del df1; del df2; del filepathin

    return

#-------------------------------------------------------------------------------------------------------------------------------------

def linear_regression( field_out_to_read, regressor_fields, regressor_fields_to_save, filepathin_training_data, filepathin_test_data, filepathin_BT_data, verbose, different_regressors, BT_data ):
    ''' This is a wrapper which calls to slightly different functions if the linear regression must be done with an ensemble of
    two methods (i.e. with two lists of regressors) or not.

    This function finds the residuals of the multiple regression (through Ordinary Least Squares) and writes it to
    files together with the regressors for all three datasets (training, test, backtesting).

    :param field_out_to_read: (string) must be either "HOMO", "LUMO" or "gap".
    :param regressor_fields: Names (labels) of the regressors to be used in the OLS (linear) regression
    :param regressor_fields_to_save: Names (labels) of the regressors to be stored in the output file
    :param filepathin_training_data: Path to the file containing training data. This file, as well as the two files below, must
           contain at least the columns which correspond to the regressors, the name of the molecule and the output value (HOMO_ren, LUMO_ren or GAP_ren).
    :param filepathin_test_data:     Path to the file containing test data
    :param filepathin_BT_data:       Path to the file containing backtesting data.
    :param verbose:  (int) True if enhanced information is to be printed to screen.
    :param different_regressors: (Boolean) False if the regressors which are used in the regression and in the ML part are the same; true otherwise.
    :param BT_data: (Boolean) True if backtesting data must be generated; False otherwise.
    '''

    if not (module_io.use_ensemble_for_ols):

        ols_multiregression( field_out_to_read, regressor_fields, regressor_fields_to_save, filepathin_training_data, filepathin_test_data, filepathin_BT_data, verbose, different_regressors, BT_data )

    else: # The linear regression must be performed as an ensemble of two methods, whose 2 outputs are averaged.

        ols_multiregression(field_out_to_read, regressor_fields[0], regressor_fields_to_save, filepathin_training_data,
                            filepathin_test_data, filepathin_BT_data, verbose, different_regressors, BT_data,"_1")
        ols_multiregression(field_out_to_read, regressor_fields[1], regressor_fields_to_save, filepathin_training_data,
                        filepathin_test_data, filepathin_BT_data, verbose, different_regressors, BT_data,"_2")
        average_files(field_out_to_read, filepathin_training_data.replace(".csv","-with_residuals.csv") )
        average_files(field_out_to_read, filepathin_test_data.replace(".csv","-with_residuals.csv") )

    return

#-------------------------------------------------------------------------------------------------------------------------------------

def ols_multiregression( field_out_to_read, regressor_fields, regressor_fields_to_save, filepathin_training_data, filepathin_test_data, filepathin_BT_data, verbose, different_regressors, BT_data, suffix="" ):
    ''' This function finds the residuals of the multiple regression (through Ordinary Least Squares) and writes it to
    files together with the regressors for all three datasets (training, test, backtesting).

    :param field_out_to_read: (string) must be either "HOMO", "LUMO" or "gap".
    :param regressor_fields: Names (labels) of the regressors to be used in the OLS (linear) regression
    :param regressor_fields_to_save: Names (labels) of the regressors to be stored in the output file
    :param filepathin_training_data: Path to the file containing training data. This file, as well as the two files below, must
           contain at least the columns which correspond to the regressors, the name of the molecule and the output value (HOMO_ren, LUMO_ren or GAP_ren).
    :param filepathin_test_data:     Path to the file containing test data
    :param filepathin_BT_data:       Path to the file containing backtesting data.
    :param verbose:  (int) True if enhanced information is to be printed to screen.
    :param different_regressors: (Boolean) False if the regressors which are used in the regression and in the ML part are the same; true otherwise.
    :param BT_data: (Boolean) True if backtesting data must be generated; False otherwise.
    '''

    from sklearn import linear_model

    # PART 0: INITIALIZATION -------------------------------------------------------------------------------------------

    regressor_fields_to_save1 = merge_lists(regressor_fields,regressor_fields_to_save)
    df_raw_data = pd.read_csv( filepathin_training_data, header=0, engine='c', usecols=["Name"] + regressor_fields_to_save1 + [field_out_to_read])
    y_frc   = pd.DataFrame(index=df_raw_data["Name"], columns=[ field_out_to_read + "_forecasted"])   # forecasted value of the independent variable
    y_resid = pd.DataFrame(index=df_raw_data["Name"], columns=[ field_out_to_read + "_residual"])
    df_raw_data.set_index("Name",inplace=True)
    X = df_raw_data[ regressor_fields ]
    X_to_save = df_raw_data[regressor_fields_to_save]
    y = df_raw_data[ field_out_to_read  ]


    # PART 1: ACTUAL REGRESSION (with multiple regressors; finding alpha and beta) -------------------------------------

    regr = linear_model.LinearRegression()
    resu = regr.fit(X.values, y.values)
    # resu_coefs = regr.coef_ ; print(resu_coefs)

    if (verbose == 2):
        print("\n     * Now calculating the residuals of the <<"+field_out_to_read+">> variable using linear regression with multiple ("+str(len(regressor_fields ))+" regressors. The regressors are:\n       "+str(regressor_fields))
        print("          Our input dataset consists of ", len(X), " data (rows).")
        print("\n          Number of regressors,R^2 of the multiple regression")
        print("          ",len(regressor_fields),",",resu.score(X,y),"\n")
        print("  TRAINING DATA: \n")
        print("i,  Yobser,  Yforcst,  Diff(meV),  Diff(%) ")
    avg_error_meV_tr = 0.0; avg_error_Perc_tr  = 0.0

    for i in range(len(X)):
        X_trial = X.iloc[i]
        y_trial = y.iloc[i]
        y_forecasted = resu.predict(np.array([ X_trial ]))[0]
        y_frc.iloc[i][field_out_to_read + "_forecasted"] = y_forecasted
        y_resid.iloc[i][field_out_to_read + "_residual"] = y_trial - y_forecasted
        if (verbose == 2):
             print(i, ",  ",y_trial, ",  ","{:.3}".format( y_forecasted), ",  ", "{:.3}".format( (y_trial - y_forecasted)), ",  ", "{:.4}".format( 100*(y_trial - y_forecasted)/y_trial ) )
        avg_error_meV_tr  += abs((y_trial - y_forecasted)); avg_error_Perc_tr  += abs( 100*(y_trial - y_forecasted)/y_trial )

    avg_error_meV_tr /= len(X) ;  avg_error_Perc_tr /= len(X)
    df_out = pd.concat( [y, y_frc], axis=1, sort=False )
    df_out = pd.concat([df_out, y_resid], axis=1, sort=False)
    df_out = pd.concat([df_out, X_to_save], axis=1, sort=False ) #old: pd.concat([df_out, X], axis=1, sort=False )
    if (different_regressors):
       df_aux = pd.read_csv( filepathin_training_data, header=0, engine='c', nrows=1)
       neoheader = list(df_aux); neoheader.remove(field_out_to_read); #neoheader.remove("HOMO ren (output)"); neoheader.remove("LUMO ren (output)"); neoheader.remove("Gap ren (output)")
       for label in regressor_fields:
           neoheader.remove(label)
       df_aux = pd.read_csv(filepathin_training_data, header=0, engine='c', usecols=neoheader)
       df_aux.set_index("Name", inplace=True)
       df_out = pd.concat([df_out, df_aux ], axis=1, sort=False)
       del df_aux
    filepathout_training_data = filepathin_training_data.replace( ".csv", "-with_residuals"+suffix+".csv" )    # filepathout = f"{module_io.execution_directory}residuals_from_multiple_regression_{quantity_to_forecast}.csv"
    df_out.to_csv( filepathout_training_data,index=True)


    # PART 2: CALCULATION OF RESIDUALS OF TEST AND BT DATA (using the alpha and beta found in Part 1) ------------------

    df_raw_data = pd.read_csv( filepathin_test_data, header=0, engine='c', usecols=["Name"] + regressor_fields_to_save1 + [field_out_to_read])
    y_resid = pd.DataFrame(index=df_raw_data["Name"], columns=[ field_out_to_read + "_residual"])
    df_raw_data.set_index("Name",inplace=True)
    X = df_raw_data[ regressor_fields ]
    X_to_save = df_raw_data[regressor_fields_to_save]
    y = df_raw_data[ field_out_to_read  ]
    if (verbose == 2):
        print("\n\n  TEST DATA: \n")
        print("i,  Yobser,  Yforcst,  Diff(meV),  Diff(%) ")
    avg_error_meV_tt =0.0; avg_error_Perc_tt = 0.0

    for i in range(len(X)):
        X_trial = X.iloc[i]
        y_trial = y.iloc[i]
        y_forecasted = resu.predict(np.array([ X_trial ]))[0]
        y_resid.iloc[i][field_out_to_read + "_residual"] = y_trial - y_forecasted
        if (verbose== 2):
            print(i, ",  ",y_trial, ",  ","{:.3}".format( y_forecasted), ",  ", "{:.3}".format( (y_trial - y_forecasted)), ",  ", "{:.4}".format( 100*(y_trial - y_forecasted)/y_trial ) )
        error_perc = abs(100 * (y_trial - y_forecasted) / y_trial)
        error_meV  = abs((y_trial - y_forecasted))
        avg_error_meV_tt  += error_meV;
        if (error_meV > 5.0): # This reduces the distortion; nonetheless, error_perc is not a very informative output.
           avg_error_Perc_tt  += error_perc


    avg_error_meV_tt /= len(X) ;  avg_error_Perc_tt /= len(X)
    df_out = pd.concat([X_to_save, y], axis=1, sort=False) #old: pd.concat([X, y], axis=1, sort=False)
    df_out = pd.concat([df_out, y_resid], axis=1, sort=False)
    if (different_regressors):
       df_aux = pd.read_csv( filepathin_test_data, header=0, engine='c', nrows=1)
       neoheader = list(df_aux); neoheader.remove(field_out_to_read); #neoheader.remove("HOMO ren (output)"); neoheader.remove("LUMO ren (output)"); neoheader.remove("Gap ren (output)")
       for label in regressor_fields:
           neoheader.remove(label)
       df_aux = pd.read_csv(filepathin_test_data, header=0, engine='c', usecols=neoheader)
       df_aux.set_index("Name", inplace=True)
       df_out = pd.concat([df_out, df_aux ], axis=1, sort=False)
       del df_aux
    filepathout_test_data = filepathin_test_data.replace( ".csv", "-with_residuals"+suffix+".csv" )
    df_out.to_csv( filepathout_test_data,index=True)

    if (BT_data):

       df_raw_data = pd.read_csv( filepathin_BT_data, header=0, engine='c', usecols=["Name"] + regressor_fields + [field_out_to_read])
       y_resid = pd.DataFrame(index=df_raw_data["Name"], columns=[ field_out_to_read + "_residual"])
       df_raw_data.set_index("Name",inplace=True)
       X = df_raw_data[ regressor_fields ]
       y = df_raw_data[ field_out_to_read  ]

       if (verbose == 2):
           print("\n\n  BACKTESTING DATA: \n")
           print("i,  Yobser,  Yforcst,  Diff(meV),  Diff(%) ")
       avg_error_meV_bt =0.0; avg_error_Perc_bt =0.0

       for i in range(len(X)):
           X_trial = X.iloc[i]
           y_trial = y.iloc[i]
           y_forecasted = resu.predict(np.array([ X_trial ]))[0]
           y_resid.iloc[i][field_out_to_read + "_residual"] = y_trial - y_forecasted
           if (verbose):
               print(i, ",  ",y_trial, ",  ","{:.3}".format( y_forecasted), ",  ", "{:.3}".format( (y_trial - y_forecasted)), ",  ", "{:.4}".format( 100*(y_trial - y_forecasted)/y_trial ) )
           avg_error_meV_bt  += abs((y_trial - y_forecasted)); avg_error_Perc_bt  += abs(100 * (y_trial - y_forecasted) / y_trial)

       avg_error_meV_bt /= len(X); avg_error_Perc_bt /= len(X)

       df_out = pd.concat([X, y], axis=1, sort=False)
       df_out = pd.concat([df_out, y_resid], axis=1, sort=False)
       filepathout_BT_data = filepathin_BT_data.replace( ".csv", "-with_residuals.csv" )
       df_out.to_csv( filepathout_BT_data,index=True)


    if (verbose > 0):
      print("\n\n TRAINING lin. regr.: avg_abs_error(meV)=", "{:.3}".format( avg_error_meV_tr)  , "; avg_abs_error(%)=", "{:.4}".format( avg_error_Perc_tr ))
      print(" TEST lin. regr.:     avg_abs_error(meV)=", "{:.3}".format( avg_error_meV_tt)  , "; avg_abs_error(%)=","{:.4}".format( avg_error_Perc_tt)  )
      if (BT_data): print(" BT lin. regr.:       avg_abs_error(meV)=","{:.3}".format( avg_error_meV_bt)  ,"; avg_abs_error(%)=","{:.4}".format( avg_error_Perc_bt  ))

    del df_out; del X; del y; del y_trial; del y_frc; del df_raw_data

    return

#-------------------------------------------------------------------------------------------------------------------------------------

def ml_regression( config, param_sw, filepath_training_data, filepath_validation_data ):
    '''
    :param method: (string) The keyword of the ML method to be used.
    :param quantity_to_forecast: (string) either HOMO_ren, LUMO_ren or gap_ren
    :param filepath_training_data: (string) path to the .csv file which contains the training data
    :param filepath_validation_data:     (string) path to the .csv file which contains the validation data (i.e. either "test" or "backtesting" data)
    :param li_regressors: (list of strings) list of the names of the regressors, which appear at both the filepath_training and filepath_validation
    :return:
    '''
    '''This function performs the regression of a set of data ("training") and predicts the T1_RETURN
     of a different set of data (validation) using the parameters obtained in the regression. '''

    method=config.ml_method
    quantity_to_forecast=config.quantity_to_forecast
    forecast_of_residuals = config.forecast_of_residuals
    li_regressors = config.regressor_fields_ml
    verbose = module_io.verbosity

    df_err_train_1 = module_io.errors_df.copy()
    df_err_test_1  = module_io.errors_df.copy()

    if (method == "Lasso"):
       from sklearn import linear_model
    elif ((method == "NN") or (method == "Neural networks")):
       from sklearn.neural_network import MLPRegressor
    elif ((method=="RF") or (method=="Random forest") or (method=="Random forests")):
       from sklearn.ensemble import RandomForestRegressor
    elif ((method == "KRR") or (method=="Kernel ridge regression") or (method=="Kernel ridge") ):
       from sklearn.kernel_ridge import KernelRidge
    elif ((method=="SVM") or (method == "Support vector machine")):
       from sklearn import svm
    elif ((method=="KNN") or (method == "K-nearest neighbours")):
       from sklearn import neighbors #KNeighborsRegressor
    elif (method=="GB"):
        from sklearn.ensemble import GradientBoostingRegressor
    else:
       print("\n   ERROR: Unrecognized ML method ",method,"; please enter 'Lasso', 'RF', 'NN', 'KNN,  'KRR' or 'SVM'.\n")
       exit(1)

    if (forecast_of_residuals):
       output_field = str( config.field_out_to_read+ "_residual") #str(quantity_to_forecast + "_residual")
    else:
       output_field = str(config.field_out_to_read)
    filepathout = f"{module_io.execution_directory}output_ML.csv"   
    f = open(filepathout, 'w')
    firstline = "Name;" + str(output_field) + ";"+str(output_field)+"_forecasted;Diff(meV);Diff(%)"
    print(firstline, file=f)
    f.close()

    #np.random.seed(0)  # To avoid random variation of the result on the same input dataset. Uncomment this if you want to compare results with the same ML method.

    listcols = ["Name"] + [output_field]
    if (forecast_of_residuals): listcols.append(output_field.replace("_residual", ""))
    listcols_train = listcols.copy().append(output_field.replace("_residual", "")+"_forecasted")

    df_output_tr = pd.read_csv( filepath_training_data, header=0, usecols = listcols_train )
    df_regressors_tr = pd.read_csv( filepath_training_data, header=0, usecols =  ["Name"] + li_regressors  )

    df_output_tr.set_index("Name", inplace=True)
    df_regressors_tr.set_index("Name", inplace=True)
    m = len(df_output_tr)
    Nregressors = len(li_regressors)

    X_train = np.zeros((m,Nregressors))  # We declare and initialize the independent variables (X, regressors) for the "training" (i.e. calculation of parameters of the regression).
    y_train = np.zeros(m)                # We declare and initialize the dependent variable (y) for the training.

    for i in range(m):
            y_train[i] = float( df_output_tr.iloc[i][output_field] )
            for j in range(len(li_regressors)):  # We avoid the 1st, 2nd and last fields because they are 'Name' (molecule label) and the output_variable (e.g. 'HOMO ren (output)_residual')
                X_train[i,j] = float( df_regressors_tr.iloc[i][str(li_regressors[j])] )




    #df_err_train_1["Y_full_ref"] = df_output_tr[output_field]

    # Actual ML training:

    if (method == "Lasso"):

      ml_variable = linear_model.Lasso(alpha=param_sw) # kernel{‘linear’, ‘poly’, ‘rbf’, ‘sigmoid’, ‘precomputed’},

    elif ((method == "NN") or (method == "Neural networks")):

       # ‘identity’, ‘logistic’, ‘tanh’, ‘relu’; 'identity' is like ols; Do NOT use 'relu', it is asymmetrical, it gives 0 for negative x.
       # 'identity' gives linear regression. And tanh(x) = 2·logistic(2x) - 1; Hence, in summary, there is only ONE practical activation function for NN. Y, por lo que leo en internet, NO es fácil usar otras (tendrías que escribirlas tú en una librería).
       # I don't recommend to use a solver other than "adam"; set <<hidden_layer_sizes=(3), activation='identity',>> to ckeck correctness (it must give the same as OLS)
       #Bests so far
       # Con 40, list_regressors_2:
       # (8.05 vs 8.32) ml_variable = MLPRegressor(hidden_layer_sizes=(int(param_sw),int(param_sw/2),), activation='logistic', solver='adam', max_iter=60000, alpha=0.01,learning_rate_init=0.002, momentum=0.6)
       # (7.89 vs 8.32) ml_variable = MLPRegressor(hidden_layer_sizes=( 200,100 ),                      activation='logistic',solver='adam', max_iter=400, alpha=0.01, learning_rate_init=0.002, momentum=0.6)
       # (7.82 vs 8.32) ml_variable = MLPRegressor(hidden_layer_sizes=( 200,100 ),                      activation='logistic',solver='adam', max_iter=1000, alpha=0.01, learning_rate_init=0.001, momentum=0.6)
       #  Con 40, list_regressors_3:
       # (7.69 vs 8.32)ml_variable = MLPRegressor(hidden_layer_sizes=( 300,200 ), activation='logistic',solver='adam', max_iter=1000, alpha=0.01, learning_rate_init=0.001, momentum=0.6)
       #  Con 40, list_regressors_4:
       # (7.69 vs 8.32) MLPRegressor(hidden_layer_sizes=( 300,200 ), activation='logistic',solver='adam', max_iter=1000, alpha=0.01, learning_rate_init=0.001, momentum=0.6)
       # Con list_regressors_4 y data woTIsymm: PEOR
       # IN: (400,200) mejor!!
       #   ml_variable = MLPRegressor(hidden_layer_sizes=( 400,300 ), activation='logistic',solver='adam', max_iter=10000, alpha=0.01, learning_rate_init=0.0015, momentum=0.6)
       # ml_variable = MLPRegressor(hidden_layer_sizes=(int(param_sw),int(param_sw/2),), activation='logistic', solver='adam', max_iter=10000, alpha=0.01, learning_rate_init=0.0015, momentum=0.6)
       #ml_variable = MLPRegressor(hidden_layer_sizes = (int(param_sw),1,), activation='logistic',solver='adam', max_iter=10000, alpha=0.01, learning_rate_init=0.0015, momentum=0.6)
       ml_variable = MLPRegressor(hidden_layer_sizes=(int(param_sw), int(param_sw / 2),), activation='logistic', solver='adam',max_iter=10000, alpha=0.01, learning_rate_init=0.0015, momentum=0.6)

    elif ((method == "KNN") or (method == "K-nearest neighbours")):

       ml_variable =  neighbors.KNeighborsRegressor(n_neighbors=param_sw,weights="uniform", p=1 ) # "uniform" or "distance"; (40,"uniform") slightly beats ols; leaf_size=1 #algorithm{‘auto’, ‘ball_tree’, ‘kd_tree’, ‘brute’

    elif ((method=="RF") or (method=="Random forest") or (method=="Random forests")):

       ml_variable = RandomForestRegressor( n_estimators=param_sw, min_samples_leaf=2, max_features=3 )#, max_features=2 , max_depth=1, min_samples_leaf=2, min_samples_split=4)

    elif (method == "GB"):

        ml_variable = GradientBoostingRegressor(learning_rate=0.05, n_estimators=70, subsample=1.0, criterion='friedman_mse', min_samples_split=2, min_samples_leaf=1, min_weight_fraction_leaf=0.0, max_depth=3, max_features=None, alpha=0.9 )

    elif ((method == "KRR") or (method=="Kernel ridge regression") or (method=="Kernel ridge") ):

        ml_variable = KernelRidge(alpha=param_sw, kernel='laplacian')  # kernel= "linear", "rbf", "laplacian", "polynomial", "exponential", "chi2", "sigmoid"

    elif ((method=="SVM") or (method == "Support vector machine")):

        ml_variable = svm.SVR(kernel='rbf',epsilon=0.06 ) #'poly' (e.g. with degree=2), 'rbf', 'linear', 'precomputed', 'sigmoid'
        # IMPORTANT: The regressors MUST be standardized!!
    else:

        print("\nERROR: Unknown ML method",method,"\n")
        exit(1)


    ml_variable.fit(X_train, y_train)

    df_err_train_1["Y_full_ref"] = df_output_tr[output_field.replace("_residual", "")] # Value of the ab initio renormalization (for the training dataset)
    if (forecast_of_residuals):
        df_err_train_1["Y_raw_train"] = df_output_tr[output_field.replace("_residual", "")+"_forecasted"] # ("forecasted" means "forecasted by linear ols"). This is the number given by ols linear regression (before applying ML)
    else:
        df_err_train_1["Y_raw_train"] = 0.0
    y_train_predicted = ml_variable.predict(X_train)
    df_err_train_1["Y_full_val"] = df_err_train_1["Y_raw_train"] + y_train_predicted
    df_err_train_1["y_ref"] = df_output_tr[output_field]
    df_err_train_1["y_val"] = df_err_train_1["Y_full_val"] - df_err_train_1["Y_raw_train"]
    df_err_train_1["y_err"] = df_err_train_1["y_ref"] - df_err_train_1["y_val"]

    df_output_valid = pd.DataFrame( pd.read_csv(filepath_validation_data, header=0, usecols = listcols  ) )
    df_regressors_valid = pd.read_csv( filepath_validation_data, header=0, usecols =  ["Name"] + li_regressors  )
    df_output_valid.set_index("Name", inplace=True)
    df_regressors_valid.set_index("Name", inplace=True)

    m_val = len(df_output_valid)
    X_val = np.zeros((m_val,Nregressors))  # We declare and initialize the independent variables (X, regressors) for the validation.
    y_ref = np.zeros(m_val)

    if (forecast_of_residuals): avg_abs_error_val_ols = 0
    else: avg_abs_error_val_ols = None

    for i in range(m_val):
            y_ref[i] = df_output_valid.iloc[i][output_field]  # yref are the data that we will try to predict
            if (forecast_of_residuals): avg_abs_error_val_ols += abs(y_ref[i])
            for j in range(len(li_regressors)):
                X_val[i, j] = df_regressors_valid.iloc[i][str(li_regressors[j])]

    y_val = ml_variable.predict(X_val)  # y_val is the value of "y" forecasted by the ML method (e.g. by the Neural Network).

    df_err_test_1 = df_output_valid.copy()
    df_err_test_1 = df_err_test_1.rename(columns={output_field.replace("_residual", ""): "Y_full_ref"})
    if (config.forecast_of_residuals):
        df_err_test_1 = df_err_test_1.rename(columns={output_field: "y_ref"})
        df_err_test_1["Y_raw_test"] = df_err_test_1["Y_full_ref"] - df_err_test_1["y_ref"]
    else:
        df_err_test_1["Y_raw_test"] = 0
        df_err_test_1["y_ref"] = df_err_test_1["Y_raw_test"]
    #df_err_test_1[ "Y_full_ref"] = df_output_valid[output_field.replace("_residual", "")]  # Value of the ab initio renormalization (for the test dataset)
    #df_err_test_1["Y_raw_test"] = df_output_valid[output_field.replace("_residual", "")] - df_output_valid[output_field] # This is the number given by ols linear regression (before applying ML)
    df_err_test_1["Y_full_val"] = df_err_test_1["Y_raw_test"] + y_val
    #df_err_test_1["y_ref"] = y_ref
    df_err_test_1["y_val"] = y_val
    if (config.forecast_of_residuals):
        df_err_test_1["y_err"] = df_err_test_1["y_ref"]  - df_err_test_1["y_val"]
    else:
        df_err_test_1["y_err"] = 0

    f = open(filepathout , 'a')
    avg_err_meV=0.0; avg_err_perc = 0.0
    for i in range(m_val):
            linea = df_output_valid.index[i] +  ";  " + str("{:.3}".format( y_ref[i])) + ";  " + str("{:.3}".format( y_val[i])) + ";  " + str(  "{:.3}".format( y_ref[i] - y_val[i]) )+ ";  " + str(  "{:.4}".format( 100*(y_ref[i] - y_val[i])/y_ref[i]) )
            avg_err_meV  += abs( y_ref[i] - y_val[i] )
            avg_err_perc += abs( (y_ref[i] - y_val[i])/y_ref[i]  )
            print(linea, file=f)

    avg_err_meV /= len(y_ref); avg_err_perc /= len(y_ref)

    if (forecast_of_residuals):
        avg_abs_error_val_ols = avg_abs_error_val_ols / len(X_val)
        text = "; avg. abs. error. linear regr.=" +str( "{:.3}".format(avg_abs_error_val_ols ) )
    else:
        text=""
    if (verbose > 0): print(" TEST ML ("+str(config.ml_method)+"):        avg_abs_error(meV)=","{:.3}".format( avg_err_meV),"; avg_abs_error(%)=","{:.3}".format( 100*avg_err_perc),text,"\n\n")

    f.close()

    df_err_train_1 = df_err_train_1.drop(['Y_raw_train'], axis=1)
    df_err_test_1  = df_err_test_1.drop(['Y_raw_test'], axis=1)
    df_err_train_1 = df_err_train_1[['Y_full_ref','Y_full_val','y_ref','y_val','y_err']]
    df_err_test_1  = df_err_test_1[['Y_full_ref', 'Y_full_val', 'y_ref', 'y_val', 'y_err']]

    #print("df_err_train_1",df_err_train_1);print("df_err_test_1",df_err_test_1);exit(0)

    del filepathout; del firstline; del ml_variable; del X_train; del X_val; del f; del df_regressors_tr; del df_output_tr; del df_output_valid; del df_regressors_valid
    if ( (not module_io.store_ml_errors) or (not config.forecast_of_residuals) ):
        df_err_train_1 = None; df_err_test_1 = None

    return   avg_err_meV, 100*avg_err_perc, avg_abs_error_val_ols, df_err_train_1, df_err_test_1

# -----------------------------------------------------------------------------------------------------------

def update_df_errors( iter, df_errors, df_err_1, filepath, forecast_of_residuals, all_cols, datatype  ):
    '''This function updates the dataframe which stores the errors.'''

    assert datatype in ['train','test']

    if ( (not module_io.store_ml_errors) or (not forecast_of_residuals)  ):
        del df_errors; del df_err_1; del filepath
        return None

    if (module_io.store_regressors):
        if (datatype=='train'):
            df_regressors = pd.read_csv( module_io.execution_directory + 'training_data.csv', usecols=all_cols )
        elif (datatype=='test'):
            df_regressors = pd.read_csv( module_io.execution_directory + 'test_data.csv', usecols=all_cols  )
        df_regressors.set_index("Name", inplace=True)
        df_aux = module_io.errors_df.copy()#pd.DataFrame( columns= list(module_io.errors_df) + list(df_regressors), index=df_regressors.index )
    else:
        df_aux =  module_io.errors_df.copy()

    df_aux["Y_full_ref"] = df_err_1["Y_full_ref"]
    df_aux["y_ref"] = df_err_1["y_ref"]
    df_aux["Y_full_val"] = df_err_1["Y_full_val"]
    df_aux["y_val"] = df_err_1["y_val"]
    df_aux["y_err"] = df_err_1["y_ref"] - df_err_1["y_val"]

    if (module_io.store_regressors):
        df_aux = pd.concat( [df_aux,df_regressors],axis=1)

    df_errors = pd.concat( [df_errors,df_aux], axis=0 )

    if (iter==100):
        df_errors.to_csv(filepath, index=False)
        df_errors = module_io.errors_df.copy()
    elif ((iter % 100) == 0):
        df_errors.to_csv(filepath, mode='a', header=False, index=False)
        df_errors = module_io.errors_df.copy()

    del df_aux

    return df_errors

#-------------------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------------------
















'''
RAVENPACK FUNCTION:



def neural_networks_regression ( name_counter, vec_names, path_in_training, path_in_validation, path_out ):

#This function performs the regression of a set of data ("training") and predicts the T1_RETURN
# of a different set of data (validation) using the parameters obtained in the regression. 

 from sklearn.neural_network import MLPRegressor


 print(" Now performing regressions.")

 filepathout = path_out + "/AllStocks-results-regression.csv"
 f = open( filepathout , 'w')
 firstline = "DATE" + "," + "RP_ENTITY_ID" + "," + "PREDICTION" + "," + "T1_RETURN"
 print(firstline, file=f) 
 f.close()
 
 
 for nc in range(name_counter):
 
   np.random.seed(0)  # To avoid random variation of the result on the same input dataset.
  
   file_name = path_in_training + "/" + str(vec_names[nc]).rstrip() + ".csv"
   print(" Now reading",file_name)
   datos = pd.read_csv( file_name,header=0, dtype={'RP_ENTITY_ID': object}  ) 
   m = datos.DATE.shape[0] 
   list_fields = list(datos)



   Nregressors = len( list_fields ) -3

   X_train = np.zeros((m, Nregressors))    # We declare and initialize the independent variables (X, regressors) for the "training" (i.e. calculation of parameters of the regression).
   y_train = np.zeros(m)                   # We declare and initialize the dependent variable (y) for the training.

   for i in range(m):
        y_train[i] = datos.iloc[i]["T1_RETURN"]
        for j in range(2,len(list_fields)-1):   # We avoid the 1st, 2nd and last fields because they are DATE, ID and T1_RETURN.
          if not isnan(datos.iloc[i][ str(list_fields[j]) ]):
             X_train[i,j-2] = datos.iloc[i][ str(list_fields[j]) ]
      
  

   ml_variable = MLPRegressor( hidden_layer_sizes=(100,300,100 ), activation='relu',solver='adam',max_iter=1000,learning_rate_init=0.01,alpha=0.01)
   ml_variable.fit(X_train, y_train)     


   file_name = path_in_validation + "/" + str(vec_names[nc]).rstrip() + ".csv"
   datos_valid = pd.read_csv( file_name,header=0, dtype={'RP_ENTITY_ID': object}  ) 
   m_val = datos_valid.DATE.shape[0] 

   X_val = np.zeros((m_val, Nregressors))    # We declare and initialize the independent variables (X, regressors) for the validation.
   y_val = np.zeros(m_val)
   y_ref = np.zeros(m_val)

   for i in range(m_val):
        y_ref[i] = datos_valid.iloc[i]["T1_RETURN"]   # yref are the data that we will try to predict
        for j in range(2,len(list_fields)-1):   
          if not isnan(datos_valid.iloc[i][ str(list_fields[j]) ]):
             X_val[i,j-2] = datos_valid.iloc[i][ str(list_fields[j]) ]

   
   y_val = ml_variable.predict( X_val ) 
   
   f = open(path_out+"/AllStocks-results-regression.csv" , 'a')           
         
   for i in range(m_val):   
           
         linea = datos_valid.iloc[i][ 'DATE' ] + "," + str(datos_valid.iloc[i][ 'RP_ENTITY_ID' ]) + "," 
         linea = linea + str( y_val[i] ) + "," + str( y_ref[i] )
         print(linea, file=f) 
   
   f.close() 


 # We sort the file of results by date:
 datos = pd.read_csv( filepathout, header=0, dtype={'RP_ENTITY_ID': object}  ) 
 datos = datos.sort_values(by=["DATE","RP_ENTITY_ID"]).copy()
 datos.to_csv( filepathout ,index=False)

 return

#-----------------------------------------------------------------------------------------------------------

'''

#-----------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------

def create_array_random_variables():
    '''This function simply creates a file whose columns are standardized random variables with uniform distribution.'''

    from numpy.random import uniform

    df0 = pd.DataFrame(columns=["rv1","rv2","rv3","rv4","rv5","rv6","rv7","rv8","rv9","rv10"])
    for i in range(10):
      vec = uniform(-0.5,0.5,148)
      media = np.mean(vec)
      stdev = np.std(vec)
      vec = (vec-media)/stdev
      df0[ "rv"+str(i+1) ] = vec

    df0.to_csv( f"{module_io.execution_directory}random_vectors.csv",index=False)

    return

# ------------------------------------------------------------------------------------------------------------------------------------


'''
    if (quantity_to_forecast == "HOMO"):
       field_out_to_read =  "HOMO ren (output)"
       regressor_fields  = ["b","Gap_PBE", "deg HOMO", "avgoccupied(eV-1)",
                           "HOMO-(HOMO-1)", "HOMO-(HOMO-2)", "HOMO-(HOMO-3)", "HOMO-(HOMO-4)","HOMO-(HOMO-5)",
                           "inv(HOMO-(HOMO-1))", "inv(HOMO-(HOMO-2))", "inv(HOMO-(HOMO-3))", "inv(HOMO-(HOMO-4))", "inv(HOMO-(HOMO-5))",
                           "inv((LUMO+1)-HOMO)","inv((LUMO+2)-HOMO)"   ]

    elif (quantity_to_forecast == "LUMO"):
       field_out_to_read =  "LUMO ren (output)"
       regressor_fields  = [ "d","Nat", "Gap_PBE", "deg LUMO",  "avgempty(eV-1)",
                           "(LUMO+1)-LUMO", "(LUMO+2)-LUMO", "(LUMO+3)-LUMO", "(LUMO+4)-LUMO", "(LUMO+5)-LUMO",
                           "inv((LUMO+1)-LUMO)","inv((LUMO+2)-LUMO)","inv((LUMO+3)-LUMO)","inv((LUMO+4)-LUMO)","inv((LUMO+5)-LUMO)",
                             "inv(LUMO-(HOMO-1))","inv(LUMO-(HOMO-2))" ]

    elif (quantity_to_forecast == "gap"):
       field_out_to_read =  "Gap ren (output)"
       regressor_fields  = ["b","d","Nat", "Gap_PBE", "deg HOMO","deg LUMO","avgoccupied(eV-1)", "avgempty(eV-1)",
                           "HOMO-(HOMO-1)", "HOMO-(HOMO-2)", "HOMO-(HOMO-3)",
                           "inv(HOMO-(HOMO-1))", "inv(HOMO-(HOMO-2))", "inv(HOMO-(HOMO-3))",
                            "(LUMO+1)-LUMO", "(LUMO+2)-LUMO", "(LUMO+3)-LUMO",
                            "inv((LUMO+1)-LUMO)", "inv((LUMO+2)-LUMO)", "inv((LUMO+3)-LUMO)" ]
    else:
       print("   ERROR: Unrecognized quantity to forecast: "+module_io.quantity_to_forecast+"; Please, check the module_io.py file.\n"); exit(1)


'''