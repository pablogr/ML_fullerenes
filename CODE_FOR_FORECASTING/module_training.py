
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

def ols_multiregression( quantity, field_out_to_read, regressor_fields, regressor_fields_to_save, filepathin_training_data, filepathout_training_data, verbose, different_regressors ):
    ''' This function finds the residuals of the multiple regression (through Ordinary Least Squares) and writes it to
    files together with the regressors for all three datasets (training, test, backtesting).

    :param field_out_to_read: (string) must be either "HOMO", "LUMO" or "gap".
    :param regressor_fields: Names (labels) of the regressors to be used in the OLS (linear) regression
    :param regressor_fields_to_save: Names (labels) of the regressors to be stored in the output file
    :param filepathin_training_data: Path to the file containing training data. This file, as well as the two files below, must
           contain at least the columns which correspond to the regressors, the name of the molecule and the output value (HOMO_ren, LUMO_ren or GAP_ren).
    :param filepathin_test_data:     Path to the file containing the regressors of the molecule to forecast (just one row after the header).
    :param verbose:  (int) True if enhanced information is to be printed to screen.
    :param different_regressors: (Boolean) False if the regressors which are used in the regression and in the ML part are the same; true otherwise.
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
    #filepathout_training_data = filepathout_training_data.replace( ".csv", "-with_residuals.csv" )    # filepathout = f"{module_io.execution_directory}residuals_from_multiple_regression_{quantity_to_forecast}.csv"
    df_out.to_csv( filepathout_training_data,index=True)


    # PART 2: FORECAST BASED ON LINEAR REGRESSION (using the alpha and beta found in Part 1) ------------------

    df_raw_data = pd.read_csv( module_io.filepath_regressors_molecule_to_forecast, header=0, engine='c', usecols=["Name"] + regressor_fields_to_save1  )
    df_raw_data.set_index("Name",inplace=True)
    X = df_raw_data[ regressor_fields ]

    y_forecasted = resu.predict(np.array([ X.iloc[0] ]))[0]
    print("   * The "+quantity+" renormalization forecasted through linear regression is","{:.6}".format( y_forecasted))

    del df_out; del X; del y; del y_trial; del y_frc; del df_raw_data

    return y_forecasted

#-------------------------------------------------------------------------------------------------------------------------------------

def ml_regression( config, param_sw, filepath_training_data, filepath_validation_data, my_forecast_lr ):
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
    forecast_of_residuals = config.forecast_of_residuals
    li_regressors = config.regressor_fields_ml

    df_err_train_1 = module_io.errors_df.copy()

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
       from sklearn import  neighbors #KNeighborsRegressor
    else:
       print("\n   ERROR: Unrecognized ML method ",method,"; please enter 'Lasso', 'RF', 'NN', 'KNN,  'KRR' or 'SVM'.\n")
       exit(1)

    if (forecast_of_residuals):
       output_field = str( config.field_out_to_read+ "_residual") #str(quantity_to_forecast + "_residual")
    else:
       output_field = str(config.field_out_to_read)

    #np.random.seed(0)  # To avoid random variation of the result on the same input dataset. Uncomment this if you want to compare results with the same ML method.

    listcols = ["Name"] + [output_field]
    if (forecast_of_residuals): listcols.append(output_field.replace("_residual", ""))
    listcols_train = listcols.copy().append(output_field.replace("_residual", "")+"_forecasted")
    listcols_forecast = listcols.remove(output_field)

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

    elif ((method == "KRR") or (method=="Kernel ridge regression") or (method=="Kernel ridge") ):

        ml_variable = KernelRidge(alpha=param_sw, kernel='laplacian')  # kernel= "linear", "rbf", "laplacian", "polynomial", "exponential", "chi2", "sigmoid"

    elif ((method=="SVM") or (method == "Support vector machine")):

        ml_variable = svm.SVR(kernel='linear')

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

    df_output_valid = pd.DataFrame( pd.read_csv(filepath_validation_data, header=0, usecols = listcols_forecast  ) )
    df_regressors_valid = pd.read_csv( filepath_validation_data, header=0, usecols =  ["Name"] + li_regressors  )
    df_output_valid.set_index("Name", inplace=True)
    df_regressors_valid.set_index("Name", inplace=True)

    m_val = len(df_output_valid)
    X_val = np.zeros((m_val,Nregressors))  # We declare and initialize the independent variables (X, regressors) for the validation.
    y_ref = np.zeros(m_val)

    for i in range(m_val):
            #y_ref[i] = df_output_valid.iloc[i][output_field]  # yref are the data that we will try to predict
            #if (forecast_of_residuals): avg_abs_error_val_ols += abs(y_ref[i])
            for j in range(len(li_regressors)):
                X_val[i, j] = df_regressors_valid.iloc[i][str(li_regressors[j])]

    y_val = ml_variable.predict(X_val)  # y_val is the value of "y" forecasted by the ML method (e.g. by the Neural Network).


    del ml_variable; del X_train; del X_val; del df_regressors_tr; del df_output_tr; del df_output_valid; del df_regressors_valid

    return y_val[0]

# -----------------------------------------------------------------------------------------------------------



