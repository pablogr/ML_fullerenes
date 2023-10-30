import numpy as np

import module_io
import pandas as pd
from os import makedirs, path, remove

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

def initialization(  config, quantity  ):
    '''
    :return: field_out_to_read (string): Either "HOMO ren (output)", "LUMO ren (output)" or  "Gap ren (output)"
    :return: regressor_fields (list of strings): The list of regressors to use in each case.
    '''

    for filepath in [module_io.execution_directory,module_io.execution_directory+"/Results",module_io.execution_directory+"/Plots"]:
        if not path.exists(filepath): makedirs(filepath)

    print("\n   PROGRAM FOR FORECASTING OF FULLERENES FROZEN-PHONON "+ config.quantity_to_forecast+" RENORMALIZATION USING MACHINE LEARNING \n")

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
    if config.forecast_of_residuals:
        print("   The program will forecast the residuals of the OLS regression.")
    else:
        print("   The program will forecast the output itself (not its residuals of an OLS regression).")
    if not (module_io.BT_data):
        print("   These calculations include no backtesting data; just training and test sets.")
    else:
        print("   These calculations include training, test and backtesting data.")

    print("   The input and output data are read from",module_io.filepath_raw_input_data)
    df_in = pd.read_csv(module_io.filepath_raw_input_data, header=0, engine='c')

    if module_io.standardize_regressors: df_in = standardize_regressors(df_in)

    N_data_test = round(len(df_in) * module_io.ratio_data_test)

    if not (module_io.BT_data):
        N_data_training = len(df_in) - N_data_test
    else:
        N_data_training = len(df_in) - 2 * N_data_test

    print("\n   Number of data for training:    ",N_data_training,"(","{:.3}".format( 100*N_data_training/(len(df_in) )),"%)")
    print("   Number of data for test:        ",N_data_test," (","{:.3}".format( 100*N_data_test/(len(df_in) )),"%)\n\n")


    #if (module_io.verbosity==0): print("\n           Avgd.     Avgd. \n   param, |error|,  |error| \n           (meV)     (%)")

    if ((module_io.ratio_data_test < 0.01) or (module_io.ratio_data_test > 0.5)):
        print("\n    ERROR: The ratio of the input data that was chosen for the test is too low (",
              module_io.ratio_data_test, "). Please, modify it at module_io.py.\n");
        exit(1)

    if (module_io.BT_data): print("   Number of data for backtesting: ",len(df_in) -N_data_training-N_data_test," (","{:.3}".format( 100*(len(df_in) -N_data_training-N_data_test)/(len(df_in) )),"%)")

    if not module_io.use_ensemble_for_ols:
        all_cols = [ "Name", quantity.replace("gap","Gap") +" ren (output)"] + list(set(config.regressor_fields_regression) | set(config.regressor_fields_ml))
    else:
        all_cols = ["Name", quantity.replace("gap", "Gap") + " ren (output)"] + list(set(config.regressor_fields_regression[0]+config.regressor_fields_regression[1]) | set(config.regressor_fields_ml))

    return  df_in, config.field_out_to_read, config.regressor_fields_regression,  config.regressor_fields_ml, all_cols, config.different_regressors, config.ml_method, module_io.N_essais, module_io.BT_data, N_data_training, N_data_test,  module_io.verbosity

#--------------------------------------------------------------------------------------------------------------------------------------

def middle_message(config, param_sw, sum_error_meV, sum_error_perc, sum_avg_abs_error_val_ols):

    print("\n\n\n======= "+config.quantity_to_forecast+" renormalization using "+config.ml_method+" ==============================\n",module_io.N_essais," random trials.")

    if not ( config.forecast_of_residuals):
        print("\n           Avgd.     Avgd. \n   param, |error|,  |error| \n           (meV)     (%)")
        print("   ", param_sw, ", ", "{:.5}".format(sum_error_meV), ", ", "{:.5}".format(sum_error_perc))
    else: # Forecast of residuals
        print("\n           Avgd.     Avgd.    Avgd.\n   param, |error|,  |error|, |error| \n          ML (meV)   ML (%)   ols linear regr (meV)")
        print("   ", param_sw, ", ", "{:.5}".format(sum_error_meV), ", ", "{:.5}".format(sum_error_perc), ", ", "{:.5}".format(sum_avg_abs_error_val_ols))

#--------------------------------------------------------------------------------------------------------------------------------------

def final_message(config):

    print("\n   Note: The Averaged absolute value of the error in percentage is just an orientation. The main output is << Avgd. |error| (meV) >>. ")
    if (config.forecast_of_residuals):
        print("\n   **** The calculations of "+config.quantity_to_forecast+" using "+config.ml_method+" and using the residuals from linear regression concluded. **********************************************************************************************************\n")
    else:
        print("\n   **** The calculations of " + config.quantity_to_forecast + " using " + config.ml_method + " (not using the residuals of linear regression) concluded. **********************************************************************************************************\n")

#--------------------------------------------------------------------------------------------------------------------------------------

#Functions for other purpose, data preparation:

def extract_data_from_fullerenes_2021( pathin ):

        #pathnames = r"/Users/pgr/Desktop/Ciencia/HamburgProjects/_Fullerenes-1-puros/extract_info_fullerenes_plus_calcs_dic2019/EXTRACT_DATA_FULLERENES_2021/calculos/"
        #pathbase = r"/Users/pgr/Desktop/Ciencia/HamburgProjects/_Fullerenes-1-puros/extract_info_fullerenes_plus_calcs_dic2019/EXTRACT_DATA_FULLERENES_2021/base/"
        pathnames = pathin + "calculos/"
        pathbase =  pathin + "base/"
        df0 = pd.read_csv( f"{pathnames}lista_fulerenos.csv",header=0)
        list_fullerenes = df0["Name"].values.tolist()

        #list_fullerenes.remove("C36-D2d-14")

        for Nmaxunoccavg in [20]:
          print(  "   --> Average for unoccupied with up to ",Nmaxunoccavg ," states.")
          df_out = pd.DataFrame()
          for i in range(len(list_fullerenes)):
              label_fullerene = list_fullerenes[i]
              print("   Now analysing", label_fullerene)
              foldername = pathnames + label_fullerene + "/"
              HOMO, LUMO, avgoccupied, avgempty, countocc, countempty, difHOMO1, difHOMO2, difHOMO3, difHOMO4, difHOMO5, difHOMO6, difLUMO1, difLUMO2, difLUMO3, difLUMO4, difLUMO5, difLUMO6 = find_avgseigenvalues( foldername, Nmaxunoccavg )
              #print("   Now analysing", label_fullerene," ; avg. occupied:", avgoccupied )
              df_onerow = pd.DataFrame( data={ "molecule":[label_fullerene],"HOMO(eV)":[HOMO], "LUMO(eV)":[LUMO], "avgoccupied(eV)":[avgoccupied], "avgempty(eV)":[avgempty], "countocc":[countocc], "countempty":[countempty],
                                        "HOMO-(HOMO-1)":[difHOMO1], "HOMO-(HOMO-2)":[difHOMO2], "HOMO-(HOMO-3)":[difHOMO3], "HOMO-(HOMO-4)":[difHOMO4], "HOMO-(HOMO-5)":[difHOMO5], "HOMO-(HOMO-6)":[difHOMO6],
                                        "(LUMO+1)-LUMO":[difLUMO1],"(LUMO+2)-LUMO":[difLUMO2], "(LUMO+3)-LUMO":[difLUMO3], "(LUMO+4)-LUMO":[difLUMO4], "(LUMO+5)-LUMO":[difLUMO5], "(LUMO+6)-LUMO":[difLUMO6] } )
              df_out = pd.concat( [df_out, df_onerow], axis=0, sort=False)

          df_out.to_csv( f"{pathin}ExtractInfoFullerenes2021_final_results-unocc{str(Nmaxunoccavg)}.csv",index=False )

        print("\n   ** The calculations finished satisfactorily. ** \n\n")

        exit()

#--------------------------------------------------------------------------------

def copy_files( list_fullerenes, pathnames, pathbase ):

        from shutil import copyfile

        # Example:
        #pathnames = r"/Users/pgr/Desktop/Ciencia/HamburgProjects/_Fullerenes-1-puros/extract_info_fullerenes_plus_calcs_dic2019/EXTRACT_DATA_FULLERENES_2021/calculos/"
        #pathbase = r"/Users/pgr/Desktop/Ciencia/HamburgProjects/_Fullerenes-1-puros/extract_info_fullerenes_plus_calcs_dic2019/EXTRACT_DATA_FULLERENES_2021/base/"
        #df0 = pd.read_csv(f"{pathnames}lista_fulerenos.csv", header=0)
        #list_fullerenes = df0["Name"].values.tolist()

        # Copying files and writing number of carbons:
        for i in range(len(list_fullerenes)):
           label_fullerene = list_fullerenes[i]
           Ncarbons = label_fullerene.partition("-")[0]
           Ncarbons = Ncarbons.replace("C","")
           stri = "Number_of_atoms = " + str(Ncarbons)
           print(i,label_fullerene, stri )
           namefolder = pathnames + label_fullerene + "/"
           if not (path.exists(namefolder)):
              makedirs( namefolder )
           copyfile( pathbase+"compile_and_run.sh",namefolder+"compile_and_run.sh" )
           copyfile(pathbase + "extract_data_fullerenes.f90", namefolder+ "extract_data_fullerenes.f90")
           copyfile(pathbase + "extract_data_fullerenes.in", namefolder+ "extract_data_fullerenes.in")
           file_object = open(namefolder+ "extract_data_fullerenes.in", 'a')
           file_object.write( stri )
           file_object.close()

#----------------------------------------------------------------------------------------------

def find_avgseigenvalues( foldername, Nmaxunoccavg ):

         #with open( pathnames+"C38-C2v-15/eigenval.xml" ) as infile, open( pathnames+"C38-C2v-15/fileaux1.csv", 'w') as outfile:
         with open( foldername + "eigenval.xml") as infile, open( foldername + "fileaux1.csv",'w') as outfile:
            outfile.write("eigenvalues\n")
            copy = False
            for line in infile:
               if ( ("<EIGENVALUES type=" in line.strip() ) ) :
                  copy = True
                  continue
               elif "</EIGENVALUES>" in line.strip():
                  copy = False
                  break
               else:
                  if ( (copy) and (line.strip() != "")):
                     outfile.write(line)

         with open( foldername +"eigenval.xml" ) as infile, open( foldername +"fileaux2.csv", 'w') as outfile:
            outfile.write("occupations\n")
            copy = False
            for line in infile:
               if ( ("<OCCUPATIONS type=" in line.strip() ) ) :
                  copy = True
                  continue
               elif "</OCCUPATIONS>" in line.strip():
                  break
               else:
                  if ( (copy) and (line.strip() != "")):
                     outfile.write(line)

         df1 = pd.read_csv(foldername +"fileaux1.csv")
         df2 = pd.read_csv(foldername +"fileaux2.csv")
         df3 = pd.concat( [df1,df2], axis=1, sort=False  )

         hartree_to_eV = 27.211385056

         avgoccupied = 0; countocc = 0
         avgempty = 0;    countempty = 0
         iHOMO = -1

         # Find HOMO and LUMO
         HOMO = -11111111; LUMO = -11111112;
         for i in range(len(df3)):
            if ( (df3.iloc[i]["occupations"] == 1) and (df3.iloc[i+1]["occupations"] == 0) ):
               HOMO = df3.iloc[i]["eigenvalues"]
               LUMO = df3.iloc[i+1]["eigenvalues"]
               iHOMO = i
               break
         if ( (HOMO < -11111110 ) or (LUMO < -11111110 ) ) :
            print("  ERROR: HOMO and LUMO not properly read.",foldername,"\n"); exit(1)

         difHOMO1 = ( HOMO - df3.iloc[i-1]["eigenvalues"] )*hartree_to_eV
         difHOMO2 = ( HOMO - df3.iloc[i - 2]["eigenvalues"] )*hartree_to_eV
         difHOMO3 = ( HOMO - df3.iloc[i - 3]["eigenvalues"] )*hartree_to_eV
         difHOMO4 = ( HOMO - df3.iloc[i - 4]["eigenvalues"] )*hartree_to_eV
         difHOMO5 = ( HOMO - df3.iloc[i - 5]["eigenvalues"] )*hartree_to_eV
         difHOMO6 = ( HOMO - df3.iloc[i - 6]["eigenvalues"] )*hartree_to_eV
         difLUMO1 = ( df3.iloc[i + 2]["eigenvalues"] - LUMO )*hartree_to_eV
         difLUMO2 = ( df3.iloc[i + 3]["eigenvalues"] - LUMO )*hartree_to_eV
         difLUMO3 = ( df3.iloc[i + 4]["eigenvalues"] - LUMO )*hartree_to_eV
         difLUMO4 = ( df3.iloc[i + 5]["eigenvalues"] - LUMO )*hartree_to_eV
         difLUMO5 = ( df3.iloc[i + 6]["eigenvalues"] - LUMO )*hartree_to_eV
         difLUMO6 = ( df3.iloc[i + 7]["eigenvalues"] - LUMO )*hartree_to_eV

         # Find averages of occupied and unoccupied states; the 0 in energy is the HOMO.
         for i in range( len(df3)):
            if df3.iloc[i]["occupations"] == 1:
               countocc +=1
               # avgoccupied += ( df3.iloc[i]["eigenvalues"] - HOMO ) #xx
               di = ( df3.iloc[i]["eigenvalues"] - HOMO )
               if ( abs(di) > 0.02/hartree_to_eV ):
                  avgoccupied += 1/di
            elif df3.iloc[i]["occupations"] == 0:
               countempty += 1
               if (countempty > Nmaxunoccavg ):  # if (countempty > countocc / 2):
                   countempty -= 1
                   break
               #avgempty += (df3.iloc[i]["eigenvalues"] - LUMO) #xx   ;   old: avgempty += ( df3.iloc[i]["eigenvalues"] - HOMO )
               di = (df3.iloc[i]["eigenvalues"] - LUMO)
               if (abs(di) > 0.02 / hartree_to_eV):
                   avgempty += 1/di
               #print(  hartree_to_eV*(df3.iloc[i]["eigenvalues"] - LUMO) )
            else:
               print("  ERROR: Check your data. \n"); exit(0)
         #avgoccupied = hartree_to_eV * avgoccupied / countocc
         #avgempty    = hartree_to_eV *  avgempty / countempty #avgempty = hartree_to_eV * avgempty / (countempty-1)
         avgoccupied = avgoccupied / ( hartree_to_eV * countocc)
         avgempty    = avgempty / (hartree_to_eV * countempty) #avgempty = hartree_to_eV * avgempty / (countempty-1)
         HOMO /= hartree_to_eV  #xx
         LUMO /= hartree_to_eV  #xx
         #print("HOMO",HOMO); print("LUMO",LUMO); print(" Average occupied:", avgoccupied); print(" Average empty:", avgempty )

         return HOMO, LUMO, avgoccupied, avgempty, countocc, countempty, difHOMO1, difHOMO2, difHOMO3, difHOMO4, difHOMO5, difHOMO6, difLUMO1, difLUMO2, difLUMO3, difLUMO4, difLUMO5, difLUMO6

#------------------------------------------------------------------------------------------

def initialize_error_files(my_configuration, method, quantity, forecast_of_residuals, param_sw, all_cols ):
    '''This function initializes the files where the errors will be stored (to do a histogram later). It also deletes
    existing versions of them.'''

    if (forecast_of_residuals): suffix = "_fromresid"
    else: suffix = "_notfromresid"
    filepath_err_train = my_configuration.filepath_errors_train + "_" + method+"_"+quantity+"_"+suffix+"_"+str(param_sw)
    filepath_err_test  = my_configuration.filepath_errors_test  + "_" + method + "_" + quantity + "_" + suffix + "_"+ str(param_sw)

    if not (module_io.alternative_regressors):
        filepath_err_train +=  ".csv"
        filepath_err_test +=  ".csv"
    else:
        suffix = ""
        if (quantity=="HOMO"):
           regressor_list =  module_io.list_regressors_HOMOr.copy()
        elif (quantity=="LUMO"):
           regressor_list =  module_io.list_regressors_LUMOr.copy()
        elif (quantity=="gap"):
           regressor_list =  module_io.list_regressors_gapr.copy()
        else:
            raise Exception('\nERROR: Unrecognized quantity '+quantity+"\n")
        regressor = regressor_list[-1]
        regressor = str(regressor)
        for spec in ['"','/','\\', ' ']:
            regressor = regressor.replace(spec,"")
        suffix += "_" + regressor
        filepath_err_train += suffix + ".csv"
        filepath_err_test  += suffix + ".csv"

    for filepath in [filepath_err_train, filepath_err_test]:
        open(filepath, 'w').close()
    '''
    if not (module_io.store_regressors):

        for filepath in [ filepath_err_train, filepath_err_test ]:
            if (path.exists(filepath)):
                remove(filepath)
            with open(filepath, 'w') as f:
                f.write('Y_full_ref,Y_full_val,y_ref,y_val,y_err\n')
            f.close()

    else:

        list_cols0 = ['Y_full_ref', 'Y_full_val', 'y_ref', 'y_val', 'y_err']
        for filepath,filedata in zip([filepath_err_train, filepath_err_test],['training_data-with_residuals.csv','test_data.csv']):
            list_cols = list_cols0 + all_cols #list(pd.read_csv(module_io.execution_directory + filedata))
            list_cols.remove("Name")
            firstrow = ""
            for mycol in list_cols: firstrow += str(mycol) + ","
            firstrow = firstrow[:len(firstrow) - 1]
            firstrow += "\n"
            if (path.exists(filepath)):
                remove(filepath)
            with open(filepath, 'w') as f:
                f.write(firstrow)
            f.close()
    '''

    #f_train = open( filepath_err_train, 'a')
    #f_test  = open( filepath_err_test, 'a')
    #  file_object.write('hello')
    # file_object.close()

    del my_configuration

    return filepath_err_train, filepath_err_test

#------------------------------------------------------------------------------------------
