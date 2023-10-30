'''
This is the module which contains the data of input parameters.
'''

from pandas import DataFrame
from os import getcwd, makedirs, path

ratio_data_test = 0.2          # This is the ratio of the raw data which is used for test and for backtesting. For the training, the ratio of data used is 1 - 2·ratio_data_test_and_BT.
N_essais = 1000                # (e.g. 1000) This is the number of times that one runs the ML algorithm, this is the number of times that one randomly selects the training, test and BT datasets and runs ML on them, collects results, etc.
verbosity = 0                  # 0, 1 or 2, for different length of the outputs to screen.
BT_data = False                # True if the whole dataset must be split into 3 parts (training, test, backtesting) and False if BT is not necessary.
standardize_regressors = False # If True then the regressors are standardized at the beginning of the algorithm (before splitting the whole dataset into training and test datasets).
store_ml_errors        = True  # If true then the in-sample and out-of-sample errors are stored to files.
store_regressors       = True  # If true and store_ml_errors also true, then the regressors are stored together with the outputs in Executions/Results.
alternative_regressors = False  # Set to True to try non-optimized lists of regressors (see below).
use_ensemble_for_ols   = True   # If True, then the linear regression (LR) is calculated using two different lists of regressors; the forecast is the average of both outputs. Set to True if you want to use Eva's regressors and the quantity is either LUMO_ren or Gap_ren

code_dir = getcwd()
execution_directory = code_dir+"/Execution/"
if not (path.exists(execution_directory)): makedirs(execution_directory)

if (not alternative_regressors):
    filepath_raw_input_data = f"Data_ab_initio/fullerenes_data_2022_wo3outliers_enhanced.csv"
else: # alternative regressors
    filepath_raw_input_data = f"Data_ab_initio/fullerenes_data_2023_wo3outliers_more_enhanced.csv"

if (not alternative_regressors):

    ### LIST OF EVA'S ALTERNATIVE REGRESSORS (SECTION OF II.B.3 OF THE ARTICLE) ############################################

    list_regressors_lin_ols_HOMO = ["avgoccupied(eV-1)", "HOMO-(HOMO-1)", "HOMO-(HOMO-2)", "inv(HOMO-(HOMO-1))","inv(HOMO-(HOMO-2))"]

    list_regressors_lin_ols_LUMO_1 = ["Gap_PBE", "avgempty(eV-1)", "HOMO-(HOMO-1)", "HOMO-(HOMO-2)"]
    list_regressors_lin_ols_LUMO_2 = ["Gap_PBE", "inv((LUMO+1)-LUMO)", "inv((LUMO+2)-LUMO)"]
    list_regressors_lin_ols_LUMO = [list_regressors_lin_ols_LUMO_1, list_regressors_lin_ols_LUMO_2]

    list_regressors_lin_ols_gap_1 = ["Gap_PBE", "avgoccupied(eV-1)", "inv((LUMO+1)-LUMO)"]
    list_regressors_lin_ols_gap_2 = ["Gap_B3LYP", "avgoccupied(eV-1)", "inv((LUMO+1)-LUMO)"]
    list_regressors_lin_ols_gapr = [list_regressors_lin_ols_gap_1, list_regressors_lin_ols_gap_2]

    # end of LIST OF EVA'S ALTERNATIVE REGRESSORS ############################################################################


    list_regressors_lin_ols_HOMO = ["Gap_PBE", "avgoccupied(eV-1)",
                     "HOMO-(HOMO-1)", "HOMO-(HOMO-2)", "HOMO-(HOMO-3)", "HOMO-(HOMO-4)", "HOMO-(HOMO-5)",
                     "inv(HOMO-(HOMO-1))", "inv(HOMO-(HOMO-2))", "inv(HOMO-(HOMO-3))","inv(HOMO-(HOMO-4))", "inv(HOMO-(HOMO-5))",
                     "inv(LUMO-(HOMO-1))", "inv(LUMO-(HOMO-2))", "inv((LUMO+1)-HOMO)","inv((LUMO+2)-HOMO)"]
    list_regressors_lin_ols_LUMO = ["Gap_PBE", "avgempty(eV-1)", "HOMO-(HOMO-1)", "inv((LUMO+2)-LUMO)", "inv((LUMO+1)-LUMO)"]
    list_regressors_lin_ols_gapr = ["Gap_PBE", "avgoccupied(eV-1)", "avgempty(eV-1)", "inv(HOMO-(HOMO-1))", "inv((LUMO+1)-LUMO)"]
    
    ''' 
    list_regressors_RF_HOMOr = list_regressors_lin_ols_HOMO.copy() #["Gap_PBE","Gap_B3LYP","avgoccupied(eV-1)" ,  "inv(HOMO-(HOMO-1))", "inv(HOMO-(HOMO-2))", "inv(HOMO-(HOMO-3))","inv(HOMO-(HOMO-4))" , "HOMO-(HOMO-1)" ,  ]
    list_regressors_RF_LUMOr = ["Gap_PBE","Gap_B3LYP","avgempty(eV-1)" ,"inv((LUMO+4)-LUMO)",  "inv((LUMO+3)-LUMO)", "inv((LUMO+2)-LUMO)","inv((LUMO+1)-LUMO)"  ,   ]
    list_regressors_RF_gapr  = ["Gap_PBE","Gap_B3LYP","avgempty(eV-1)" ,"avgoccupied(eV-1)" ,  "inv(HOMO-(HOMO-1))", "inv(HOMO-(HOMO-2))", "inv(HOMO-(HOMO-3))","inv(HOMO-(HOMO-4))", "inv((LUMO+4)-LUMO)",  "inv((LUMO+3)-LUMO)", "inv((LUMO+2)-LUMO)","inv((LUMO+1)-LUMO)"  ,   ]

    list_regressors_NN_HOMOr = ["Gap_B3LYP","Gap_PBE","avgoccupied(eV-1)","avgempty(eV-1)","deg HOMO","deg LUMO","ElecEigvalMean","ElecEigvalVariance","ElecEigvalSkewness","ElecEigvalKurtosis","HOMO-(HOMO-1)","HOMO-(HOMO-2)","HOMO-(HOMO-3)","HOMO-(HOMO-4)","HOMO-(HOMO-5)","HOMO-(HOMO-6)","(LUMO+1)-LUMO","(LUMO+2)-LUMO","(LUMO+3)-LUMO","(LUMO+4)-LUMO","(LUMO+5)-LUMO","(LUMO+6)-LUMO","inv(HOMO-(HOMO-1))","inv(HOMO-(HOMO-2))","inv(HOMO-(HOMO-3))","inv(HOMO-(HOMO-4))","inv(HOMO-(HOMO-5))","inv(HOMO-(HOMO-6))","inv((LUMO+1)-LUMO)","inv((LUMO+2)-LUMO)","inv((LUMO+3)-LUMO)","inv((LUMO+4)-LUMO)","inv((LUMO+5)-LUMO)","inv((LUMO+6)-LUMO)","inv(LUMO-(HOMO-1))","inv(LUMO-(HOMO-2))","inv((LUMO+1)-HOMO)","inv((LUMO+2)-HOMO)"]
    list_regressors_NN_LUMOr = ["Gap_PBE","avgempty(eV-1)" ,"avgoccupied(eV-1)" ,"(LUMO+1)-LUMO",  "inv((LUMO+2)-LUMO)","inv((LUMO+1)-LUMO)"  ,   ] #'ratio_bonds_hybridized',  'ratio_bonds_with_order_till_1p1_eq_nonhybridized',
    list_regressors_NN_gapr  = ["Gap_PBE","Gap_B3LYP","avgempty(eV-1)", "avgoccupied(eV-1)",  "inv(HOMO-(HOMO-1))", "inv(HOMO-(HOMO-2))", "inv(HOMO-(HOMO-3))", "inv((LUMO+3)-LUMO)", "inv((LUMO+2)-LUMO)","inv((LUMO+1)-LUMO)"  ,   ] #'ratio_bonds_hybridized',  'ratio_bonds_with_order_till_1p1_eq_nonhybridized',

    list_regressors_KNN_HOMOr = [ "avgoccupied(eV-1)", "avgempty(eV-1)" ]
    list_regressors_KNN_LUMOr = [ "avgoccupied(eV-1)", "avgempty(eV-1)" ]
    list_regressors_KNN_gapr  = [ "avgoccupied(eV-1)", "avgempty(eV-1)" ]
    '''

else: # alternative_regressors

    # list_regr_alt = ["Gap_PBE", "avgoccupied(eV-1)", "avgempty(eV-1)","(LUMO+1)-LUMO", "(LUMO+2)-LUMO", "(LUMO+3)-LUMO", "(LUMO+4)-LUMO","(LUMO+5)-LUMO","inv((LUMO+1)-LUMO)","inv((LUMO+2)-LUMO)", "inv((LUMO+3)-LUMO)", "inv((LUMO+4)-LUMO)", "inv((LUMO+5)-LUMO)","inv(LUMO-(HOMO-1))", "inv(LUMO-(HOMO-2))", "inv((LUMO+1)-HOMO)","inv((LUMO+2)-HOMO)"]

    # list_regr_alt = ["Gap_B3LYP","Gap_PBE","avgoccupied(eV-1)","avgempty(eV-1)","deg HOMO","deg LUMO","ElecEigvalMean","ElecEigvalVariance","ElecEigvalSkewness","ElecEigvalKurtosis"]
    # list_regr_alt = ["Gap_B3LYP","Gap_PBE","avgoccupied(eV-1)","avgempty(eV-1)","deg HOMO","deg LUMO","ElecEigvalMean","ElecEigvalVariance","ElecEigvalSkewness","ElecEigvalKurtosis","HOMO-(HOMO-1)","HOMO-(HOMO-2)","HOMO-(HOMO-3)","HOMO-(HOMO-4)","HOMO-(HOMO-5)","HOMO-(HOMO-6)","(LUMO+1)-LUMO","(LUMO+2)-LUMO","(LUMO+3)-LUMO","(LUMO+4)-LUMO","(LUMO+5)-LUMO","(LUMO+6)-LUMO","inv(HOMO-(HOMO-1))","inv(HOMO-(HOMO-2))","inv(HOMO-(HOMO-3))","inv(HOMO-(HOMO-4))","inv(HOMO-(HOMO-5))","inv(HOMO-(HOMO-6))","inv((LUMO+1)-LUMO)","inv((LUMO+2)-LUMO)","inv((LUMO+3)-LUMO)","inv((LUMO+4)-LUMO)","inv((LUMO+5)-LUMO)","inv((LUMO+6)-LUMO)","inv(LUMO-(HOMO-1))","inv(LUMO-(HOMO-2))","inv((LUMO+1)-HOMO)","inv((LUMO+2)-HOMO)"]
    # list_regr_alt = ["1stPhFreq","LastPhFreq","PhFreqMean","PhFreqVariance","PhFreqSkewness","PhFreqKurtosis"]

    list_regr_alt = ["bond_length_angstrom_average","bond_length_angstrom_standard_deviation","bond_length_skewness","bond_length_excess_Kurtosis","bond_length_1st_lowest_value","bond_length_2nd_lowest_value","bond_length_3rd_lowest_value","bond_length_4th_lowest_value","bond_length_4th_highest_value","bond_length_3rd_highest_value","bond_length_2nd_highest_value","bond_length_1st_highest_value"]
    #list_regr_alt += ["Mayer_bond_order_average","Mayer_bond_order_standard_deviation","Mayer_bond_order_skewness","Mayer_bond_order_excess_Kurtosis","Mayer_bond_order_1st_lowest_value","Mayer_bond_order_2nd_lowest_value","Mayer_bond_order_3rd_lowest_value","Mayer_bond_order_4th_lowest_value","Mayer_bond_order_4th_highest_value","Mayer_bond_order_3rd_highest_value","Mayer_bond_order_2nd_highest_value","Mayer_bond_order_1st_highest_value"]
    #list_regr_alt += ["average_GJ_bond_order","standard_deviation_GJ_bond_order","skewness_GJ_bond_order","excess_Kurtosis_GJ_bond_order","GJ_Nbond_with_order_till_1p1_eq_nonhybridized","GJ_Nbond_with_order_till_1p2","GJ_Nbond_with_order_till_1p3","GJ_Nbond_with_order_till_1p4","GJ_Nbond_with_order_till_1p5","GJ_Nbond_with_order_over_1p5","GJ_Nbond_hybridized","GJ_ratio_bonds_with_order_till_1p1_eq_nonhybridized","GJ_ratio_bonds_with_order_till_1p2","GJ_ratio_bonds_with_order_till_1p3","GJ_ratio_bonds_with_order_till_1p4","GJ_ratio_bonds_with_order_till_1p5","GJ_ratio_bonds_with_order_over_1p5","GJ_ratio_bonds_hybridized","1st_lowest_value_GJ_bond_order","2nd_lowest_value_GJ_bond_order","3rd_lowest_value_GJ_bond_order","4th_lowest_value_GJ_bond_order","4th_highest_value_GJ_bond_order","3rd_highest_value_GJ_bond_order","2nd_highest_value_GJ_bond_order","1st_highest_value_GJ_bond_order"]
    #list_regr_alt += ["average_NM1_bond_order","standard_deviation_NM1_bond_order","skewness_NM1_bond_order","excess_Kurtosis_NM1_bond_order","NM1_Nbond_with_order_till_1p1_eq_nonhybridized","NM1_Nbond_with_order_till_1p2","NM1_Nbond_with_order_till_1p3","NM1_Nbond_with_order_till_1p4","NM1_Nbond_with_order_till_1p5","NM1_Nbond_with_order_over_1p5","NM1_Nbond_hybridized","NM1_ratio_bonds_with_order_till_1p1_eq_nonhybridized","NM1_ratio_bonds_with_order_till_1p2","NM1_ratio_bonds_with_order_till_1p3","NM1_ratio_bonds_with_order_till_1p4","NM1_ratio_bonds_with_order_till_1p5","NM1_ratio_bonds_with_order_over_1p5","NM1_ratio_bonds_hybridized","1st_lowest_value_NM1_bond_order","2nd_lowest_value_NM1_bond_order","3rd_lowest_value_NM1_bond_order","4th_lowest_value_NM1_bond_order","4th_highest_value_NM1_bond_order","3rd_highest_value_NM1_bond_order","2nd_highest_value_NM1_bond_order","1st_highest_value_NM1_bond_order"]
    #list_regr_alt += ["average_NM2_bond_order","standard_deviation_NM2_bond_order","skewness_NM2_bond_order","excess_Kurtosis_NM2_bond_order","NM2_Nbond_with_order_till_1p1_eq_nonhybridized","NM2_Nbond_with_order_till_1p2","NM2_Nbond_with_order_till_1p3","NM2_Nbond_with_order_till_1p4","NM2_Nbond_with_order_till_1p5","NM2_Nbond_with_order_over_1p5","NM2_Nbond_hybridized","NM2_ratio_bonds_with_order_till_1p1_eq_nonhybridized","NM2_ratio_bonds_with_order_till_1p2","NM2_ratio_bonds_with_order_till_1p3","NM2_ratio_bonds_with_order_till_1p4","NM2_ratio_bonds_with_order_till_1p5","NM2_ratio_bonds_with_order_over_1p5","NM2_ratio_bonds_hybridized","1st_lowest_value_NM2_bond_order","2nd_lowest_value_NM2_bond_order","3rd_lowest_value_NM2_bond_order","4th_lowest_value_NM2_bond_order","4th_highest_value_NM2_bond_order","3rd_highest_value_NM2_bond_order","2nd_highest_value_NM2_bond_order","1st_highest_value_NM2_bond_order"]
    list_regr_alt += ["average_NM3_bond_order","standard_deviation_NM3_bond_order","skewness_NM3_bond_order","excess_Kurtosis_NM3_bond_order","NM3_Nbond_with_order_till_1p1_eq_nonhybridized","NM3_Nbond_with_order_till_1p2","NM3_Nbond_with_order_till_1p3","NM3_Nbond_with_order_till_1p4","NM3_Nbond_with_order_till_1p5","NM3_Nbond_with_order_over_1p5","NM3_Nbond_hybridized","NM3_ratio_bonds_with_order_till_1p1_eq_nonhybridized","NM3_ratio_bonds_with_order_till_1p2","NM3_ratio_bonds_with_order_till_1p3","NM3_ratio_bonds_with_order_till_1p4","NM3_ratio_bonds_with_order_till_1p5","NM3_ratio_bonds_with_order_over_1p5","NM3_ratio_bonds_hybridized","1st_lowest_value_NM3_bond_order","2nd_lowest_value_NM3_bond_order","3rd_lowest_value_NM3_bond_order","4th_lowest_value_NM3_bond_order","4th_highest_value_NM3_bond_order","3rd_highest_value_NM3_bond_order","2nd_highest_value_NM3_bond_order","1st_highest_value_NM3_bond_order"]

    print(" ==> USING ("+str(len(list_regr_alt))+") ALTERNATIVE REGRESSORS: ",list_regr_alt )

    list_regressors_HOMOr = list_regr_alt.copy()
    list_regressors_LUMOr = list_regr_alt.copy()
    list_regressors_gapr  = list_regr_alt.copy()

    list_regressors_lin_ols_HOMO = list_regressors_HOMOr.copy()
    list_regressors_lin_ols_LUMO = list_regressors_LUMOr.copy()
    list_regressors_lin_ols_gapr = list_regressors_gapr.copy()
    

    list_regressors_RF_HOMOr = list_regressors_HOMOr.copy()
    list_regressors_RF_LUMOr = list_regressors_LUMOr.copy()
    list_regressors_RF_gapr  = list_regressors_gapr.copy()

    list_regressors_NN_HOMOr = list_regressors_HOMOr.copy()
    list_regressors_NN_LUMOr = list_regressors_LUMOr.copy()
    list_regressors_NN_gapr = list_regressors_gapr.copy()

    list_regressors_KNN_HOMOr = list_regressors_HOMOr.copy()
    list_regressors_KNN_LUMOr = list_regressors_LUMOr.copy()
    list_regressors_KNN_gapr = list_regressors_gapr.copy()


list_regressors_all = ['Gap ren (output)', 'Nat', 'Gap_B3LYP', 'Gap_PBE', 'avgoccupied(eV-1)', 'avgempty(eV-1)', 'deg HOMO', 'deg LUMO', 'ElecEigvalMean', 'ElecEigvalVariance', 'ElecEigvalSkewness', 'ElecEigvalKurtosis', 'HOMO-(HOMO-1)', 'HOMO-(HOMO-2)', 'HOMO-(HOMO-3)', 'HOMO-(HOMO-4)', 'HOMO-(HOMO-5)', 'HOMO-(HOMO-6)', '(LUMO+1)-LUMO', '(LUMO+2)-LUMO', '(LUMO+3)-LUMO', '(LUMO+4)-LUMO', '(LUMO+5)-LUMO', '(LUMO+6)-LUMO', 'inv(HOMO-(HOMO-1))', 'inv(HOMO-(HOMO-2))', 'inv(HOMO-(HOMO-3))', 'inv(HOMO-(HOMO-4))', 'inv(HOMO-(HOMO-5))', 'inv(HOMO-(HOMO-6))', 'inv((LUMO+1)-LUMO)', 'inv((LUMO+2)-LUMO)', 'inv((LUMO+3)-LUMO)', 'inv((LUMO+4)-LUMO)', 'inv((LUMO+5)-LUMO)', 'inv((LUMO+6)-LUMO)', 'inv(LUMO-(HOMO-1))', 'inv(LUMO-(HOMO-2))', 'inv((LUMO+1)-HOMO)', 'inv((LUMO+2)-HOMO)', 'num hexagons', 'Volume', 'Area', 'Area/volume ratio', '1stPhFreq', 'LastPhFreq', 'PhFreqMean', 'PhFreqVariance', 'PhFreqSkewness', 'PhFreqKurtosis', 'bond_length_angstrom_average', 'bond_length_angstrom_standard_deviation', 'bond_length_skewness', 'bond_length_excess_Kurtosis', 'bond_length_1st_lowest_value', 'bond_length_2nd_lowest_value', 'bond_length_3rd_lowest_value', 'bond_length_4th_lowest_value', 'bond_length_4th_highest_value', 'bond_length_3rd_highest_value', 'bond_length_2nd_highest_value', 'bond_length_1st_highest_value', 'average_NM3_bond_order', 'standard_deviation_NM3_bond_order', 'skewness_NM3_bond_order', 'excess_Kurtosis_NM3_bond_order', 'Nbond_with_order_till_1p1_eq_nonhybridized', 'Nbond_with_order_till_1p2', 'Nbond_with_order_till_1p3', 'Nbond_with_order_till_1p4', 'Nbond_with_order_till_1p5', 'Nbond_with_order_over_1p5', 'Nbond_hybridized', 'ratio_bonds_with_order_till_1p1_eq_nonhybridized', 'ratio_bonds_with_order_till_1p2', 'ratio_bonds_with_order_till_1p3', 'ratio_bonds_with_order_till_1p4', 'ratio_bonds_with_order_till_1p5', 'ratio_bonds_with_order_over_1p5', 'ratio_bonds_hybridized', '1st_lowest_value_NM3_bond_order', '2nd_lowest_value_NM3_bond_order', '3rd_lowest_value_NM3_bond_order', '4th_lowest_value_NM3_bond_order', '4th_highest_value_NM3_bond_order', '3rd_highest_value_NM3_bond_order', '2nd_highest_value_NM3_bond_order', '1st_highest_value_NM3_bond_order']

errors_df  = DataFrame(columns=["Y_full_ref","Y_full_val","y_ref","y_val","y_err"])
errors_df0 = DataFrame( )

# ['Name', 'Symmetry (generic)', 'Symmetry (specific)', 'HOMO ren (output)', 'LUMO ren (output)', 'Gap ren (output)', 'Nat', 'Gap_B3LYP', 'Gap_PBE', 'avgoccupied(eV-1)', 'avgempty(eV-1)', 'deg HOMO', 'deg LUMO', 'ElecEigvalMean', 'ElecEigvalVariance', 'ElecEigvalSkewness', 'ElecEigvalKurtosis', 'HOMO-(HOMO-1)', 'HOMO-(HOMO-2)', 'HOMO-(HOMO-3)', 'HOMO-(HOMO-4)', 'HOMO-(HOMO-5)', 'HOMO-(HOMO-6)', '(LUMO+1)-LUMO', '(LUMO+2)-LUMO', '(LUMO+3)-LUMO', '(LUMO+4)-LUMO', '(LUMO+5)-LUMO', '(LUMO+6)-LUMO', 'inv(HOMO-(HOMO-1))', 'inv(HOMO-(HOMO-2))', 'inv(HOMO-(HOMO-3))', 'inv(HOMO-(HOMO-4))', 'inv(HOMO-(HOMO-5))', 'inv(HOMO-(HOMO-6))', 'inv((LUMO+1)-LUMO)', 'inv((LUMO+2)-LUMO)', 'inv((LUMO+3)-LUMO)', 'inv((LUMO+4)-LUMO)', 'inv((LUMO+5)-LUMO)', 'inv((LUMO+6)-LUMO)', 'inv(LUMO-(HOMO-1))', 'inv(LUMO-(HOMO-2))', 'inv((LUMO+1)-HOMO)', 'inv((LUMO+2)-HOMO)', 'num hexagons', 'Volume', 'Area', 'Area/volume ratio', '1stPhFreq', 'LastPhFreq', 'PhFreqMean', 'PhFreqVariance', 'PhFreqSkewness', 'PhFreqKurtosis', 'bond_length_angstrom_average', 'bond_length_angstrom_standard_deviation', 'bond_length_skewness', 'bond_length_excess_Kurtosis', 'bond_length_1st_lowest_value', 'bond_length_2nd_lowest_value', 'bond_length_3rd_lowest_value', 'bond_length_4th_lowest_value', 'bond_length_4th_highest_value', 'bond_length_3rd_highest_value', 'bond_length_2nd_highest_value', 'bond_length_1st_highest_value', 'average_NM3_bond_order', 'standard_deviation_NM3_bond_order', 'skewness_NM3_bond_order', 'excess_Kurtosis_NM3_bond_order', 'Nbond_with_order_till_1p1_eq_nonhybridized', 'Nbond_with_order_till_1p2', 'Nbond_with_order_till_1p3', 'Nbond_with_order_till_1p4', 'Nbond_with_order_till_1p5', 'Nbond_with_order_over_1p5', 'Nbond_hybridized', 'ratio_bonds_with_order_till_1p1_eq_nonhybridized', 'ratio_bonds_with_order_till_1p2', 'ratio_bonds_with_order_till_1p3', 'ratio_bonds_with_order_till_1p4', 'ratio_bonds_with_order_till_1p5', 'ratio_bonds_with_order_over_1p5', 'ratio_bonds_hybridized', '1st_lowest_value_NM3_bond_order', '2nd_lowest_value_NM3_bond_order', '3rd_lowest_value_NM3_bond_order', '4th_lowest_value_NM3_bond_order', '4th_highest_value_NM3_bond_order', '3rd_highest_value_NM3_bond_order', '2nd_highest_value_NM3_bond_order', '1st_highest_value_NM3_bond_order']

#-----------------------------------------------------------------------------------------------------------------------

class Configuration:
    '''This class contains the information on the regressors and other variables necessary as inputs for the calculation'''

    def __init__(self, quantity_to_forecast, ml_method, forecast_of_residuals=True):

        self.quantity_to_forecast = quantity_to_forecast    # Either "HOMO","LUMO" or "gap"
        self.ml_method = ml_method                          # "Lasso", "Neural networks", "Random forests", "Kernel ridge regression", "KNN" or "Support vector machine"
        self.forecast_of_residuals = forecast_of_residuals  # forecast_of_residuals = True   # If True, application of ML methods on top of linear regression (i.e. using residuals of the regression as output)

        if (quantity_to_forecast == "HOMO"):
            self.field_out_to_read = "HOMO ren (output)"
            self.regressor_fields_regression = list_regressors_lin_ols_HOMO
            if ((ml_method == "NN") or (ml_method == "Neural networks")):
                self.regressor_fields_ml = list_regressors_NN_HOMOr
            elif ((ml_method == "KNN")):
                self.regressor_fields_ml = list_regressors_KNN_HOMOr #self.regressor_fields_regression
            else:
                self.regressor_fields_ml = list_regressors_RF_HOMOr

        elif (quantity_to_forecast == "LUMO"):
            self.field_out_to_read = "LUMO ren (output)"
            self.regressor_fields_regression = list_regressors_lin_ols_LUMO
            if ((ml_method == "NN") or (ml_method == "Neural networks")):
                self.regressor_fields_ml = list_regressors_NN_LUMOr
            elif ((ml_method == "KNN")):
                self.regressor_fields_ml = list_regressors_KNN_LUMOr #self.regressor_fields_regression
            else:
                self.regressor_fields_ml = list_regressors_RF_LUMOr

        elif (quantity_to_forecast == "gap"):
            self.field_out_to_read = "Gap ren (output)"
            self.regressor_fields_regression = list_regressors_lin_ols_gapr
            if ((ml_method == "NN") or (ml_method == "Neural networks")):
                self.regressor_fields_ml = list_regressors_NN_gapr
            elif ((ml_method == "KNN")):
                self.regressor_fields_ml = list_regressors_KNN_gapr
            else:
                self.regressor_fields_ml = list_regressors_RF_gapr

        else:
            print("   ERROR: Unrecognized quantity to forecast: " + quantity_to_forecast + "; Please, check the module_io.py file.\n");
            exit(1)

        if (self.regressor_fields_ml == self.regressor_fields_regression):
            self.different_regressors = True
        else:
            self.different_regressors = False

        # The variable << param_sweep >> is the list of values to be considered in the sweeping. These values determine how the ML algorithm is defined.
        if (ml_method == "Lasso"):
            self.param_sweep = [0.000001, 0.1, 0.5, 1.0, 1.5]  # alpha
        elif ((ml_method == "NN") or (ml_method == "Neural networks")):
            self.param_sweep = [400]  # hidden_layer_sizes [300, 200, 100]
        elif ((ml_method == "RF") or (ml_method == "Random forest") or (ml_method == "Random forests")):
            self.param_sweep = [  250 ]  #  n_estimators  [ 3, 10, 50, 100, 200, 300]
        elif ((ml_method == "KRR") or (ml_method == "Kernel ridge regression") or (ml_method == "Kernel ridge")):
            self.param_sweep = [1, 0.000001, 0.1, 0.5, 1.0, 2]  # alpha
        elif ((ml_method == "SVM") or (ml_method == "Support vector machine")):
            self.param_sweep = [1]
        elif ((ml_method == "KNN") or (ml_method == "K-nearest neighbours")):
            self.param_sweep = [22]#[ 2, 4, 10, 15, 18, 20, 22, 24, 26, 28, 30, 32, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 90, 100, 120 ]  # [ 40 ]
        else:
            self.param_sweep = [1]

        self.filepath_errors_train = f"{execution_directory}Results/ml_training_errors"
        self.filepath_errors_test = f"{execution_directory}Results/ml_test_errors"

        if (forecast_of_residuals): text = "fromresid_"
        else: text = "NOTfromresid"

        if (len(self.param_sweep)==1):

            self.boxplot_prefix = "ml_test_errors_" + ml_method + "_" +quantity_to_forecast + "__" + text + str(self.param_sweep[0]) #box_plot_quantiles(  "ml_test_errors_RF_HOMO__fromresid_250"

            if (alternative_regressors):
                regressor_list = list_regr_alt.copy()
            else:
                if ((ml_method == "NN") or (ml_method == "SVM")):
                    if (quantity_to_forecast == "HOMO"): regressor_list = list_regressors_NN_HOMOr.copy()
                    if (quantity_to_forecast == "LUMO"): regressor_list = list_regressors_NN_LUMOr.copy()
                    if (quantity_to_forecast == "gap"): regressor_list  = list_regressors_NN_gapr.copy()
                elif (ml_method == "KNN"):
                    if (quantity_to_forecast == "HOMO"): regressor_list = list_regressors_KNN_HOMOr.copy()
                    if (quantity_to_forecast == "LUMO"): regressor_list = list_regressors_KNN_LUMOr.copy()
                    if (quantity_to_forecast == "gap"): regressor_list  = list_regressors_KNN_gapr.copy()
                else: #(ml_method in ["RF","GB", "SVR", "Lasso"]) :
                    if (quantity_to_forecast == "HOMO"): regressor_list = list_regressors_RF_HOMOr.copy()
                    if (quantity_to_forecast == "LUMO"): regressor_list = list_regressors_RF_LUMOr.copy()
                    if (quantity_to_forecast == "gap"): regressor_list  = list_regressors_RF_gapr.copy()

            regressor = regressor_list[-1]
            regressor = str(regressor)
            for spec in ['"', '/', '\\', ' ']:
                regressor = regressor.replace(spec, "")
            self.boxplot_prefix += "_" + regressor

            del  regressor_list; del regressor; del spec

        else:
            self.boxplot_prefix=None




#-----------------------------------------------------------------------------------------------------------------------
