# ML_fullerenes_renormalization_code_paper
This repository contains the code used for the ML-regression calculations of the article 'Electron-vibrational renormalization in fullerenes through ab initio and machine learning methods' by Pablo García-Risueno et al. The code is easy to run, one just needs to download it and set the desired values of the input parameters in the file module_io.py. In the first block of this file one finds some parameters with their explanations in the comments. Then one finds the lists of regressors to be used in the calculations. 
The user can also make small modifications in the main.py file. In it he will find the three lines that follow:

    for quantity in ["LUMO","gap"]: # ["HOMO","LUMO","gap"]

        for forecast_of_residuals in [True,False]:  # [True,False]
    
            for method in ["RF","NN", "KNN"]: # ["RF","NN", "KNN","GB"]

The first one specifies which quantities are analysed (HOMO, LUMO or gap, and any combination of them). The second line specifies whether or not the Machine Learning-based forecast must be performed on top of a linear (ols) regression, or directly on the input data (set to True or to False, respectively, for these possibilities). The third line specifies which machine learning method will be used. You can use Random Forests (RF), Neural Networks (NN), k-nearest Neigbours (KNN), Gradient Boosting (GB), Lasso regression ("Lasso"), Kernel Ridge Regression (KRR) or Support Vector Machines (SVM). The parameters for the corresponding ML regression method must be specified in module_io.py, in the class called Configuration. The variable called self.param_sweep specifies the parameters for all solvers, see module_training.py to see how this is used. For example, if the method is RF, then this parameter specifies the list of trees to be used ([400] indicates that just one value of the parameter is considered -400 trees-). Further parameters of the ML solver can be modified in module_training.py, search for "ml_variable".
