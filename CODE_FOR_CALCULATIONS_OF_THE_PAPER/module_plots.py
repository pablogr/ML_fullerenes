
import pandas as pd
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from os import makedirs, path
from scipy.interpolate import make_interp_spline
#pd.set_option('max_columns', 20); pd.set_option('max_rows', 99)
from scipy.stats import norm, nct

#-----------------------------------------------------------------------------------------------------------------------
import module_io


def plot_2D( datafile_path, output_dir, plot_name, quantity_x, quantity_y, quantity_y2=None, rescaling_y2=1, mytitle=None  ):
    '''This function makes a 2D plot whose variable x is "quantity_x" and whos variable y is the optimum (maximum or
     minimum of "quantity_y").'''

    list_product_labels_name = output_dir.replace("Time_series/Spreads/Spreads_","");  list_product_labels_name = list_product_labels_name.replace("/Plots","")
    if ( list_product_labels_name=="cryptocurrencies" ): my_interval=91#365
    else: my_interval =252

    lt1 = plot_name.replace("spr_",""); lt2 = lt1.split("_")[1] ;  lt1 = lt1.split("_")[0]
    axes_text = {"Spread_y_vs_x":"Spread","Spread_x_vs_y":"Spread", "Sharpe_ratio":"$SR$", "Sharpe_ratio_with_semideviation":"$SR'$", "Sharpe_ratio_from_semideviation":"$SR'$", "profit_taking_param":"profit-taking ($)"}
    if (quantity_y=="Spread_y_vs_x"):
        legend_text = lt2+"-vs-"+lt1
    elif (quantity_y=="Spread_x_vs_y"):
        legend_text = lt1+"-vs-"+lt2
    else:
        try:
            legend_text = axes_text[quantity_y]
        except KeyError:
            legend_text = quantity_y

    # Reading the data
    df0 = pd.read_csv( datafile_path, header=0 , usecols=[quantity_x, quantity_y] )
    print(datafile_path,"\n",df0)
    arr_x = df0[quantity_x]; arr_y = df0[quantity_y]
    text=""

    # Actual plotting
    fig, ax = plt.subplots()
    try:
        label_x = axes_text[quantity_x]
    except KeyError:
        label_x = quantity_x
    plt.xlabel( label_x, fontsize=16)
    try:
        plt.ylabel( axes_text[quantity_y], fontsize=16) #plt.ylabel( "Sharpe ratio", fontsize=16)
    except KeyError:
        plt.ylabel( quantity_y, fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=11)
    ax.tick_params(axis='both', which='minor', labelsize=11)
    plt.xticks(rotation=45, ha='right')

    if (quantity_x=="Date"):
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=my_interval))
        #ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))

    color1 = 'b'
    color2 = '#8b0000'

    #ax.plot( arr_x, arr_y, '-', label=legend_text, color=color1, linewidth=1.2, markersize=2 )
    #ax.plot( [arr_x[0],arr_x[660]], [1,1], '-',   color='#ff4500', linewidth=3.5 )
    #ax.plot([arr_x[720], arr_x[2380]], [1.18, 1.18], '-',  color='#ff4500', linewidth=3.5)
    #ax.plot([arr_x[2510], arr_x[len(arr_x)-1]], [1.412, 1.412], '-',  color='#ff4500', linewidth=3.5)
    ax.plot(arr_x, arr_y, "o", label=legend_text, color=color1, linewidth=2, markersize=5)
    text = "-and-" + str(quantity_y2)
    try:
        x_smooth = np.linspace(arr_x.min(), arr_x.max(), 300)
        y_smooth = np.array(make_interp_spline(arr_x, arr_y)(x_smooth))
        ax.plot(x_smooth, y_smooth, '-', color=color1, linewidth=2.2)
        del x_smooth; del y_smooth
    except ValueError:
        print("\nWARNING: Unable to perform the interpolation of " + quantity_y2 + "-vs-" + quantity_x + "\n")


    if (quantity_y2!=None):
        df1 = pd.read_csv( datafile_path, header=0, usecols=[quantity_x, quantity_y2] )
        print(datafile_path,"\n",df1)
        arr_x2 = df1[quantity_x]; arr_y2 = df1[quantity_y2]
        arr_y2 *= rescaling_y2
        if (quantity_y2 == "Spread_y_vs_x"):
            legend_text2 = lt2 + "-vs-" + lt1
        elif (quantity_y2 == "Spread_x_vs_y"):
            legend_text2 = lt1 + "-vs-" + lt2
        else:
            try:
                legend_text2 = axes_text[quantity_y2]
            except KeyError:
                legend_text2 = quantity_y2
        #if (quantity_y2 == "Sharpe_ratio_with_semideviation"): my_label2 = ("Sharpe ratio from\nsemi-deviation "+r'($\times 1/\sqrt{2}$' + ")")
        ax.plot(arr_x2, arr_y2, "s", label=legend_text2, color=color2, linewidth=2, markersize=5)
        text = "-and-" + str(quantity_y2)
        try:
            x2_smooth = np.linspace(arr_x2.min(), arr_x2.max(), 300)
            y2_smooth = np.array( make_interp_spline(arr_x2, arr_y2)(x2_smooth) )
            ax.plot(x2_smooth, y2_smooth, '--', color=color2, linewidth=2)
            del x2_smooth; del y2_smooth
        except ValueError:
            print("\nWARNING: Unable to perform the interpolation of " + quantity_y2 + "-vs-" + quantity_x + "\n")

        del arr_x2; del arr_y2; del df1

    if (mytitle!=None):
        ax.set_title(mytitle,fontsize=16 )
    plt.legend(loc='best', labelspacing=1.5) # loc='lower right', loc='best', bbox_to_anchor=(0.5, 0.1, 0.5, 0.5),
    plot_path = output_dir+"/"+plot_name +"_"+quantity_y+".pdf"
    plt.savefig(plot_path, format="pdf", bbox_inches="tight")
    plt.clf(); plt.cla(); plt.close('all')

    del quantity_x; del quantity_y; del ax; del fig; del arr_y; del arr_x;

    return

#-----------------------------------------------------------------------------------------------

def my_hist( bins, array_in ):
    '''This function calculates the count for a histogram'''

    #array_in = np.sort( array_in )
    my_hist_x = np.zeros(len(bins)-1,dtype=float)
    my_hist_y = np.zeros(len(bins) - 1,dtype=float)
    size = len(array_in)

    if ((array_in[0]<bins[0]) or (array_in[-1]>bins[-1]) ):
        print("\nWARNING: The limits of the bins ("+str(bins)+") are exceeded by the input array ("+str(array_in)+").\nPlease, modify the code to make it suitable to plotting this; the plot is omitted.\n")
        del bins; del array_in;  del size
        return my_hist_x, my_hist_y, False

    for i in range(1,len(bins)):
        bin0 = bins[i-1]
        bin1 = bins[i]
        count = 0
        for j in range(len(array_in)):
            if (( array_in[j] >= bin0) and ( array_in[j] < bin1) ): count+=1
        my_hist_y[i-1] = count

    for i in range(len(my_hist_y)):
        my_hist_x[i] = (bins[i]+bins[i+1])/2
        my_hist_y[i] /= size
        #print(bins[i],bins[i+1],":", my_hist_x[i], my_hist_y[i])

    del bins; del array_in; del count; del bin0; del bin1; del i; del j; del size

    return my_hist_x, my_hist_y, True

#-----------------------------------------------------------------------------------------------

def integrated_plot( filepath_train, filepath_test, quantity, method, forecast_of_residuals, dir_out=module_io.execution_directory  ):
    '''This function makes a plot of y_train-vs-x_train and y_test-vs-x_test, including histograms with the vec_train and vec_test
    (which contain the errors).'''

    from math import ceil
    import importlib
    importlib.import_module('mpl_toolkits').__path__
    #from mpl_toolkits.axes_grid.inset_locator import inset_axes
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    truncation = 300

    df_train = pd.read_csv(filepath_train,header=0)
    df_test  = pd.read_csv(filepath_test, header=0)

    x_train = np.array( df_train["Y_full_ref"] ); x_train = x_train[0:min(truncation,len(x_train)-1)]
    y_train = np.array( df_train["Y_full_val"] ); y_train = y_train[0:min(truncation,len(y_train)-1)]
    vec_train = np.array( df_train["y_err"] )


    x_test = np.array(df_test["Y_full_ref"]); x_test = x_test[0:min(int(truncation/2),len(x_test)-1)]
    y_test = np.array(df_test["Y_full_val"]);
    y_test = y_test[0:min(int(truncation/2),len(y_test)-1)];
    vec_test = np.array(df_test["y_err"])

    '''
    x_train = 20 + 180 * uniform.rvs(size=150)
    y_train = x_train + norm.rvs(size=150, loc=0, scale=15)
    x_test = 20 + 180 * uniform.rvs(size=35)
    y_test = x_test + norm.rvs(size=35, loc=0, scale=25)
    vec_train = norm.rvs(size=30, loc=-1, scale=10)
    vec_test = norm.rvs(size=150, loc=1, scale=10)
    '''

    if (quantity == "HOMO"):
        color1 = '#0000ff'
        color2 = '#00ced1'
    elif (quantity == "LUMO"):
        color1 = '#228b22'
        color2 = '#00ff00'
    else:
        color1 = '#9400d3'#'b'
        color2 = 'r'#'#00ff00'#7f 'r'

    fig, ax = plt.subplots()
    ax.set_box_aspect(1)  # For square frame


    # CENTRAL PLOT

    ax.plot([-1000, 1000], [-1000, 1000], "-", color='k', linewidth=1.2, alpha=0.6)               # Diagonal line
    ax.scatter( x_train, y_train, s=41, marker='s', label='Training', color=color1,    alpha=0.6) # Data train
    ax.scatter( x_test,  y_test,  s=37, marker='o', label='Test', color=color2,   alpha=0.5)      # Data tests

    ax.set_xlabel( 'Ab initio (meV)', fontsize=14)
    ax.set_ylabel('ML-predicted (meV)', fontsize=14)
    if (method=="RF"):
        text = "Random forests"
    elif (method=="NN"):
        text = "Neural networks"
    elif (method=="KNN"):
        text = "k-nearest neighbors"
    if (quantity=="gap"):
        ax.set_title('Gap renormalization - '+text, fontsize=15)
    else:
        ax.set_title( quantity + ' renormalization - '+text, fontsize=15)
    #legend = plt.legend(loc='best', bbox_to_anchor=(0.28, 0.50, 0.5, 0.5),edgecolor="gray")  # , labelspacing=1.5 loc='lower right', loc='best', bbox_to_anchor=(0.5, 0.1, 0.5, 0.5),
    #legend = plt.legend(loc='best', bbox_to_anchor=(0.23, 0.0, 0.5, 0.5), edgecolor="gray")
    legend = plt.legend(loc='best', bbox_to_anchor=(0.262, 0.49, 0.5, 0.5), edgecolor="gray")
    legend.get_frame().set_alpha(None)


    # Plot limits

    roundto = 25
    lim0a = int(1.02 * min( min(x_train), min(y_train) ))
    lim0a -= (lim0a % roundto)
    lim1a = ceil(1.02 * max( max(x_train), max(y_train)))
    lim1a += (lim1a % roundto)
    if ((lim1a%roundto )!=0): lim1a = lim1a +roundto  - (lim1a%roundto )

    lim0b = int(1.02 * min( min(x_test), min(y_test) ))
    lim0b -= (lim0b % roundto)
    lim1b = ceil(1.02 * max( max(x_test), max(y_test)))
    if ((lim1b%roundto )!=0): lim1b = lim1b +roundto  - (lim1b%roundto )

    lim0 = min(lim0a,lim0b)
    lim1 = min(lim1a,lim1b)

    plt.xlim([lim0,lim1]); plt.ylim([lim0,lim1])
    my_tics = np.arange( lim0, lim1+1, roundto )
    ax.set_xticks(my_tics)
    ax.set_yticks(my_tics)
    ax.tick_params(axis='both', which='major', labelsize=11)
    ax.tick_params(axis='both', which='minor', labelsize=11)
    plt.xticks(rotation=45, ha='right')


    # PLOT OF HISTOGRAMS

    lim_bins = 32.5
    my_tics_hist = np.arange(-30, 30.000000001, 10 ) #old: np.arange(-30, 31, 10)
    bins = np.arange(-lim_bins,lim_bins+0.000000001, 5)
    width = 0.8 * (bins[1] - bins[0])

    my_hist_x_train, my_hist_y_train, ok = my_hist(bins, vec_train)
    if (not ok): return
    my_hist_x_test,  my_hist_y_test, ok  = my_hist(bins, vec_test)
    if (not ok): return

    ax_hist1 = inset_axes( ax, width="30%", height="30%", loc='upper left', bbox_to_anchor=(0.13,-0.01,1,1), bbox_transform=ax.transAxes)  #loc='best', bbox_to_anchor=(0.5, 0.1, 0.5, 0.5)
    ax_hist1.bar( my_hist_x_train, my_hist_y_train, width, label='Histogram', color=color1)
    ax_hist1.tick_params(axis='both', which='major', labelsize=9)
    ax_hist1.set_xlim([-33,33])
    ax_hist1.set_xticks(my_tics_hist)
    ax_hist1.set_xlabel('Train. errors (meV) $\quad$    ', fontsize=10); ax_hist1.set_ylabel('Probability', fontsize=10)
    plt.xticks(rotation=90 )

    ax_hist2 = inset_axes(ax, width="30%", height="30%", loc='lower right', bbox_to_anchor=(0, 0.14, 1, 1),bbox_transform=ax.transAxes)  # loc='best', bbox_to_anchor=(0.5, 0.1, 0.5, 0.5)
    ax_hist2.bar( my_hist_x_test, my_hist_y_test, width, label='Histogram', color=color2)
    ax_hist2.tick_params(axis='both', which='major', labelsize=9)
    ax_hist2.set_xticks(my_tics_hist)
    ax_hist2.set_xlim([-33, 33])
    ax_hist2.set_xlabel('Test errors (meV)  ', fontsize=10); ax_hist2.set_ylabel('Probability', fontsize=10)
    plt.xticks(rotation=90 )

    if (forecast_of_residuals): text='fromresid'
    else: text='NOTfromresid'
    pathout = dir_out + '/Plots/'+quantity+'_renorm_'+method+'_'+text+'.pdf'
    plt.savefig(pathout, format="pdf", bbox_inches="tight")

    plt.clf(); plt.cla(); plt.close('all'); del ax; del fig

    return

'''
# Example run:
##np.random.seed(0)
x_train  = 20 + 180*uniform.rvs(size=150)
y_train  = x_train + norm.rvs(size=150, loc=0, scale=15 )
x_test   = 20 + 180*uniform.rvs(size=35)
y_test   = x_test + norm.rvs(size=35, loc=0, scale=25 )
train_error  = norm.rvs(size=30, loc=-1, scale=10 )
test_error = norm.rvs(size=150, loc=1, scale=10  )
integrated_plot( x_train, y_train, x_test, y_test, train_error, test_error, "HOMO" )
'''
#-----------------------------------------------------------------------------------------------

def plot_histogram_orig(  fitted_time_series, plot_several_curves=False, params_several_curves=None ):
    ''' This function plots a histogram.
    '''


    # Define the bins
    len_to_remove = int(0.01 * len(fitted_time_series.ts_to_fit))
    lim_bins =  4*np.std( fitted_time_series.ts_to_fit[len_to_remove:-len_to_remove] )
    num_bins = 11
    bins = [-lim_bins  ] # OLD: [-lim_bins + loc_in ]; Si el loc se calculó mal, esto hace que el histograma tb se pinte mal.
    labels = []
    for i in range(1,num_bins+1):
        bins   += [ -lim_bins + i*(2*lim_bins)/num_bins ] #old:  [ -lim_bins + loc_in  + i*(2*lim_bins)/num_bins ]
        if (i%2 == 0):
          if (lim_bins >= 1):
            labels += [ "{:.4f}".format( -lim_bins  + (i-0.5)*(2*lim_bins)/num_bins  ) ] # old: labels += [ "{:.4f}".format( -lim_bins + loc_in  + (i-0.5)*(2*lim_bins)/num_bins  ) ]
          else:
            labels += ["{:.5f}".format(-lim_bins  + (i - 0.5) * (2 * lim_bins) / num_bins)] # old:  ["{:.5f}".format(-lim_bins + loc_in  + (i - 0.5) * (2 * lim_bins) / num_bins)]
        else:
            labels += [None]

    # Calculate values of the "float" histogram
    counts0, bins = np.histogram( fitted_time_series.ts_to_fit, bins )
    counts = np.zeros(len(counts0))
    for i in range(len(counts0)):
        counts[i] = counts0[i]/len(fitted_time_series.ts_to_fit)

    # Plot the histogram
    x = np.array( [ (bins[i]+bins[i+1])/2 for i in range(len(labels)) ] ) # old: np.arange(len(labels))

    width = 0.91*(bins[i+1]-bins[i])
    fig, ax = plt.subplots()
    ax.bar(x , counts, width, label='Histogram (normalized)', color="darkblue") #'mediumblue'

    # Plot fitting function:
    for fitting_function in list_sweep:
        #loc_plot = (len(labels)-1)/2
        #scale_plot = scale_in * len(labels) / (2*lim_bins)
        x_cont = x
        x_cont0 = np.array( [ x[0] + i*(x[-1]-x[0])/(12*num_bins-1) for i in range(12*num_bins) ]  ) #np.linspace(0, len(labels) - 1, num=num_bins*10  )
        name_plot = fitted_time_series.filename; name_plot = name_plot.replace(".csv","")

        loc_in, scale_in, skew_in, tail1_in, tail2_in = provide_curve_parameters(fitting_function, fitted_time_series, plot_several_curves, params_several_curves)
        if (None in [loc_in, scale_in]): continue
        if ( fitting_function =="norm"):
           y_cont0 = norm.pdf(x_cont0, loc=loc_in, scale=scale_in ) * (bins[i + 1] - bins[i] )
           y_cont  = norm.pdf(x_cont,  loc=loc_in, scale=scale_in ) * (bins[i + 1] - bins[i] )
           mycolor  = "red"
           ax.plot(x_cont0, y_cont0, color='black', lw=5.2)
           if not (plot_several_curves): ax.plot(x_cont, y_cont, '.', color='black', markersize=17)
           ax.plot(x_cont0, y_cont0, label='Fit to normal', color=mycolor, lw=3)
           pathout = f"{fitted_time_series.directory_output_plots}/histofit_{name_plot}_normal.pdf"
        elif (fitting_function == "nct"):
           y_cont0 = nct.pdf(x_cont0, loc=loc_in, scale=scale_in, nc=skew_in, df=tail1_in)*(bins[i+1]-bins[i])
           y_cont  = nct.pdf( x_cont,  loc=loc_in, scale=scale_in, nc=skew_in, df=tail1_in)*(bins[i+1]-bins[i])
           mycolor = "lime"#'#00ced1'#
           ax.plot(x_cont0, y_cont0, color='black', lw=5.2)
           if not (plot_several_curves): ax.plot(x_cont, y_cont, '.', color='black', markersize=17)
           ax.plot(x_cont0, y_cont0, label='Fit to t-student (nct)', color=mycolor, lw=3)
           pathout = f"{fitted_time_series.directory_output_plots}/histofit_{name_plot}_nct.pdf"
        elif ( fitting_function =="levy_stable"):
           y_cont0 = levy_stable.pdf( x_cont0 , loc=loc_in, scale=scale_in,  beta=skew_in, alpha=tail1_in )*(bins[i+1]-bins[i])
           y_cont  = levy_stable.pdf( x_cont,   loc=loc_in, scale=scale_in,  beta=skew_in, alpha=tail1_in)*(bins[i+1]-bins[i])
           mycolor = '#ffa500'#'#daa520'
           ax.plot(x_cont0, y_cont0, color='black', lw=5.2)
           if not (plot_several_curves): ax.plot(x_cont, y_cont, '.', color='black', markersize=17)
           ax.plot(x_cont0, y_cont0, label='Fit to stable', color=mycolor, lw=3)
           pathout = f"{fitted_time_series.directory_output_plots}/histofit_{name_plot}_stable.pdf"
        elif ( fitting_function=="genhyperbolic"):
           y_cont0 = genhyperbolic.pdf(x_cont0, loc=loc_in, scale=scale_in,  b=skew_in, a=tail1_in, p=tail2_in )*(bins[i+1]-bins[i])
           y_cont  = genhyperbolic.pdf( x_cont, loc=loc_in, scale=scale_in,  b=skew_in, a=tail1_in, p=tail2_in )*(bins[i+1]-bins[i])
           mycolor = "#ff1493"#"#FF00FF"
           ax.plot(x_cont0, y_cont0, color='black', lw=5.2)
           if not (plot_several_curves): ax.plot(x_cont, y_cont, '.', color='black', markersize=17)
           if ( fitted_time_series.consider_p):
               ax.plot(x_cont0, y_cont0, label='Fit to hyperbolic (with p)', color=mycolor, lw=3)
               pathout = f"{fitted_time_series.directory_output_plots}/histofit_{name_plot}_hyperbolicwithp.pdf"
           else:
               ax.plot(x_cont0, y_cont0, label='Fit to hyperbolic', color=mycolor, lw=3)
               pathout = f"{fitted_time_series.directory_output_plots}/histofit_{name_plot}_hyperbolicwop.pdf"
        if not (plot_several_curves):
            ax.plot(x_cont, y_cont, '.', label='Fit to ' + str( fitting_function),  color=mycolor,  markersize=14)

        del x_cont; del x_cont0; del y_cont0; del y_cont; del scale_in; del skew_in; del tail1_in; del tail2_in; del mycolor;

    # Add labels to the plot:
    if not (plot_several_curves):
        plt.ylim([0,  max(counts) *1.3 ] )
    else:
        plt.ylim([0, max(counts) * 1.5])
    ax.set_ylabel('Probability',fontsize=15)
    ax.set_xlabel('Spread residual',fontsize=14)
    ax.yaxis.set_label_coords(-0.12, 0.5)
    mytitle = str(fitted_time_series.filename)
    mytitle = mytitle.replace("Yield","Price"); mytitle = mytitle.replace("YIELD","Price"); mytitle = mytitle.replace("yield","Price"); # We change this because in the class definition we have calculated returns from prices, not from yields.
    mytitle = rewrite_title( mytitle )
    ax.set_title( mytitle ,fontsize=16 )
    ax.set_xticks(x)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.tick_params(axis='both', which='minor', labelsize=10)
    plt.xticks(rotation=45, ha='right')
    #ax.set_xticklabels(labels)
    handles, labels = plt.gca().get_legend_handles_labels()
    if (loc_in != None):
        if not (plot_several_curves):
            order = [2, 0]
            plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order], prop={'size': 11})
        else:
            order = [4,0,1,2,3]
            plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order], prop={'size': 11}, ncol=2, labelspacing=0.4)

    #plt.legend(loc='upper right' )
    fig.tight_layout()

    del len_to_remove; del lim_bins; del num_bins; del labels; del x; del width; del counts; del counts0

    #plt.show()

    # Save the plot:
    if (plot_several_curves): pathout=pathout.replace(".pdf","_all.pdf")
    plt.savefig( pathout, format="pdf", bbox_inches="tight")
    plt.clf(); plt.cla(); plt.close('all')
    del ax; del fig

    return

#-----------------------------------------------------------------------------------------------


def compile_files_for_box_plots(  prefix  ):

    from os import path

    directory = module_io.execution_directory + "/Results/"
    list_suffix = [ "",  "ElecEigvalKurtosis","inv((LUMO+2)-HOMO)" ,"PhFreqKurtosis",  "Areavolumeratio",  "bond_length_1st_highest_value"  , "Mayer_bond_order_1st_highest_value", "1st_highest_value_GJ_bond_order", "1st_highest_value_NM1_bond_order", "1st_highest_value_NM2_bond_order", "1st_highest_value_NM3_bond_order" ]
    dict_labels = { "":"Elec. str. selected","1st_highest_value_GJ_bond_order":"Bond orders - GJ", "1st_highest_value_NM1_bond_order":"Bond orders - NM1", "1st_highest_value_NM2_bond_order":"Bond orders - NM2", "1st_highest_value_NM3_bond_order":"Bond orders - NM3", "Areavolumeratio":"Geometric", "ElecEigvalKurtosis":"Elec. struct. (few)", "Mayer_bond_order_1st_highest_value":"Bond orders - Mayer", "PhFreqKurtosis":"Phonons", "bond_length_1st_highest_value":"Bond lengths", "inv((LUMO+2)-HOMO)":"Elec. str. (all)" }


    df_all = pd.DataFrame()

    for suffix in list_suffix:

        if (suffix==""): filepath = directory + prefix +".csv"
        else:            filepath = directory + prefix + "_" + suffix + ".csv"
        if not (path.exists(filepath)): print("\nWARNING:",filepath,"not found; omitting it for plotting.")
        else: print(" Now reading",filepath)
        df0 = pd.read_csv(filepath,header=0,usecols=["y_err"])
        #df0["y_err"] = -df0["y_err"]
        df0.columns = [dict_labels[suffix]]
        df_all = pd.concat([df_all, df0], axis=1)
        #print(suffix,"\n",df0,"weams\n",df_all,"\n\n\n")


    df_all = df_all.dropna()


    del df0; del prefix; del suffix; del dict_labels; del directory; del list_suffix

    return df_all

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def box_plot_quantiles( prefix ):

    df_to_plot = compile_files_for_box_plots( prefix )

    #x = pd.read_csv(filepath,header=0,usecols=[colname])

    fig, ax = plt.subplots()

    my_title = "Errors in the "
    for type in ["HOMO","LUMO","gap"]:
        if (type in prefix):
            my_title += type
    my_title += " renormalization - "
    my_title = my_title.replace("gap","Gap")
    if ("RF" in prefix): my_title+="Random forests"
    elif (("kNN" in prefix) or ("KNN" in prefix)):my_title += "k-nearest neighbors"
    elif ("NN" in prefix): my_title += "Neural networks"


    if ("HOMO" in prefix):
        #plt.ylim([-80,60])
        #plt.yticks([-80, -60, -40, -20, 0, 20, 40, 60])
        plt.ylim([-50, 50])
        plt.yticks([-50, -25, 0, 25, 50])
    elif ("LUMO" in prefix):
        plt.ylim([-50,50])
        plt.yticks([-50,-25,0,25,50])
    elif ("gap" in prefix):
        plt.ylim([-100,100])
        plt.yticks([-50, -25, 0, 25, 50])



    ax.set_title( my_title, fontsize=17 )
    meanlineprops = dict(linestyle=':', linewidth=1.5, color='purple')
    medianlineprops = dict(linestyle='-', linewidth=1.5, color='b')

    #ax = df_to_plot [ list(df_to_plot ) ].plot(kind='box', title='boxplot')
    ax.boxplot(df_to_plot [ list(df_to_plot ) ], medianprops=medianlineprops, meanprops=meanlineprops, meanline=True, showmeans=True, showfliers=False,whis=(5, 95))

    plt.xticks([1, 2, 3,4,5,6,7,8,9,10,11], ["Elec. str. selected","Elec. struct. (few)","Elec. str. (all)" ,"Phonons","Geometric","Bond lengths","Bond orders - Mayer","Bond orders - GJ","Bond orders - NM1","Bond orders - NM2","Bond orders - NM3"])
    plt.xticks(rotation=45, ha='right')
    plt.ylabel("Error (meV)", fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.tick_params(axis='both', which='minor', labelsize=12)

    plot_path = module_io.execution_directory + "/Plots/" +"/boxplot_"+prefix + ".pdf"
    plt.savefig(plot_path, format="pdf", bbox_inches="tight")
    #plt.show()
    plt.clf(); plt.cla(); plt.close('all'); del ax; del fig

    return

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def plot_feature_importance_old( ml_variable, feature_names, quantity, method ):
    ''''''


    importances = ml_variable.feature_importances_
    std = np.std([tree.feature_importances_ for tree in ml_variable.estimators_], axis=0)
    ml_importances = pd.Series(importances, index=feature_names)

    print("Weamos\n",ml_importances)

    fig, ax = plt.subplots()
    ml_importances.plot.bar(yerr=std, ax=ax)
    ax.set_title("Feature importances using MDI")
    ax.set_ylabel("Mean decrease in impurity")
    fig.tight_layout()
    plt.show()

    plot_path = module_io.execution_directory + "/Plots/" + "/feature_importance_" +  quantity + "_" +  method + ".pdf"
    plt.savefig(plot_path, format="pdf", bbox_inches="tight")

    plt.clf(); plt.cla(); plt.close('all'); del ax; del fig

    return

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def plot_feature_importance( features_weights, quantity, method, forecast_of_residuals  ):
    ''''''

    dict_methods = { "RF":"Random Forests","NN":"Neural Networks","KNN":"k-Nearest Neighbors" }

    my_subtitle = quantity+" renormalization - "+ dict_methods[method]
    filename = 'FI_MDI_' + quantity + "_" + method
    if (forecast_of_residuals):
        my_subtitle += " (resid.)"
        filename += "_resid"
    else:
        filename += "_woresid"

    ax = plt.subplot(111)
    plt.suptitle("Feature Importance - MDI", fontsize=17)
    plt.title(my_subtitle, fontsize=14)
    plt.ylabel("Mean decrease in impurity", fontsize=16)
    plt.xticks(rotation=45,ha='right', fontsize=12)
    plt.yticks([0.05, 0.10, 0.15, 0.20, 0.25, 0.3, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70],["", 0.10, "", 0.20, "", 0.3, "", 0.40, "", 0.50, "", 0.60, "", 0.70], fontsize=12)
    ax.bar(features_weights.keys(), features_weights.values(), width=0.7, align='center')  # color='blue',
    plt.ylim([0, 0.05 + 0.05 * int(max(features_weights.values())/0.05) ]);

    plot_path = module_io.execution_directory + "/Plots/" + "/" + filename + ".pdf"
    plt.savefig(plot_path, format="pdf", bbox_inches="tight")

    plt.savefig("/Users/pgr/Desktop/"+filename+"-datasubset.pdf", format="pdf", bbox_inches="tight")

    plt.clf(); plt.cla(); plt.close('all'); del ax


    return

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def plot_feature_importance_alldata_MDI( my_configuration, field_out_to_read, regressor_fields_lr, regressor_fields_ml, param_sw, N_essais, forecast_of_residuals  ):
    '''This function draws feature importance using Mean Decreasing in Impurity. To this end, all the
    data are used, this is here ther is no separation between training and test data. If you want to plot
    doing such separation, use set the corresponding flag (plot_FI_MDI_separating_train) to True in module_io.py.'''

    from sklearn.inspection import permutation_importance

    ''' 
    quantity = my_configuration.quantity_to_forecast
    method = my_configuration.ml_method
    features_weights = {'avgoccupied(eV-1)': 0.24245583897385756, 'Gap_PBE': 0.02026678190740172, 'HOMO-(HOMO-1)': 0.08475717436917594, 'HOMO-(HOMO-2)': 0.08112798549919549, 'HOMO-(HOMO-3)': 0.06676476406478968, 'HOMO-(HOMO-4)': 0.056106953091398994, 'HOMO-(HOMO-5)': 0.036441726146585716, 'inv(HOMO-(HOMO-1))': 0.06408082370818045, 'inv(HOMO-(HOMO-2))': 0.08843977177163773, 'inv(HOMO-(HOMO-3))': 0.08693205270128751, 'inv(HOMO-(HOMO-4))': 0.06881151928556523, 'inv(HOMO-(HOMO-5))': 0.038230577808364294, 'inv(LUMO-(HOMO-1))': 0.010418599568986164, 'inv(LUMO-(HOMO-2))': 0.00987107325264252, 'inv((LUMO+1)-HOMO)': 0.02194201071224471, 'inv((LUMO+2)-HOMO)': 0.023352347138686357}
    make_plot_feature_importances( features_weights, "Feature Importance - MDI", forecast_of_residuals,  quantity, method )
    features_weights_from_permut = {'avgoccupied(eV-1)': 0.23002306708707498, 'Gap_PBE': 0.008897227388372411, 'HOMO-(HOMO-1)': 0.04006272564517631, 'HOMO-(HOMO-2)': 0.0367421719736308, 'HOMO-(HOMO-3)': 0.023609985151293134, 'HOMO-(HOMO-4)': 0.01900128545570482, 'HOMO-(HOMO-5)': 0.011050236757638504, 'inv(HOMO-(HOMO-1))': 0.042094701186923045, 'inv(HOMO-(HOMO-2))': 0.032612788902985586, 'inv(HOMO-(HOMO-3))': 0.03386814283240922, 'inv(HOMO-(HOMO-4))': 0.02462444304309945, 'inv(HOMO-(HOMO-5))': 0.011464583338866652, 'inv(LUMO-(HOMO-1))': 0.0048202992036918794, 'inv(LUMO-(HOMO-2))': 0.004874443718317562, 'inv((LUMO+1)-HOMO)': 0.012287723666235366, 'inv((LUMO+2)-HOMO)': 0.011918892703338008}
    make_plot_feature_importances( features_weights_from_permut, "Feature Importance - Permutations",  forecast_of_residuals, quantity, method)
    exit(0)
    '''

    # PART 1: Calculation of mean decreases in impurity.

    quantity      = my_configuration.quantity_to_forecast
    method        = my_configuration.ml_method
    li_regressors = my_configuration.regressor_fields_ml

    if (forecast_of_residuals):
        from module_training import ols_multiregression
        output_field = str(my_configuration.field_out_to_read + "_residual")  # str(quantity_to_forecast + "_residual")
        ols_multiregression(field_out_to_read, regressor_fields_lr, regressor_fields_ml, module_io.filepath_raw_input_data, module_io.filepath_raw_input_data, module_io.filepath_raw_input_data, 0, my_configuration.different_regressors, False )
        filepath_data = module_io.filepath_raw_input_data.replace(".csv", "-with_residuals.csv")
    else:
        output_field = str(my_configuration.field_out_to_read)
        filepath_data = module_io.filepath_raw_input_data

    listcols = ["Name"] + [output_field]
    if (forecast_of_residuals): listcols.append(output_field.replace("_residual", ""))
    listcols_train = listcols.copy().append(output_field.replace("_residual", "") + "_forecasted")
    df_output_all     = pd.read_csv( filepath_data, header=0, usecols=listcols_train)
    df_regressors_all = pd.read_csv( filepath_data, header=0, usecols=["Name"] + li_regressors)
    df_output_all.set_index("Name", inplace=True)
    df_regressors_all.set_index("Name", inplace=True)
    m = len(df_output_all)
    Nregressors = len(li_regressors)

    X_all = np.zeros((m,Nregressors))  # We declare and initialize the independent variables (X, regressors) for the "training" (i.e. calculation of parameters of the regression).
    y_all = np.zeros(m)                # We declare and initialize the dependent variable (y) for the training.

    for i in range(m):
            y_all[i] = float( df_output_all.iloc[i][output_field] )
            for j in range(len(li_regressors)):  # We avoid the 1st, 2nd and last fields because they are 'Name' (molecule label) and the output_variable (e.g. 'HOMO ren (output)_residual')
                X_all[i,j] = float( df_regressors_all.iloc[i][str(li_regressors[j])] )

    # Actual ML training:

    if ((method == "NN") or (method == "Neural networks")):
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

    if ((method == "NN") or (method == "Neural networks")):
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
        print("\nERROR: Unknown ML method",method,"\n"); exit(1)

    features_weights             = dict(zip(my_configuration.regressor_fields_ml, [0.0] * (len(my_configuration.regressor_fields_ml))))
    features_weights_from_permut =  features_weights.copy()

    for i in range(N_essais):

        ml_variable.fit(X_all, y_all)
        importances = ml_variable.feature_importances_
        ml_importances = pd.Series(importances, index=li_regressors)
        importances_from_permut    = permutation_importance( ml_variable, X_all, y_all, n_repeats=10 )
        ml_importances_from_permut = pd.Series( importances_from_permut.importances_mean , index=li_regressors )
        for regressor_name in li_regressors:
            features_weights[regressor_name]             += ml_importances[regressor_name]
            features_weights_from_permut[regressor_name] += ml_importances_from_permut[regressor_name]

    for regressor_name in my_configuration.regressor_fields_ml:
        features_weights[regressor_name]             /= N_essais
        features_weights_from_permut[regressor_name] /= N_essais

    print("Features weights:\n",features_weights)
    print("Features weights from permutation:\n", features_weights_from_permut, "\n")


    # PART 2: Actual plotting

    make_plot_feature_importances( features_weights, "Feature Importance - MDI", forecast_of_residuals,  quantity, method )
    make_plot_feature_importances( features_weights_from_permut, "Feature Importance - Permutations",  forecast_of_residuals, quantity, method)

    return

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def make_plot_feature_importances( features_weights, my_suptitle, forecast_of_residuals, quantity, method ):

    #dict_methods = {"RF": "Random Forests", "NN": "Neural Networks", "KNN": "k-Nearest Neighbors"}

    my_subtitle = quantity + " renormalization"
    if ("MDI" in my_suptitle):
        filename = 'FI_MDI_' + quantity + "_" + method
        my_title = "Mean decrease in impurity"
    else:
        filename = 'FI_permut_'+ quantity + "_" + method
        my_title = "Mean accuracy decrease"

    if (forecast_of_residuals):
        my_subtitle += " (resid.)"
        filename += "_resid"
    else:
        filename += "_woresid"

    fig, ax = plt.subplots(figsize=(max(6.4,0.1*6.4*(len(features_weights.keys()))),4.8))
    plt.suptitle(my_suptitle, fontsize=17)
    plt.title(my_subtitle, fontsize=14)
    plt.ylabel(my_title, fontsize=16)
    plt.xticks(rotation=45,ha='right', fontsize=12)
    ax.bar(features_weights.keys(), features_weights.values(), width=0.7, align='center')  # color='blue',
    my_ylim =  0.05 * int(max(features_weights.values())/0.05) + 0.05
    if ((my_ylim-max(features_weights.values() ) ) < 0.02):
        my_ylim += 0.05
    if (my_ylim < 0.40001):
        plt.yticks([0.05, 0.10, 0.15, 0.20, 0.25, 0.3, 0.35, 0.40], [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40],fontsize=12)
    else:
        plt.yticks([0.05, 0.10, 0.15, 0.20, 0.25, 0.3, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0],["", 0.1, "", 0.2, "", 0.3, "", 0.4, "", 0.5, "", 0.6, "", 0.7,"", 0.8, "", 0.9, "", 1.0], fontsize=12)
    plt.ylim([0, my_ylim]);
    #f = plt.figure()
    #f.set_figwidth()
    #f.set_figheight(4.8)
    #print("La length es ",len(features_weights.keys()))

    plot_path = module_io.execution_directory + "/Plots/" + "/" + filename + ".pdf"
    plt.savefig(plot_path, format="pdf", bbox_inches="tight")

    plt.savefig("/Users/pgr/Desktop/"+filename+".pdf", format="pdf", bbox_inches="tight")#xxx

    plt.clf(); plt.cla(); plt.close('all'); del ax

    return

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
