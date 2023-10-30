
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
