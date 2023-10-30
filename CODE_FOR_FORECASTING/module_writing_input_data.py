''' This module contains the functions which use Quantum Espresso's output files
and generate the data which must be used as input for the code for forecasting, this is the numerical values of the
second row of the file input_regressors_molecule_to_forecast.csv.
'''

import pandas as pd
from os import path, makedirs
from numpy import mean,  var
from scipy.stats import skew, kurtosis

#-------------------------------------------------------------------------------------------------------------------------

def reinv( x ):
    if (abs(x)<0.00000001): return 0
    else: return 1/x

#-------------------------------------------------------------------------------------------------------------------------

def extract_data_from_fullerenes( name_fullerene_to_forecast="my_fullerene" ):

        print("\n *** NOW GENERATING THE FILE WHICH CONTAINS THE REGRESSORS (input_regressors_molecule_to_forecast.csv) ***")
        print("To this end the eigenvalues from the eigenval.xml file are used.")

        for Nmaxunoccavg in [20]:
          print(  "   --> Average for unoccupied with up to ",Nmaxunoccavg ," states.")
          df_out = pd.DataFrame()
          print("   Now analysing", name_fullerene_to_forecast )
          foldername = "Quantum_Espresso_output/"
          gap_PBE, avgoccupied, avgempty, countocc, countempty, difHOMO1, difHOMO2, difHOMO3, difHOMO4, difHOMO5, difHOMO6, difLUMO1, difLUMO2, difLUMO3, difLUMO4, difLUMO5, difLUMO6, difLH1, difLH2, difL1H, difL2H, mymean, myvar, myskew, mykurt = find_avgseigenvalues( foldername, Nmaxunoccavg )
          #print("   Now analysing", label_fullerene," ; avg. occupied:", avgoccupied )
          df_onerow = pd.DataFrame( data={ "molecule":[name_fullerene_to_forecast],"Gap_PBE":[gap_PBE], "avgoccupied(eV-1)":[avgoccupied], "avgempty(eV-1)":[avgempty], "countocc":[countocc], "countempty":[countempty],
                                    "HOMO-(HOMO-1)":[difHOMO1], "HOMO-(HOMO-2)":[difHOMO2], "HOMO-(HOMO-3)":[difHOMO3], "HOMO-(HOMO-4)":[difHOMO4], "HOMO-(HOMO-5)":[difHOMO5], "HOMO-(HOMO-6)":[difHOMO6],
                                    "inv(HOMO-(HOMO-1))":[reinv(difHOMO1)], "inv(HOMO-(HOMO-2))":[reinv(difHOMO2)], "inv(HOMO-(HOMO-3))":[reinv(difHOMO3)], "inv(HOMO-(HOMO-4))":[reinv(difHOMO4)], "inv(HOMO-(HOMO-5))":[reinv(difHOMO5)], "inv(HOMO-(HOMO-6))":[reinv(difHOMO6)],
                                    "(LUMO+1)-LUMO":[difLUMO1],"(LUMO+2)-LUMO":[difLUMO2], "(LUMO+3)-LUMO":[difLUMO3], "(LUMO+4)-LUMO":[difLUMO4], "(LUMO+5)-LUMO":[difLUMO5], "(LUMO+6)-LUMO":[difLUMO6],
                                    "inv((LUMO+1)-LUMO)":[reinv(difLUMO1)],"inv((LUMO+2)-LUMO)":[reinv(difLUMO2)], "inv((LUMO+3)-LUMO)":[reinv(difLUMO3)], "inv((LUMO+4)-LUMO)":[reinv(difLUMO4)], "inv((LUMO+5)-LUMO)":[reinv(difLUMO5)], "inv((LUMO+6)-LUMO)":[reinv(difLUMO6)],
                                    "inv(LUMO-(HOMO-1))": [reinv(difLH1)], "inv(LUMO-(HOMO-2))":[reinv(difLH2)], "inv((LUMO+1)-HOMO)":[reinv(difL1H)], "inv((LUMO+2)-HOMO)":[reinv(difL2H)],
                                    "ElecEigvalMean": [mymean],"ElecEigvalVariance": [myvar],"ElecEigvalSkewness": [myskew],"ElecEigvalKurtosis":[mykurt]       } )
          df_out = pd.concat( [df_out, df_onerow], axis=0, sort=False)

          df_out.to_csv( f"input_regressors_molecule_to_forecast.csv",index=False )

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

         difLH1 =  ( LUMO - df3.iloc[i-1]["eigenvalues"] )*hartree_to_eV #(LUMO-(HOMO-1))"
         difLH2 =  ( LUMO - df3.iloc[i-2]["eigenvalues"] )*hartree_to_eV #LUMO-(HOMO-2)
         difL1H =  ( df3.iloc[i + 2]["eigenvalues"] - HOMO )*hartree_to_eV #(LUMO+1)-HOMO
         difL2H =  ( df3.iloc[i + 3]["eigenvalues"] - HOMO )*hartree_to_eV #(LUMO+2)-HOMO


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


         mymean = ( mean(df3[0:iHOMO]["eigenvalues"] -HOMO ))*hartree_to_eV
         myvar = (var(df3[0:iHOMO ]["eigenvalues"]  ))*hartree_to_eV*hartree_to_eV
         myskew = (skew(df3[0:iHOMO]["eigenvalues"]))
         mykurt = (3+kurtosis(df3[0:iHOMO]["eigenvalues"]))


         return ((LUMO-HOMO)*hartree_to_eV*1000), avgoccupied, avgempty, countocc, countempty, difHOMO1, difHOMO2, difHOMO3, difHOMO4, difHOMO5, difHOMO6, difLUMO1, difLUMO2, difLUMO3, difLUMO4, difLUMO5, difLUMO6, difLH1, difLH2, difL1H, difL2H, mymean, myvar, myskew, mykurt

#------------------------------------------------------------------------------------------

extract_data_from_fullerenes()