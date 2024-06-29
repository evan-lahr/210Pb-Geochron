################################################################################
###                                 LEADTOOLS.py                             ###
#         3 functions: spe_to_counts, DET_MATCH_SUM, counts_to_activity        #
#           Function descriptions are given directly above the code.           #
###                                Evan Lahr 2020                            ###
################################################################################

################################################################################
###                                spe_to_counts                             ###
#                     this script performs the following:                      #
#        1. Open a folder of SPE files and compute/write counting info to df       
#        2. Match a csv of lab measurements to the df from step 1. 
#        3. Write a csv with the columns specified below (see "OUTPUT")
#
#  ------INPUT #1:  SPE FOLDER--------      --------------OUTPUT-----------------
#  . type: string path               .      . Z_midpt (cm)	                    .
#  . content: folder of .SPE files   .      . Δt_in_counting (sec)	            .
#  -----------------------------------      . ΔZ (cm)	M_pan (g)	            .
#                   \\                      . M_WetSed+Pan (g)	                .
#                    \\                     . M_DrySed+Pan (g)	                .
#                    ====================>> . siltclay (volfrac)	            .
#                    //                     . M_WetChemSed (g)	                .
#                   //                      . Plating_StartDate (DD/MM/YYYY)    .
#  ------INPUT #2:  LAB DATA CSV------      . Plating_StartTime (HH:MM:SS)	    .
#  . type: string path               .      . detID	                            .
#  . content: csv                    .      . 209Po_decays (counts)	            .
#  . columns:                        .      . 210Po_decays (counts)	            .
#  .   CoreID                        .      . Counting_StartDate+Time	        .
#  .   Z_upper (cm)                  .      . Counting_StartTime	            .
#  .   Z_lower (cm)                  .      . Counting_StartDate                .
#  .   Plating_StartDate (DD/MM/YYYY).      -------------------------------------
#  .   Plating_StartTime (HH:MM:SS)  .
#  .   M_pan (g)                     .
#  .   M_WetSed+Pan (g)              .
#  .   M_DrySed+Pan (g)              .
#  .   M_WetChemSed (g)              .
#  .   siltclay (volfrac)            .
#  -----------------------------------
#
###                                                                          ###
################################################################################



def spe_to_counts(SPEs_path, labsheet_path, fout, **PlotSPEs):
    print("|------------------------  spe_to_counts STARTED  ----------------------|")

    # import modules
    import pandas as pd
    import numpy as np
    import glob

    # create list of all the spe files to open using the input folder path
    files = glob.glob(SPEs_path)
    df_files = pd.DataFrame(files)

    # print statement for verification
    print(f"||    Reading {len(files)} spe files at path:            {SPEs_path}")

    # read and compute all relevant info from each spe file, write into dfs
    # create empty dfs to be added to in for loop
    spe = pd.DataFrame()
    speDet = pd.DataFrame()
    speDate = pd.DataFrame()
    speCounts = pd.DataFrame()
    spectra_sum = pd.DataFrame()
    counts = pd.DataFrame()
    # loop through each individual spe file and read characters into dataframe
    for i in range(len(files)):
        spe_raw = pd.read_csv(files[i])
        # read relevant regions of the spe file
        spe_new = pd.DataFrame(spe_raw[11:2059].astype(int))  # raw counts data
        speDet_new = spe_raw[2:3].astype(str)  # name of detector ID
        speDate_new = spe_raw[6:7].astype(str)  # date of counting
        speCounts_new = spe_raw[8:9].astype(str)  # total number of counts
        # add new values to existing dataframe
        spe = pd.concat([spe, spe_new], axis=1, ignore_index=True)
        speDet = pd.concat([speDet, speDet_new], axis=1, ignore_index=True)
        speDate = pd.concat([speDate, speDate_new], axis=1, ignore_index=True)
        speCounts = pd.concat(
            [speCounts, speCounts_new], axis=1, ignore_index=True
        )
        # Match spes to detectors, sum α-decays with the function "det_match_sum"
        spectra_sum_new = pd.DataFrame(det_match_sum(spe[i], speDet[i][2],files[i], PlotSPEs['PlotSPEs'])).T
        spectra_sum = pd.concat(
            [spectra_sum, spectra_sum_new], axis=0, ignore_index=True
        )

    # create df named 'counts' to store final values, begin adding computed columns
    counts = pd.DataFrame()
    # the depth midpoint at section i (temp, will be overwritten)
    counts["Z_midpt (cm)"] = pd.DataFrame(files)
    # total elapsed counting time in seconds
    counts["Δt_in_counting (sec)"] = speCounts.T

    # Compute a few more values and concat
    midpt = np.zeros(len(files))  # midpoint of the section interval (cm bsf)
    totcts = np.zeros(len(files))  # total number of counting seconds
    depInt = np.zeros(
        len(files)
    )  # the vertical thickness of the section analyzed
    for i in range(len(files)):
        midpt[i] = (
            int(counts["Z_midpt (cm)"][i][22:25])
            + int(counts["Z_midpt (cm)"][i][26:29])
        ) / 2
        totcts[i] = int(counts["Δt_in_counting (sec)"][i].split(" ")[0])
        depInt[i] = int(counts["Z_midpt (cm)"][i][26:29]) - int(
            counts["Z_midpt (cm)"][i][22:25]
        )
    
    # add computed values to df
    # section i depth midpoint
    counts["Z_midpt (cm)"] = midpt
    # section i vertical thickness
    counts["ΔZ (cm)"] = depInt
    # the detector bin ID that the sample was counted in
    counts["detID"] = spectra_sum[:][0]
    # the elapsed time between the start and end of counting
    counts["Δt_in_counting (sec)"] = totcts
    # the total number of 209Po α-decays detected
    counts["209Po_decays (counts)"] = spectra_sum[:][1]
    # the total number of 210Po α-decays detected
    counts["210Po_decays (counts)"] = spectra_sum[:][2]
    # the date and time of α-counting
    counts["Counting_StartDate+Time"] = speDate.T
    # spe format is hh:mm:ss
    counts["Counting_StartTime"] = (
        counts["Counting_StartDate+Time"].astype(str).str[11:])
    # spe format is MM/DD/YYYY
    counts["Counting_StartDate"] = (
        counts["Counting_StartDate+Time"].astype(str).str[:10] )
    # sort the dataframe by section depth
    counts = counts.sort_values(by=["Z_midpt (cm)"], ignore_index=True)
    
    
    
    # read in weight, grain size, and time data obtained in the lab
    labsheet = pd.read_csv(labsheet_path, header=0).sort_values(by=["Z_upper (cm)"], ignore_index=True)
    print( f"||    Read {len(labsheet)} rows of lab csv data at path: {labsheet_path}")    
    # mass of the pan used to dry sediment
    counts["M_pan (g)"] = labsheet["M_pan (g)"]
    # mass of the pan plus wet sediment
    counts["M_WetSed+Pan (g)"] = labsheet["M_WetSed+Pan (g)"]
    # mass of the pan plus dried sediment
    counts["M_DrySed+Pan (g)"] = labsheet["M_DrySed+Pan (g)"]
    # the volume fraction of mud expressed as a value between 0 (0%) and 1 (100%)
    counts["siltclay (volfrac)"] = labsheet["siltclay (volfrac)"]
    # the mass of the crushed sediment used in 210Pb analysis.
    counts["M_WetChemSed (g)"] = labsheet["M_WetChemSed (g)"]
    # the date of polonium plating onto the planchet
    counts["Plating_StartDate (DD/MM/YYYY)"] = labsheet[
        "Plating_StartDate (DD/MM/YYYY)"
    ]
    # the time of polonium plating onto the planchet
    counts["Plating_StartTime (HH:MM:SS)"] = labsheet[
        "Plating_StartTime (HH:MM:SS)"
    ]
    

    print(f"||    Writing data to csv at path:             {fout}")
    counts.to_csv(f"{fout}", index=False)
    print(
        f"|-------------------------  SPE_READER FINISHED  -----------------------|"
    )
    print(" ")
    return counts


################################################################################
###                              DET_MATCH_SUM                               ###
#  Uses header info from an SPE file to match it to a particular detector.     #
#  Sums all alpha decays within a detector-specific energy window     
#                                                                            
#      INPUTS  : "COUNTS": SPECTRAL DATA FROM A SINGLE SPE FILE
#                    "NAME"  : DETECTOR NAME AS OBTAINED VIA speName
#      PERFORMS:  MATCHES SPECTRAL DATA TO A SPECIFIC DETECTOR, AND
#                     SUMS DATA FOR 209Po, 210Po FOR EACH UNIQUE DETECTOR
#      OUTPUTS :  SUMMED 209Po, 210Po DATA FOR EACH COLUMN
###                                                                          ###
################################################################################


def det_match_sum(counts,name,sampleID,PlotSPEs):
    import numpy as np
    import matplotlib.pyplot as plt

    #SET ACTIVE CHANNEL BOUNDS FOR 209Po & 210Po
    #format is [209Po lower, 209Po upper, 210Po lower, 210Po upper]
    det1 = [608, 789, 789, 970]
    det2 = [608, 789, 789, 970]
    det3 = [608, 789, 789, 970]
    det4 = [658, 839, 839, 1020]
    det5 = [608, 789, 789, 970]
    det6 = [608, 789, 789, 970]
    det7 = [608, 789, 789, 970]
    det8 = [608, 789, 789, 970]
    
    
    if name == "DET# 1":
        po209 = np.sum(counts[det1[0] : det1[1]])
        po210 = np.sum(counts[det1[2] : det1[3]])
        detID = "EnsembleInput1"
        if PlotSPEs == True:
            fig, ax = plt.subplots(figsize=(10,4))
            ax.fill_between(counts.index.values[det1[0] : det1[1]],0, counts[det1[0] : det1[1]],color='red', zorder=3)
            ax.fill_between(counts.index.values[det1[2] : det1[3]],0, counts[det1[2] : det1[3]],color='blue', zorder=4)
            ax.fill_between(counts.index.values[det1[0]-100 : det1[3]+100],0, counts[det1[0]-100 : det1[3]+100],color='0.7', zorder=1)
            ax.text(.025, .92, 'SPECTRUM INTEGRATION RANGES', transform=ax.transAxes, fontsize='medium', fontweight='bold')
            ax.text(.025, .85, f'Detector: {detID}', transform=ax.transAxes, fontsize='small')
            ax.text(.025, .78, f'File: {sampleID}', transform=ax.transAxes, fontsize='small')
            ax.text(.885, .922, f'[{det1[0]}:{det1[1]}]', transform=ax.transAxes, fontsize='medium', color='red',zorder=6)
            ax.text(.885, .86, f'[{det1[2]}:{det1[3]}]', transform=ax.transAxes, fontsize='medium', color='blue',zorder=6)
            ax.legend(['209Po', '210Po', 'not counted         '])
            ax.set_ylabel('counts')
            ax.set_xlabel('decay energy channels')
            
      

    elif name == "DET# 2":
        po209 = np.sum(counts[det2[0] : det2[1]])
        po210 = np.sum(counts[det2[2] : det2[3]])
        detID = "EnsembleInput2"
        if PlotSPEs == True:
            fig, ax = plt.subplots(figsize=(10,4))
            ax.fill_between(counts.index.values[det2[0] : det2[1]],0, counts[det2[0] : det2[1]],color='red', zorder=3)
            ax.fill_between(counts.index.values[det2[2] : det2[3]],0, counts[det2[2] : det2[3]],color='blue', zorder=4)
            ax.fill_between(counts.index.values[det2[0]-100 : det2[3]+100],0, counts[det2[0]-100 : det2[3]+100],color='0.7', zorder=1)
            ax.text(.025, .92, 'SPECTRUM INTEGRATION RANGES', transform=ax.transAxes, fontsize='medium', fontweight='bold')
            ax.text(.025, .85, f'Detector: {detID}', transform=ax.transAxes, fontsize='small')
            ax.text(.025, .78, f'File: {sampleID}', transform=ax.transAxes, fontsize='small')
            ax.text(.885, .922, f'[{det2[0]}:{det2[1]}]', transform=ax.transAxes, fontsize='medium', color='red',zorder=6)
            ax.text(.885, .86, f'[{det2[2]}:{det2[3]}]', transform=ax.transAxes, fontsize='medium', color='blue',zorder=6)
            ax.legend(['209Po', '210Po', 'not counted         '])
            ax.set_ylabel('counts')
            ax.set_xlabel('decay energy channels')
        
        
    elif name == "DET# 3":
        po209 = np.sum(counts[det3[0] : det3[1]])
        po210 = np.sum(counts[det3[2] : det3[3]])
        detID = "EnsembleInput3"
        if PlotSPEs == True:
            fig, ax = plt.subplots(figsize=(10,4))
            ax.fill_between(counts.index.values[det3[0] : det3[1]],0, counts[det3[0] : det3[1]],color='red', zorder=3)
            ax.fill_between(counts.index.values[det3[2] : det3[3]],0, counts[det3[2] : det3[3]],color='blue', zorder=4)
            ax.fill_between(counts.index.values[det3[0]-100 : det3[3]+100],0, counts[det3[0]-100 : det3[3]+100],color='0.7', zorder=1)
            ax.text(.025, .92, 'SPECTRUM INTEGRATION RANGES', transform=ax.transAxes, fontsize='medium', fontweight='bold')
            ax.text(.025, .85, f'Detector: {detID}', transform=ax.transAxes, fontsize='small')
            ax.text(.025, .78, f'File: {sampleID}', transform=ax.transAxes, fontsize='small')
            ax.text(.885, .922, f'[{det3[0]}:{det3[1]}]', transform=ax.transAxes, fontsize='medium', color='red',zorder=6)
            ax.text(.885, .86, f'[{det3[2]}:{det3[3]}]', transform=ax.transAxes, fontsize='medium', color='blue',zorder=6)
            ax.legend(['209Po', '210Po', 'not counted         '])
            ax.set_ylabel('counts')
            ax.set_xlabel('decay energy channels')
        
        
    elif name == "DET# 4":
        po209 = np.sum(counts[det4[0] : det4[1]])
        po210 = np.sum(counts[det4[2] : det4[3]])
        detID = "EnsembleInput4"
        if PlotSPEs == True:
            fig, ax = plt.subplots(figsize=(10,4))
            ax.fill_between(counts.index.values[det4[0] : det4[1]],0, counts[det4[0] : det4[1]],color='red', zorder=3)
            ax.fill_between(counts.index.values[det4[2] : det4[3]],0, counts[det4[2] : det4[3]],color='blue', zorder=4)
            ax.fill_between(counts.index.values[det4[0]-100 : det4[3]+100],0, counts[det4[0]-100 : det4[3]+100],color='0.7', zorder=1)
            ax.text(.025, .92, 'SPECTRUM INTEGRATION RANGES', transform=ax.transAxes, fontsize='medium', fontweight='bold')
            ax.text(.025, .85, f'Detector: {detID}', transform=ax.transAxes, fontsize='small')
            ax.text(.025, .78, f'File: {sampleID}', transform=ax.transAxes, fontsize='small')
            ax.text(.885, .922, f'[{det4[0]}:{det4[1]}]', transform=ax.transAxes, fontsize='medium', color='red',zorder=6)
            ax.text(.885, .86, f'[{det4[2]}:{det4[3]}]', transform=ax.transAxes, fontsize='medium', color='blue',zorder=6)
            ax.legend(['209Po', '210Po', 'not counted         '])
            ax.set_ylabel('counts')
            ax.set_xlabel('decay energy channels')
        
        
    elif name == "DET# 5":
        po209 = np.sum(counts[det5[0] : det5[1]])
        po210 = np.sum(counts[det5[2] : det5[3]])
        detID = "EnsembleInput5"
        if PlotSPEs == True:
            fig, ax = plt.subplots(figsize=(10,4))
            ax.fill_between(counts.index.values[det5[0] : det5[1]],0, counts[det5[0] : det5[1]],color='red', zorder=3)
            ax.fill_between(counts.index.values[det5[2] : det5[3]],0, counts[det5[2] : det5[3]],color='blue', zorder=4)
            ax.fill_between(counts.index.values[det5[0]-100 : det5[3]+100],0, counts[det5[0]-100 : det5[3]+100],color='0.7', zorder=1)
            ax.text(.025, .92, 'SPECTRUM INTEGRATION RANGES', transform=ax.transAxes, fontsize='medium', fontweight='bold')
            ax.text(.025, .85, f'Detector: {detID}', transform=ax.transAxes, fontsize='small')
            ax.text(.025, .78, f'File: {sampleID}', transform=ax.transAxes, fontsize='small')
            ax.text(.885, .922, f'[{det5[0]}:{det5[1]}]', transform=ax.transAxes, fontsize='medium', color='red',zorder=6)
            ax.text(.885, .86, f'[{det5[2]}:{det5[3]}]', transform=ax.transAxes, fontsize='medium', color='blue',zorder=6)
            ax.legend(['209Po', '210Po', 'not counted         '])
            ax.set_ylabel('counts')
            ax.set_xlabel('decay energy channels')
        
        
    elif name == "DET# 6":
        po209 = np.sum(counts[det6[0] : det6[1]])
        po210 = np.sum(counts[det6[2] : det6[3]])
        detID = "EnsembleInput6"
        if PlotSPEs == True:
            fig, ax = plt.subplots(figsize=(10,4))
            ax.fill_between(counts.index.values[det6[0] : det6[1]],0, counts[det6[0] : det6[1]],color='red', zorder=3)
            ax.fill_between(counts.index.values[det6[2] : det6[3]],0, counts[det6[2] : det6[3]],color='blue', zorder=4)
            ax.fill_between(counts.index.values[det6[0]-100 : det6[3]+100],0, counts[det6[0]-100 : det6[3]+100],color='0.7', zorder=1)
            ax.text(.025, .92, 'SPECTRUM INTEGRATION RANGES', transform=ax.transAxes, fontsize='medium', fontweight='bold')
            ax.text(.025, .85, f'Detector: {detID}', transform=ax.transAxes, fontsize='small')
            ax.text(.025, .78, f'File: {sampleID}', transform=ax.transAxes, fontsize='small')
            ax.text(.885, .922, f'[{det6[0]}:{det6[1]}]', transform=ax.transAxes, fontsize='medium', color='red',zorder=6)
            ax.text(.885, .86, f'[{det6[2]}:{det6[3]}]', transform=ax.transAxes, fontsize='medium', color='blue',zorder=6)
            ax.legend(['209Po', '210Po', 'not counted         '])
            ax.set_ylabel('counts')
            ax.set_xlabel('decay energy channels')
        
        
    elif name == "DET# 7":
        po209 = np.sum(counts[det7[0] : det7[1]])
        po210 = np.sum(counts[det7[2] : det7[3]])
        detID = "EnsembleInput7"
        if PlotSPEs == True:
            fig, ax = plt.subplots(figsize=(10,4))
            ax.fill_between(counts.index.values[det7[0] : det7[1]],0, counts[det7[0] : det7[1]],color='red', zorder=3)
            ax.fill_between(counts.index.values[det7[2] : det7[3]],0, counts[det7[2] : det7[3]],color='blue', zorder=4)
            ax.fill_between(counts.index.values[det7[0]-100 : det7[3]+100],0, counts[det7[0]-100 : det7[3]+100],color='0.7', zorder=1)
            ax.text(.025, .92, 'SPECTRUM INTEGRATION RANGES', transform=ax.transAxes, fontsize='medium', fontweight='bold')
            ax.text(.025, .85, f'Detector: {detID}', transform=ax.transAxes, fontsize='small')
            ax.text(.025, .78, f'File: {sampleID}', transform=ax.transAxes, fontsize='small')
            ax.text(.885, .922, f'[{det7[0]}:{det7[1]}]', transform=ax.transAxes, fontsize='medium', color='red',zorder=6)
            ax.text(.885, .86, f'[{det7[2]}:{det7[3]}]', transform=ax.transAxes, fontsize='medium', color='blue',zorder=6)
            ax.legend(['209Po', '210Po', 'not counted         '])
            ax.set_ylabel('counts')
            ax.set_xlabel('decay energy channels')
        
        
    elif name == "DET# 8":
        po209 = np.sum(counts[det8[0] : det8[1]])
        po210 = np.sum(counts[det8[2] : det8[3]])
        detID = "EnsembleInput8"
        if PlotSPEs == True:
            fig, ax = plt.subplots(figsize=(10,4))
            ax.fill_between(counts.index.values[det8[0] : det8[1]],0, counts[det8[0] : det8[1]],color='red', zorder=3)
            ax.fill_between(counts.index.values[det8[2] : det8[3]],0, counts[det8[2] : det8[3]],color='blue', zorder=4)
            ax.fill_between(counts.index.values[det8[0]-100 : det8[3]+100],0, counts[det8[0]-100 : det8[3]+100],color='0.7', zorder=1)
            ax.text(.025, .92, 'SPECTRUM INTEGRATION RANGES', transform=ax.transAxes, fontsize='medium', fontweight='bold')
            ax.text(.025, .85, f'Detector: {detID}', transform=ax.transAxes, fontsize='small')
            ax.text(.025, .78, f'File: {sampleID}', transform=ax.transAxes, fontsize='small')
            ax.text(.885, .922, f'[{det8[0]}:{det8[1]}]', transform=ax.transAxes, fontsize='medium', color='red',zorder=6)
            ax.text(.885, .86, f'[{det8[2]}:{det8[3]}]', transform=ax.transAxes, fontsize='medium', color='blue',zorder=6)
            ax.legend(['209Po', '210Po', 'not counted         '])
            ax.set_ylabel('counts')
            ax.set_xlabel('decay energy channels')
    else:
        print("ERROR det_match_sum: no match found")
    return detID, po209, po210


################################################################################
###                           counts_to_acivity                              ###
#   calculates unsupported 210Pb activity in sediments from alpha decay counts #
#     corrections: decay during sample processing, salt weight, %mud           #
#
#                                    INPUTS:                               
#  INPUT #1: csv produced by the "spe_to_counts" function (see above col info
#  INPUT #2: csv of detector background activity with the following column order
#  INPUT #3: csv of
#
#
#                                   RETURNS:                                                 
#                  A pd.dataframe with the following columns:                                
#
#     detID	..................................................
#     Z_midpt (cm)	..........................................
#     ΔZ (cm)	..............................................
#
#     Δt_in_counting (sec)	..................................
#     Δt_Plate2Count (min)	..................................
#     Δt_Collect2Count (min)	..............................
#     Δt_SpikeCal2Count (min)	..............................
#     Δt_in_counting (min)	..................................
#     t_platingstart	......................................
#     t_countingstart	......................................
#     Counting_StartDate+Time	..............................
#
#     M_pan (g)	..............................................
#     M_WetSed+Pan (g)	......................................
#     M_DrySed+Pan (g)	......................................
#     M_WetChemSed (g) .......................................
#     M_WetChemSed_SaltCorrected (g)	......................
#     siltclay (volfrac)	..................................
#
#     WeightFrac_Water+Salt	..................................
#     WeightFrac_Sed+Salt	..................................
#     VolFrac_water+salt	..................................
#     VolFrac_sed	..........................................
#     Φ_uncorrected (volfrac)	..............................
#     Φ_saltcor (volfrac)	..................................
#     ρ_bulk_wet (g/cm3)	..................................
#     ρ_bulk_dry (g/cm3)	..................................
#
#     209Po_decays (counts)	..................................
#     210Po_decays (counts)	..................................
#     209Po_decays_minus_bkg (counts)	......................
#     210Po_decays_minus_bkg (counts)	......................
#     210Po_DecayCor_Plate2Count	..........................
#     210Pb_DecayCor_Collect2Plate	..........................
#     209Po_DecayCor_SpikeCal2Count	..........................
#
#     BKG counts Pb209	.......................................
#     BKG counts Pb210	.......................................
#     BKG counts (sec)	.......................................
#     209Po_detector_background_activity (cpm)	...............
#     210Po_detector_background_activity (cpm)	...............
#
#     radioisotope_yield (%)	...............................
#     C_i at collection (dpm/g)		...........................
#     C_i at collection, salt correction (dpm/g)	...........
#     C_i excess at collection, salt correction (dpm/g)	.......
#     C_i excess at collection, salt+mud correction (dpm/g)	...
#
#     Po209_counting_error	...................................
#     Po210_counting_error	...................................
#     pipette_error	...........................................
#     spike_error	...........................................
#     Error_total	...........................................
#     Error	Error_SaltCorr	...................................
#     Error_MudSaltCorr (Xdir_error)	.......................
###                                                                          ###
################################################################################


def counts_to_activity(counts_fname, bkg_fname, supLvl):
    # Imports
    import matplotlib.pyplot as plt
    from datetime import datetime
    import numpy as np
    import pandas as pd

    ################################################################################
    ###                             DEFINE CONSTANTS.                            ###
    
    # core collection date
    t_collection_yCE = datetime.strptime("10/15/2021", "%m/%d/%Y")  
    # spike calibration date
    t_spikeCal = datetime.strptime("08/15/2016", "%m/%d/%Y")  
    # volume of spike used per sample, ml
    spike_volume_ml = 0.998     
    # spike activity at the time of calibration, dpm/ml
    C_spike_atCal_dpmml = 12.0469862348134  
    # the uncertainty associated with spike activity
    u_C_spike_atCal_dpmml = 0.4  
    
    # decay constant of 210Pb, in min^-1
    λ_210Pb_min = 0.00000005914  
    # decay constant of 210Po, in min^-1
    λ_210Po_min = 0.000003472848 
    # decay constant of 209Po, in min^-1
    λ_209Po_min = 0.00000001292  
    
    # density of porewater, g/cm^3
    ρ_porewater_gcm3 = 1.025     
    # density of sediment, g/cm^3
    ρ_particle_gcm3 = 2.65       
    # mass salt fraction of seawater
    porewater_saltFrac = 0.025   
    
    ###                              END OF CONSTANTS                            ###
    ################################################################################

    print(
        "|----------------------  COUNTS2ACTIVITY STARTED  ----------------------|")
    print(
        f"||    Data from file:       {counts_fname}                               "
    )
    print(
        f"||    Detector backgrounds: {bkg_fname}                                  "
    )
    print(
        f"||    Supported level:      {supLvl} dpm/g                               "
    )

    # read in csv files specified by user
    cts = pd.read_csv(counts_fname, float_precision="round_trip", header=0)
    bkg = pd.read_csv(bkg_fname)

    # CALCULATE SEDIMENT BULK DENSITY AND POROSITY
    # weight fraction of water and salt in preprocessed aliquot
    cts["WeightFrac_Water+Salt"] = (
        (cts["M_WetSed+Pan (g)"] - cts["M_pan (g)"])
        - (cts["M_DrySed+Pan (g)"] - cts["M_pan (g)"])
    ) / (cts["M_WetSed+Pan (g)"] - cts["M_pan (g)"])
    # weight fraction of sediment and salt in preprocessed aliquot
    cts["WeightFrac_Sed+Salt"] = 1 - cts["WeightFrac_Water+Salt"]
    # weight of sediment aliquot used in wet chemistry
    cts["M_WetChemSed_SaltCorrected (g)"] = cts["M_WetChemSed (g)"] - (
        (cts["WeightFrac_Water+Salt"] * porewater_saltFrac)
        * (cts["M_WetChemSed (g)"] / cts["WeightFrac_Sed+Salt"])
    )
    # volume fraction of water and salt using density assumptions
    cts["VolFrac_water+salt"] = (
        cts["WeightFrac_Water+Salt"] / (1 - porewater_saltFrac)
    ) * (1 / ρ_porewater_gcm3)
    # volume fraction of sediment only using density assumptions
    cts["VolFrac_sed"] = (
        cts["WeightFrac_Sed+Salt"]
        - (
            cts["WeightFrac_Water+Salt"]
            * (porewater_saltFrac / (1 - porewater_saltFrac))
        )
    ) / (ρ_particle_gcm3)
    # porosity of sediment with no salt correction
    cts["Φ_uncorrected (volfrac)"] = (
        cts["WeightFrac_Water+Salt"] * ρ_particle_gcm3
    ) / (
        (cts["WeightFrac_Water+Salt"] * ρ_particle_gcm3)
        + (1 - cts["WeightFrac_Water+Salt"]) * ρ_porewater_gcm3
    )
    # porosity after salt correction
    cts["Φ_saltcor (volfrac)"] = (cts["VolFrac_water+salt"]) / (
        cts["VolFrac_water+salt"] + cts["VolFrac_sed"]
    )
    # material bulk density before drying
    cts["ρ_bulk_wet (g/cm3)"] = (1 - cts["Φ_saltcor (volfrac)"]) * (
        ρ_particle_gcm3
    ) + (cts["Φ_saltcor (volfrac)"]) * (ρ_porewater_gcm3)
    # material bulk density after drying
    cts["ρ_bulk_dry (g/cm3)"] = (1 - cts["Φ_saltcor (volfrac)"]) * (
        ρ_particle_gcm3
    )

    
    
    
    
    # TIME CONVERSIONS & CALCULATIONS
    # elapsed minutes plated planchets spent in counting
    cts["Δt_in_counting (min)"] = cts["Δt_in_counting (sec)"] / 60
    # elapsed minutes spent counting with no sample to measure bkg decays
    bkg["Δt_in_counting (min)"] = bkg["counting time (sec)"] / 60
    # time of plating in datetime format
    cts["t_platingstart"] = (
        cts["Plating_StartDate (DD/MM/YYYY)"]
        + " "
        + cts["Plating_StartTime (HH:MM:SS)"]
    )
    cts["t_platingstart"] = pd.to_datetime(
        cts["t_platingstart"], infer_datetime_format=True
    )
    # time of counting in datetime format
    cts["t_countingstart"] = (
        cts["Counting_StartDate"] + " " + cts["Counting_StartTime"]
    )
    cts["t_countingstart"] = pd.to_datetime(
        cts["t_countingstart"], infer_datetime_format=True
    )

    # read in a csv of background activity with the following column names
    # "Detector Name", "counts Po209", "counts Po210", "counting time (sec)"
    # 209Po background activity of the α-counting cubby
    bkg["209Po_detector_background_activity (cpm)"] = (
        bkg["counts Po209"] / bkg["Δt_in_counting (min)"]
    )
    # 210Po background activity of the α-counting cubby
    bkg["210Po_detector_background_activity (cpm)"] = (
        bkg["counts Po210"] / bkg["Δt_in_counting (min)"]
    )

   






    # CORRECT FOR THE BACKGROUND ACTIVITY OF EACH DETECTOR
    a = []
    for i in range(len(cts)):
        a.append(bkg.loc[bkg["Detector Name"] == cts["detID"][i]].values[0])
    a = pd.DataFrame(a)
    a = a.drop(labels=[0], axis=1)
    a = a.rename(
        columns={
            1: "BKG counts Pb209",
            2: "BKG counts Pb210",
            3: "BKG counts (sec)",
            5: "209Po_detector_background_activity (cpm)",
            6: "210Po_detector_background_activity (cpm)",
        }
    )
    cts = pd.concat([a, cts], axis=1)
    # total 209Po α-counts from planchet only (background decays removed)
    cts["209Po_decays_minus_bkg (counts)"] = cts["209Po_decays (counts)"] - (
        cts["Δt_in_counting (min)"]
        * cts["209Po_detector_background_activity (cpm)"]
    )
    # total 210Po α-counts from planchet only (background decays removed)
    cts["210Po_decays_minus_bkg (counts)"] = cts["210Po_decays (counts)"] - (
        cts["Δt_in_counting (min)"]
        * cts["210Po_detector_background_activity (cpm)"]
    )
    print(f"||   ")
    print(f"||    ...calculating sample activites... ")

    # CALCULATE CORRECTION VALUES FOR DECAY OF ISOTOPES
    # elapsed time between plating and counting
    cts["Δt_Plate2Count (min)"] = (
        cts["t_countingstart"] - cts["t_platingstart"]
    ) / np.timedelta64(1, "m")
    # elapsed time between plating and counting
    cts["Δt_Collect2Count (min)"] = (
        cts["t_platingstart"] - t_collection_yCE
    ) / np.timedelta64(1, "m")
    # elapsed time between spike calibration and counting
    cts["Δt_SpikeCal2Count (min)"] = (
        cts["t_countingstart"] - t_spikeCal
    ) / np.timedelta64(1, "m")
    # fraction of 210Po remaining after the time elapsed between plating and α-counting
    cts["210Po_DecayCor_Plate2Count"] = np.exp(
        -λ_210Po_min * cts["Δt_Plate2Count (min)"]
    )
    # fraction of 210Pb remaining after the time elapsed between collection and plating
    cts["210Pb_DecayCor_Collect2Plate"] = np.exp(
        -λ_210Pb_min * cts["Δt_Collect2Count (min)"]
    )
    # fraction of 209Po remaining after the time elapsed between spike calibration and α-counting
    cts["209Po_DecayCor_SpikeCal2Count"] = np.exp(
        -λ_209Po_min * cts["Δt_SpikeCal2Count (min)"]
    )

    # 210Pb concentration of sediment section i at the date of collection
    cts["C_i at collection (dpm/g)"] = (
        (
            cts["210Po_decays_minus_bkg (counts)"]
            / (cts["M_WetChemSed (g)"] * cts["210Po_DecayCor_Plate2Count"])
        )
        * (  # 210Po activity of mud
            (
                spike_volume_ml
                * C_spike_atCal_dpmml
                * cts["209Po_DecayCor_SpikeCal2Count"]
            )
            / cts["209Po_decays_minus_bkg (counts)"]
        )
        * (  # scaled by 209Po yield
            1 / cts["210Pb_DecayCor_Collect2Plate"]
        )  # scaled by 210Pb decay
    )
    # 210Pb concentration of sediment section i at the date of collection after correcting for the mass of salt
    cts["C_i at collection, salt correction (dpm/g)"] = (
        cts["C_i at collection (dpm/g)"]
        * cts["M_WetChemSed (g)"]
        / cts["M_WetChemSed_SaltCorrected (g)"]
    )
    # 210Pb excess concentration of sediment section i at the date of collection after correcting for the mass of salt
    cts["C_i excess at collection, salt correction (dpm/g)"] = (
        cts["C_i at collection, salt correction (dpm/g)"] - supLvl
    )
    # 210Pb excess concentration of sediment section i at the date of collection after correcting for the mass of salt and the fraction of silt present
    cts["C_i excess at collection, salt+mud correction (dpm/g)"] = (
        cts["C_i excess at collection, salt correction (dpm/g)"]
        / cts["siltclay (volfrac)"]
    )

    # the percentage of the spike successfully measured in the alpha counter
    cts["radioisotope_yield (%)"] = (
        cts["209Po_decays_minus_bkg (counts)"]
        / (
            C_spike_atCal_dpmml
            * spike_volume_ml
            * cts["Δt_in_counting (sec)"]
            / 60
        )
        * 100
    )

    # calculate error
    print(
        f"||    ...calculating error...                                            "
    )
    # ☑
    cts["Po209_counting_error"] = (
        (cts["209Po_decays (counts)"]) ** (1 / 2)
    ) / (cts["209Po_decays (counts)"])
    # ☑
    cts["Po210_counting_error"] = (
        (cts["210Po_decays (counts)"]) ** (1 / 2)
    ) / (cts["210Po_decays (counts)"])
    # ☑
    cts["pipette_error"] = 0.003 / spike_volume_ml
    # ☑
    cts["spike_error"] = u_C_spike_atCal_dpmml / C_spike_atCal_dpmml
    # ☑
    cts["Error_total"] = (
        ((cts["pipette_error"]) ** 2)
        + ((cts["spike_error"]) ** 2)
        + ((cts["Po210_counting_error"]) ** 2)
        + ((cts["Po209_counting_error"]) ** 2)
    ) ** (1 / 2)
    # ☑
    cts["Error"] = cts["C_i at collection (dpm/g)"] * cts["Error_total"]
    # ☑
    cts["Error_SaltCorr"] = np.abs(
        (cts["Error"] / cts["C_i at collection (dpm/g)"])
        * (cts["C_i at collection, salt correction (dpm/g)"])
    )
    # ☑
    cts["Error_MudSaltCorr (Xdir_error)"] = (
        cts["Error_SaltCorr"] / cts["siltclay (volfrac)"]
    )

    print(
        "||    ...cleaning df...                                                  "
    )
    # DROP REDUNDANT COLS
    cts = cts.drop(
        labels=[
            "Counting_StartDate",
            "Counting_StartTime",
            "Plating_StartDate (DD/MM/YYYY)",
            "Plating_StartTime (HH:MM:SS)",
        ],
        axis=1,
    )

    print(
        '||    ...created df named "cts".                                          '
    )
    print(
        "|---------------------  COUNTS2ACTIVITY FINISHED  -----------------------|"
    )
    print("   ")
    return cts
