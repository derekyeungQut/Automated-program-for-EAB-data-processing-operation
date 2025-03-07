import tkinter as tk
import customtkinter as ctk
from tkinter import filedialog, OptionMenu
import pandas as pd
import os
from pathlib import Path
from scipy.signal import savgol_filter
import numpy as np
import matplotlib.pyplot as plt
import sys
import re
ctk.set_appearance_mode("dark")
ctk.set_default_color_theme("dark-blue")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ initializing global variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

minON = 0 # minimum ON current value (from entry box - Calculation Page)
minOFF = 0 # minimum OFF current value (from entry box - Calculation Page)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ initializing global variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#------------------------------------  Helper functions  ---------------------------------------------#

# Helper Function (called by openfolder function): 
# 1) Find all subfolders in the given directory
# 2) Find all .txt files in each subfolder
# 3) Return the list of .txt files in the ON and OFF subfolders
def find_txt_files(folderpath):
    subfolders = [f.path for f in os.scandir(folderpath) if f.is_dir()]
    subfolder_names = [os.path.basename(f.path) for f in os.scandir(folderpath) if f.is_dir()]

    on_files = []
    on_files_names = []
    off_files = []
    off_files_names = []

    if not subfolders:  # Check if there are no subfolders
        print("No subfolders found in the given directory.")
    else:
        print(f"Subfolders found: {subfolder_names}")  # Debugging line to check subfolder list

        for subfolder in subfolders:
            input_dir = Path(subfolder)
            # Use Path.glob to find all .txt files in the subfolder
            txt_files = [str(file) for file in input_dir.glob("*.txt")]

            subfolder_name = os.path.basename(subfolder)

            if "ON" in subfolder_name:
                on_files.extend(txt_files)
                on_files_names.extend([os.path.basename(file) for file in txt_files])
                print(f"Found {len(txt_files)} text files in: {subfolder_name}")
            elif "OFF" in subfolder_name:
                off_files.extend(txt_files)
                off_files_names.extend([os.path.basename(file) for file in txt_files])
                print(f"Found {len(txt_files)} text files in: {subfolder_name}")

    return on_files, off_files, on_files_names, off_files_names

# Helper Function  (called by openfolder function):
# 1) Sort the list of files by concentration
def sort_files_by_concentration(file_list):
    return sorted(file_list, key=get_concentration)

# Helper function (called by generate_dataframes function):
# 1) Extract the time from the file name
def get_time(file_path):
    # Extract the file name from the file path
    file_name = file_path.split('\\')[-1]
    
    # Use regex to find the first four digits in the file name
    match = re.search(r'(\d{2})(\d{2})', file_name)
    if match:
        hours = match.group(1)
        minutes = match.group(2)
        return f"{hours}:{minutes}"
    else:
        return None

# Helper function (called by generate_dataframes function):
# 1) Extract the concentration from the file name
def get_concentration(file_name):
    # First, check if PBS is present
    if "PBS" in file_name:
        return 0  # PBS is treated as 0 molar
    
    # Extract numerical concentration using regex
    match = re.search(r'(\d+\.?\d*)uM', file_name)  # Matches values like 2.5uM, 100uM, etc.
    if match:
        concentration = float(match.group(1))  # Convert to float for numerical comparison
        return concentration
    else:
        return float('inf')  # If no concentration is found, return a high value

# Helper function (called by generate_calibration_dataframes function):
# 1) Calculate the average of every n elements in an array
def average(array,n):

    avg = []
    averaged_array = []
    count = 0  

    for i in array:
        avg.append(i)
        count += 1
        if count == n:
            average_value = sum(avg) / len(avg)
            averaged_array.append(average_value)
            avg = []
            count = 0
    # If there are remaining elements in avg, average them as well
    if avg:
        average_value = sum(avg) / len(avg)
        averaged_array.append(average_value)

    return averaged_array

# Helper function (called by plot_data function):
# 1) Parse the time string into a float
def parse_time(t):
    hours, minutes = map(int, t.split(':'))
    return hours + (minutes / 60)

#---------------------------------- End of Helper functions -------------------------------------------#


#=======================================  Main functions  ==============================================#

# BUTTON Function (from Main Page): 
# 1) Open the folder and store the list of files in sorted_on_files and sorted_off_files
# 2) Call the calculation page
def openfolder():
    global folder_path
    folder_path = filedialog.askdirectory()
    print(f"Selected folder: {folder_path}")
    
    on_txt_files, off_txt_files, on_filenames, off_filenames = find_txt_files(folder_path)
    
    global sorted_on_files, sorted_off_files
    sorted_on_files = sort_files_by_concentration(on_txt_files)
    sorted_off_files = sort_files_by_concentration(off_txt_files)

    calculationPage()
    
# CALCULATION Function (called by sub-functions of generate_dataframes):
# 1) Extract peak currents and potentials for a single given file path
def MultichannelPeakheight(file_path):

    ################## FILE DIRECTORY AND DATAFRAME SET UP ############################
    # Load the file and skip initial irrelevant lines until the header
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Find the line with the header
    start_line = next(i for i, line in enumerate(lines) if "Potential/V" in line)

    # Read the relevant lines into a DataFrame
    rawdf = pd.read_csv(file_path, skiprows=start_line, sep=',')

    # Clean column names by stripping leading and trailing spaces
    rawdf.columns = rawdf.columns.str.strip()
    ################## FILE DIRECTORY AND DATAFRAME SET UP ############################
    
    channels_i_A = []

    for i in rawdf.columns:
        if 'd/A' in i:
            filtered_i = savgol_filter(rawdf[i],window_length = 10,polyorder = 1)
            channels_i_A.append(filtered_i)
        elif 'Diff' in i:
            filtered_i = savgol_filter(rawdf[i],window_length = 10,polyorder = 1)
            channels_i_A.append(filtered_i)

    data = {'Potential/V': rawdf['Potential/V']}

    for i in range(len(channels_i_A)):
        data[f'i{i+1}d_filtered/A'] = channels_i_A[i]

    df = pd.DataFrame(data)

    v = df['Potential/V']
    # Split the data into 8 parts
    fit_half = round(len(v) / 2)
    fit_quarter = round(len(v) / 4)
    fit_eighth = round(len(v) / 8)

    pk = []
    pk_potential= []


    for i in df.columns[1:]:
        I_diff_filtered = df[i]
        mean_1st = np.mean(I_diff_filtered[fit_eighth:-fit_quarter- fit_half]) # mean from 2nd eigth
        mean_last = np.mean(I_diff_filtered[fit_quarter + fit_half:-fit_eighth]) # mean from 2nd last eigth

        # Find the index of the closest value to mean_last and mean_1st in eval_regression
        index_mean1st = (np.abs(np.array(I_diff_filtered[fit_eighth:-fit_quarter-fit_half]) - mean_1st)).argmin()  + fit_eighth
        index_meanlast = (np.abs(np.array(I_diff_filtered[fit_quarter + fit_half:-fit_eighth]) - mean_last)).argmin() + fit_half + fit_quarter 

        # Get the corresponding potential values
        potential_meanfirstC = v[index_mean1st]
        potential_meanlastC = v[index_meanlast]
        ################## DATA SMOOTHING AND MANIPULATION ################################
        # Find the index of the maximum value in the 'Diff(i/A)' column
        max_index = abs(df[i][:fit_half+fit_quarter+fit_eighth]).idxmax()


        baseline_v = list(np.arange(potential_meanfirstC , potential_meanlastC -0.001, -0.001))
        baseline_v = [round(num, 3) for num in baseline_v]

        

        baseline_i = np.linspace(mean_1st, mean_last, len(baseline_v)).tolist()


        if (max_index - index_mean1st) < 0 :
            peakheight = np.nan
            pk.append(peakheight)
            pk_potential.append(peakheight)
        elif (max_index - index_mean1st) >= len(baseline_i):
            peakheight = np.nan
            pk.append(peakheight)
            pk_potential.append(peakheight)
        else:
            peakheight = abs(baseline_i[max_index - index_mean1st] - I_diff_filtered[max_index])
            pk.append(peakheight)
            pk_potential.append(abs(v[max_index]))

    return pk, pk_potential

# BUTTON Function (From Calculation Page):
# 1) recall selected mode from main page
# 2) determine which dataframes to generate
def generate_dataframes():
    mode = selected_modes.get()
    print(f"Selected mode: {mode}")
    if mode == 'Degradation Curve':
        generate_timeframe_dataframes()
    elif mode == 'Calibration Curve':
        generate_calibration_dataframes()

# Sub-Function :
# 1) Generate dataframes for the 'Extract Peaks Over Time' mode
def generate_timeframe_dataframes():
    num_channels = selected_channels.get()
    on_data = {'Time': []}
    for i in range(num_channels):
        on_data[f'ipeak{i+1}'] = []
        on_data[f'vpeak{i+1}'] = []


    off_data = {'Time': []}
    for i in range(num_channels):
        off_data[f'ipeak{i+1}'] = []
        off_data[f'vpeak{i+1}'] = []

    for file in sorted_on_files:
        pk, pk_potential = MultichannelPeakheight(file)
        if len(pk) != num_channels:
            print(f"Error: {file} has {len(pk)} channels instead of {num_channels}")
        else:
            time = get_time(file)
            on_data['Time'].append(time)
            for i in range(len(pk)):
                on_data[f'ipeak{i+1}'].append(pk[i])
                on_data[f'vpeak{i+1}'].append(pk_potential[i])

    global on_df, off_df

    on_data['Time'] = [parse_time(t) for t in on_data['Time']]

    on_df = pd.DataFrame(on_data)
    print('On Data: \n',on_df)

    for file in sorted_off_files:
        pk, pk_potential = MultichannelPeakheight(file)
        if len(pk) != num_channels:
            print(f"Error: {file} has {len(pk)} channels instead of {num_channels}")
        else:
            time = get_time(file)
            off_data['Time'].append(time)
            for i in range(len(pk)):
                off_data[f'ipeak{i+1}'].append(pk[i])
                off_data[f'vpeak{i+1}'].append(pk_potential[i])

    off_data['Time'] = [parse_time(t) for t in off_data['Time']]
    off_df = pd.DataFrame(off_data)
    print('Off Data: \n',off_df)

# Sub-Function :
# 1) Generate dataframes for the 'Calibrate Over Concentration' mode
def generate_calibration_dataframes():
    num_channels = selected_channels.get()
    on_data = {'Concentration': []}
    off_data = {'Concentration': []}

    for i in range(num_channels):
        on_data[f'ipeak{i+1}'] = []
        on_data[f'vpeak{i+1}'] = []
        off_data[f'ipeak{i+1}'] = []
        off_data[f'vpeak{i+1}'] = []

    for file in sorted_on_files:
        pk, pk_potential = MultichannelPeakheight(file)
        if len(pk) != num_channels:
            print(f"Error: {file} has {len(pk)} channels instead of {num_channels}")
        else:
            conc = get_concentration(file)
            on_data['Concentration'].append(conc)
            for i in range(len(pk)):
                on_data[f'ipeak{i+1}'].append(pk[i])
                on_data[f'vpeak{i+1}'].append(pk_potential[i])

    on_df = pd.DataFrame(on_data)

    for file in sorted_off_files:
        pk, pk_potential = MultichannelPeakheight(file)
        if len(pk) != num_channels:
            print(f"Error: {file} has {len(pk)} channels instead of {num_channels}")
        else:
            conc = get_concentration(file)
            off_data['Concentration'].append(conc)
            for i in range(len(pk)):
                off_data[f'ipeak{i+1}'].append(pk[i])
                off_data[f'vpeak{i+1}'].append(pk_potential[i])

    off_df = pd.DataFrame(off_data)
    
    global peakheights_df, averaged_df, calibration_df
    peakheights_df = pd.merge(on_df, off_df, left_index=True, right_index=True, suffixes=('_ON', '_OFF'))
    

    if minON > 0 or minOFF > 0:
        print('Failsafe applied, values below the minimum will be set to NaN')
        for column in peakheights_df.columns:
            if 'ON' in column:
                peakheights_df[column] = peakheights_df[column].apply(lambda x: np.nan if x < minON else x)
            elif 'OFF' in column:
                peakheights_df[column] = peakheights_df[column].apply(lambda x: np.nan if x < minOFF else x)
    else:
        print('No failsafe applied')
    
    print('All Peak Currents and Potentials Data: \n',peakheights_df)
    
    averaged_df = pd.DataFrame()
    for i in peakheights_df.columns:
        averaged_df[i] = average(peakheights_df[i],3)
    print('Data averaged over', selected_trials.get() ,' trials: \n',averaged_df)


    calibration_data = {'Concentration': averaged_df['Concentration_ON']}

    for i in range(num_channels):
        calibration_data[f'signalONchange{i+1}'] = []
        calibration_data[f'signalOFFchange{i+1}'] = []
        calibration_data[f'KDM{i+1}'] = []

    def calibrate(electrode ,df):
        iON = df[f'ipeak{electrode+1}_ON']
        iOFF = df[f'ipeak{electrode+1}_OFF']
        kdm = 100* ((iON / iON[0]) - (iOFF / iOFF[0]))
        onschange = 100* (iON - iON[0]) / iON[0]
        offchange = 100* (iOFF - iOFF[0]) / iOFF[0]
        calibration_data[f'KDM{i+1}'] = kdm
        calibration_data[f'signalONchange{i+1}'] = onschange
        calibration_data[f'signalOFFchange{i+1}'] = offchange
        

    for i in range(num_channels):
        calibrate(i, averaged_df)

    calibration_df = pd.DataFrame(calibration_data)
    print('Calibration Data: \n',calibration_df)

# BUTTON Function (From Calculation Page):
# 1) Recall the selected mode and number of channels
# 2) Plot the data 
def plotdata():
    # Disable the topmost attribute before showing the plot
    calculationWindow.attributes('-topmost', False)

    mode = selected_modes.get()
    num_channels = selected_channels.get()
    if mode == 'Degradation Curve':
        # Convert time strings to timedelta objects
        # on_time_values = [parse_time(t) for t in on_df['Time']]
        # off_time_values = [parse_time(t) for t in off_df['Time']]
        fig1, axs1 = plt.subplots(2, 1, figsize=(15, 15))
        for i in range(num_channels):
            axs1[0].plot(on_df['Time'], on_df[f'ipeak{i+1}'], label=f'Electrode {i+1} ON')
            axs1[0].plot(off_df['Time'], off_df[f'ipeak{i+1}'], label=f'Electrode {i+1} OFF')
            axs1[0].set_title('Peak current')
            axs1[0].set_xlabel('Time')
            axs1[0].set_ylabel('Peak currents')
            axs1[0].legend()

            
            axs1[1].plot(on_df['Time'], on_df[f'vpeak{i+1}'], label=f'Electrode {i+1} ON')
            axs1[1].plot(off_df['Time'], off_df[f'vpeak{i+1}'], label=f'Electrode {i+1} OFF')
            axs1[1].set_title('Peak potential')
            axs1[1].set_xlabel('Time')
            axs1[1].set_ylabel('Peak potentials')
            axs1[1].legend()

        plt.show()

    elif mode == 'Calibration Curve':
        fig, ax = plt.subplots(1, num_channels, figsize=(6*num_channels, 5))
        if num_channels == 1:
            ax.plot(calibration_df['Concentration'], calibration_df['signalONchange1'], label='Channel 1 ON')
            ax.plot(calibration_df['Concentration'], calibration_df['signalOFFchange1'], label='Channel 1 OFF')
            ax.plot(calibration_df['Concentration'], calibration_df['KDM1'], label='Channel 1 KDM')
            ax.set_xlabel('Concentration (uM)')
            ax.set_xscale('log')
            ax.set_ylabel('Signal Change (%)')
            ax.set_title('Channel 1 Signal Change vs Concentration')
            ax.legend()
        else:
            for i in range(num_channels):
                ax[i].plot(calibration_df['Concentration'], calibration_df[f'signalONchange{i+1}'], label=f'Channel {i+1} ON')
                ax[i].plot(calibration_df['Concentration'], calibration_df[f'signalOFFchange{i+1}'], label=f'Channel {i+1} OFF')
                ax[i].plot(calibration_df['Concentration'], calibration_df[f'KDM{i+1}'], label=f'Channel {i+1} KDM')
                ax[i].set_xlabel('Concentration (uM)')
                ax[i].set_xscale('log')
                ax[i].set_ylabel('Signal Change (%)')
                ax[i].set_title(f'Channel {i+1} Signal Change vs Concentration')
                ax[i].legend()
        plt.show()

# BUTTON Function (From Calculation Page):
# 1) Update the failsafe values
# Will set values below the minimum to NaN 
# (won't affect original dataframes if button is not pressed)
def failsafeupdate():
    global minON, minOFF
    if minimumON.get() == '':
        minON = 0
    else:
        try:
            minON = float(minimumON.get())
        except:
            print("Error: Minimum ON value must be a number.")
            return
    if minimumOFF.get() == '':
        minOFF = 0
    else:
        try:
            minOFF = float(minimumOFF.get())
        except:
            print("Error: Minimum OFF value must be a number.")
            return
        
    print(f"Failsafe updated: Minimum ON = {minON}, Minimum OFF = {minOFF}")

# BUTTON Function (From Calculation Page):
# 1) save the dataframes as an Excel file
def savefile():
    # Disable the topmost attribute before showing the plot
    calculationWindow.attributes('-topmost', False)
    
    #########################################################################################
    # Let the user pick where to save the Excel file
    root = tk.Tk()
    root.withdraw()  # Hide the root window
    file_path = filedialog.asksaveasfilename(defaultextension=".xlsx", filetypes=[("Excel files", "*.xlsx"), ("All files", "*.*")])

    mode = selected_modes.get()

    if file_path:
        with pd.ExcelWriter(file_path) as writer:
            if mode == 'Extract Peaks Over Time':
                on_df.to_excel(writer, sheet_name='ON Peaks', index=False)
                off_df.to_excel(writer, sheet_name='OFF Peaks', index=False)
            elif mode == 'Calibrate Over Concentration':
                calibration_df.to_excel(writer, sheet_name='Calibration Table', index=False)
                averaged_df.to_excel(writer, sheet_name='Averaged Table', index=False)
                peakheights_df.to_excel(writer, sheet_name='All peak data', index=False)

        print(f"DataFrame exported to {file_path}")
    else:
        print("Save operation cancelled.")
    #########################################################################################

# Button Function (From Calculation Page):    
# 1) close the program
def closeprogram():   
    # Close the main window to terminate the program
    InitialWindow.destroy()
    sys.exit()

#=========================================== End of Main Functions ================================================#



############################################## GUI Windows #########################################################

#__________________________________ Calculation window ________________________________#
# * called by openfolder function    
def calculationPage():

    global calculationWindow 
    calculationWindow = ctk.CTkToplevel(InitialWindow)
    calculationWindow.geometry("600x450")
    calculationWindow.title("calculation window")
    # Set the new window to stay on top
    calculationWindow.lift()
    calculationWindow.attributes('-topmost', True)

    generatedataframe = ctk.CTkButton(calculationWindow, text = "Generate Data Frames", command = generate_dataframes)
    generatedataframe.pack(padx =20, pady = 20)

    plotbutton = ctk.CTkButton(calculationWindow, text = "Plot Data", command = plotdata)
    plotbutton.pack(padx =20, pady = 20)

    frame = ctk.CTkFrame(master=calculationWindow)
    frame.pack(padx=20, pady=10, fill = 'both' , expand = True)

    minimumONlabel = ctk.CTkLabel(master=frame, text = "Minimum ON current(A)", font = ('Helvetica', 12))
    minimumONlabel.grid(row=0, column=0, padx = 40)
    global minimumON
    minimumON = tk.StringVar()
    minimumONentry = ctk.CTkEntry(master=frame, textvariable=minimumON)
    minimumONentry.grid(row=1, column=0, padx = 40)

    minimumOFFlabel = ctk.CTkLabel(master=frame, text = "Minimum OFF current(A)", font = ('Helvetica', 12))
    minimumOFFlabel.grid(row=0, column=2, padx = 40)
    global minimumOFF
    minimumOFF = tk.StringVar()
    minimumOFFentry = ctk.CTkEntry(master=frame, textvariable=minimumOFF)
    minimumOFFentry.grid(row=1, column=2, padx = 40)

    failsafebutton = ctk.CTkButton(master=frame, text = "Apply Fail-safe", command = failsafeupdate)
    failsafebutton.grid(row=2, column=1)

    savebutton = ctk.CTkButton(calculationWindow, text = "Save as Excel file", command = savefile)
    savebutton.pack(padx =20)

    closebutton = ctk.CTkButton(calculationWindow, text = "exit program", command = closeprogram)
    closebutton.pack(padx =20, pady = 10)

    calculationWindow.mainloop()
#__________________________________ Calculation window ________________________________#

#____________________________ Main Window ______________________________________#
InitialWindow =  ctk.CTk()
InitialWindow.geometry("750x300")
InitialWindow.title("Automated Data Extraction Program")
Title = ctk.CTkLabel(InitialWindow, text = "Data Extraction Program", font = ('Helvetica', 18))
Title.grid(row=0, column=1, padx=60, pady=20)

# Drop down menu to select the mode of operation
mode_label = ctk.CTkLabel(InitialWindow, text = "Select the mode of operation", font = ('Helvetica', 12))
mode_label.grid(row=2, column=0, padx=40)
selected_modes = tk.StringVar() # selected mode saved as string (from drop down menu - Main Window)
modedrop = ctk.CTkOptionMenu(InitialWindow, variable= selected_modes, values=['Degradation Curve', 'Calibration Curve'])
modedrop.grid(row=3, column=0, padx=0)

# Drop down menu to select the number of channels 
channellabel = ctk.CTkLabel(InitialWindow, text = "Select the number of channels", font = ('Helvetica', 12))
channellabel.grid(row=2, column=1)
selected_channels = tk.IntVar() # selected no. of channels saved as integer (from drop down menu - Main Window)
channeldrop = ctk.CTkOptionMenu(InitialWindow, variable= selected_channels, values=['1', '2', '3', '4', '5', '6'])
channeldrop.grid(row=3, column=1)

# Drop down menu to select the number of trials
trialslabel = ctk.CTkLabel(InitialWindow, text = "Select the number of trials", font = ('Helvetica', 12))
trialslabel.grid(row=2, column=2)
selected_trials = tk.IntVar() # selected no. of trials saved as integer (from drop down menu - Main Window)
channeldrop = ctk.CTkOptionMenu(InitialWindow, variable= selected_trials, values=['1', '2', '3', '4', '5', '6'])
channeldrop.grid(row=3, column=2)

#  Button to open the folder, store list of files and call the calculation page
openfolderbutton = ctk.CTkButton(InitialWindow, text = "Select Main Folder", command = openfolder)
openfolderbutton.grid(row=4, column=1, pady=20)

InitialWindow.mainloop()
#____________________________ Main Window ______________________________________#

#3################################ End of GUI Windows ######################################################## 