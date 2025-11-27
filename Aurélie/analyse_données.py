import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime as dt
from warnings import filterwarnings


filterwarnings("ignore")    # pour enlever des warnings de np.nanmean() gossants et inutiles

#lire les données
data = pd.read_csv("Dataset of weighing station temperature measurements.csv",delimiter=";")
time_str = data["Time"]
temp = data["Outdoor temperature [deg. C]"]
humid = data["Outdoor relative humidity [%]"]
low = data.loc[:,"T[degC]-Low-S1":"T[degC]-Low-S29"].to_numpy()
mid = data.loc[:,"T[degC]-Mid-S1":"T[degC]-Mid-S29"].to_numpy()
top = data.loc[:,"T[degC]-Top-S1":"T[degC]-Top-S29"].to_numpy()

format = "%Y-%m-%d %H:%M"
time = np.array([dt.strptime(t, format) for t in time_str])

cut_time = dt(2024,1,26,hour=0,minute=0)



def check_stratification(plot=True, time_range="all", sensor_range=None):
    '''time_range:
            "all" (or anything): take whole range
            "before": data before 2024-01-26
            "after": data after 2024-01-26
        sensor_range: tuple of 1st and last sensor to check'''
    
    # si on veut limiter les données à des périodes ou senseurs précis
    data1, data2, data3 = low, mid, top
    xdata = time
    cut = np.where(time == cut_time)[0][0]
    if time_range == "before":
        data1, data2, data3 = data1[:cut,:], data2[:cut,:], data3[:cut,:]
        xdata = xdata[:cut]
    elif time_range == "after":
        data1, data2, data3 = data1[cut:,:], data2[cut:,:], data3[cut:,:]
        xdata = xdata[cut:]
    if sensor_range is not None:
        data1, data2, data3 = data1[:,sensor_range[0]-1:sensor_range[1]], data2[:,sensor_range[0]-1:sensor_range[1]], data3[:,sensor_range[0]-1:sensor_range[1]]
    
    # vérification des écarts
    diff1 = data2 - data1
    diff2 = data3 - data2
    diff3 = data3 - data1
    print(f"mean = {round(np.nanmean(diff1),2)}, std = {round(np.nanstd(diff1),2)}, median = {round(np.nanmedian(diff1),2)}, "
          f"max = {round(np.nanmax(diff1),2)}, min = {round(np.nanmin(diff1),2)}")
    print(f"mean = {round(np.nanmean(diff2),2)}, std = {round(np.nanstd(diff2),2)}, median = {round(np.nanmedian(diff2),2)}, "
          f"max = {round(np.nanmax(diff2),2)}, min = {round(np.nanmin(diff2),2)}")
    print(f"mean = {round(np.nanmean(diff3),2)}, std = {round(np.nanstd(diff3),2)}, median = {round(np.nanmedian(diff3),2)}, "
          f"max = {round(np.nanmax(diff3),2)}, min = {round(np.nanmin(diff3),2)}")

    if plot:
        plt.rcParams['date.converter'] = 'concise'
        plt.plot(xdata, np.nanmean(data1, axis=1), label='Low')
        plt.plot(xdata, np.nanmean(data2, axis=1), label='Mid')
        plt.plot(xdata, np.nanmean(data3, axis=1), label='Top')
        plt.xlabel("Date [-]")
        plt.ylabel("Température [°C]")
        plt.legend()
        plt.show()


def plot_plates():
    plt.rcParams['date.converter'] = 'concise'
    plates = [0 for i in range(6)]
    s_range = {1:(1,5), 2:(5,9), 3:(9,13), 4:(13,19), 5:(19,23), 6:(23,29)}
    for i in range(1,7):
        plates[i-1] = np.nansum(np.dstack((data.loc[:, f"T[degC]-Low-S{s_range[i][0]}":f"T[degC]-Low-S{s_range[i][1]}"], data.loc[:, f"T[degC]-Mid-S{s_range[i][0]}":f"T[degC]-Mid-S{s_range[i][1]}"])),2)
        plates[i-1] = np.nansum(np.dstack((plates[i-1], data.loc[:, f"T[degC]-Top-S{s_range[i][0]}":f"T[degC]-Top-S{s_range[i][1]}"])),2) / 3
        plt.plot(time, np.nanmean(plates[i-1], axis=1), label=f"Plate {i}")
    plt.legend()
    plt.xlabel("Date [-]")
    plt.ylabel("Température moyenne d'une région [°C]")


def temp_evolution():
    plot_plates()
    plt.axvline(x=cut_time, color='black', linestyle='--')
    plt.show()


def control_rules():
    plot_plates()
    plt.axhline(y=0, color='black', linestyle='--')
    plt.axhline(y=3, color='red', linestyle='--')
    plt.show()




if __name__ == "__main__":
    check_stratification()
    temp_evolution()
    control_rules()
