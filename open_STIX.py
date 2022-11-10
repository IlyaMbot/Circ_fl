import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.dates as mdates
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)


names = ["time", "4-10keV", "10-15keV", "15-25keV", "25-50keV", "50-84keV"]

df = pd.read_csv("STIX/stix_quick_look_light_curves.csv", header = 0, names = names, sep = r"\,")
df = df.sort_values(by = "time")

time = np.array(df["time"][160:350], dtype = "datetime64")

jet_time, jet_int = np.load("./jet.npy", allow_pickle = True)
flux_time, flux_int = np.load("./flare.npy", allow_pickle = True)

jet_time = np.array(jet_time, dtype= "datetime64")
flux_time = np.array(flux_time, dtype= "datetime64")

time_str = [str(t) for t in time]

size = 16
myFmt = mdates.DateFormatter('%H:%M:%S')

fig, ax1 = plt.subplots(figsize = (7, 6))
ax1.set_title("STIX counts & SRH I profiles 22 May 2021 ", size = size * 1.1, weight = "bold")
ax1.set_xlabel("Time, UT", size = size, weight = "bold")
ax1.set_ylabel("Counts, [log]", size = size, weight = "bold")
ax1.set_yscale("log")
ax1.plot(time, df[names[1]][160:350], label = names[1], lw = 2, color = "C1" )
ax1.plot(time, df[names[3]][160:350], label = names[3], lw = 2, color = "C2" )

ax2 = ax1.twinx()
ax2.plot(jet_time, jet_int / np.max(jet_int), color = "Red", lw = 3, label = "Jet")
ax2.plot(flux_time, flux_int / np.max(flux_int), color = "blue", lw = 3, label = "Jet")
ax2.set_ylim(0.3,1.01)
ax2.xaxis.set_major_formatter(myFmt)
plt.ylabel("Intensity, [normalized]", size = size, weight = "bold")



ax1.legend()
plt.tight_layout()
plt.savefig(f"./STIX_and_SRH_I.png", transparent = False, dpi = 400, bbox_inches = "tight")
# plt.show()
