import seaborn as sea
import matplotlib.pyplot as plt
import csv
import os
import glob
import pandas
import numpy
sea.set(palette="muted", style="whitegrid", font_scale=0.75)
sea.despine()

folder = "../results/statistics/enright128/"

os.chdir(folder)
files = glob.glob("*.csv")

f, axes = plt.subplots(1, 2, figsize=(12, 6))

total_df = pandas.DataFrame()
cellAdv_df = pandas.DataFrame()
global_df = pandas.DataFrame()
# Total times and plot
for file in files:
    dataname = file[:-8]
    table = pandas.read_csv(file, delimiter=', ', engine="python")

    time = table.loc[:, table.columns != 'extractSurface']
    time = time.loc[:, time.columns != 'total']
#     time = time.loc[:, time.columns != 'redistance']

    total = table.loc[:, table.columns == 'total']
    df = total.rename(columns={"total": dataname})
    global_df = pandas.concat([global_df, df], axis=1, sort=False)

    total = time.sum(axis=1)
    df = pandas.DataFrame({dataname: total})
    total_df = pandas.concat([total_df, df], axis=1, sort=False)

    cellAdv = table[["cellAdvection"]]
    df = cellAdv.rename(columns={"cellAdvection": dataname})
    cellAdv_df = pandas.concat(
        [cellAdv_df, df], axis=1, sort=False)

print(global_df.head())

axisCount = 0

sea.lineplot(data=global_df, dashes=False, ax=axes[axisCount])
axes[axisCount].set(ylabel="seconds", xlabel="Frames")
axes[axisCount].set_title("Total time")
axisCount += 1

# sea.lineplot(data=total_df, dashes=False, ax=axes[axisCount])
# axes[axisCount].set(ylabel="seconds", xlabel="Frames")
# axes[axisCount].set_title("Average time")
# axisCount+=1;

sea.lineplot(data=cellAdv_df, dashes=False, ax=axes[axisCount])
axes[axisCount].set(ylabel="seconds", xlabel="timesteps")
axes[axisCount].set_title("Cell advection average")

plt.show()
