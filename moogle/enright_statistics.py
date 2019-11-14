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
# Total times and plot
for file in files:
    dataname = file[:-8]
    table = pandas.read_csv(file, delimiter=', ', engine="python")

    time = table.loc[:, table.columns != 'total']
    time = time.loc[:, time.columns != 'redistance']
    total = time.sum(axis=1)
    df = pandas.DataFrame({dataname: total})
    total_df = pandas.concat([total_df, df], axis=1, sort=False)

    cellAdv = table[["cellAdvection"]]
    df = cellAdv.rename(columns={"cellAdvection": dataname})
    cellAdv_df = pandas.concat(
        [cellAdv_df, df], axis=1, sort=False)

sea.lineplot(data=total_df, dashes=False, ax=axes[0])
axes[0].set(ylabel="seconds", xlabel="timesteps")
axes[0].set_title("Total time")

sea.lineplot(data=cellAdv_df, dashes=False, ax=axes[1])
axes[1].set(ylabel="seconds", xlabel="timesteps")
axes[1].set_title("Cell advection")

plt.show()
