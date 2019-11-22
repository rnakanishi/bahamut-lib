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
    dataname = file[6:-4]
    table = pandas.read_csv(file, delimiter=', ', engine="python")

    bandCells = table.loc[:, table.columns == 'bandCells']
    fluidCells = table.loc[:, table.columns == 'fluidCells']
    cfl = table.loc[:, table.columns == 'cflSteps']

    total = table.loc[:, table.columns == 'total']

    redistance = table.loc[:, table.columns == 'redistance']
    redistance = redistance.rename(
        columns={'redistance': 'total'})

    surfaceExtraction = table.loc[:, table.columns == 'extractSurface']
    surfaceExtraction = surfaceExtraction.rename(
        columns={'extractSurface': 'total'})

    # Removing some statistics that may interfere on real time computation
    total = total - surfaceExtraction
    total = total - redistance  # May count depending on situation

    df = total.rename(columns={"total": dataname})
    global_df = pandas.concat([global_df, df], axis=1, sort=False)

    # Each time substeps take to advect cells
    cellAdv = table[["cellAdvection"]]
    cfl = cfl.rename(columns={'cflSteps': 'cellAdvection'})
    cellAdv = cellAdv/cfl
    df = cellAdv.rename(columns={"cellAdvection": dataname})
    cellAdv_df = pandas.concat(
        [cellAdv_df, df], axis=1, sort=False)

axisCount = 0
sea.lineplot(data=global_df, dashes=False, ax=axes[axisCount])
axes[axisCount].set(ylabel="Seconds", xlabel="Frames")
axes[axisCount].set_title("Total time")
axisCount += 1

# sea.lineplot(data=total_df, dashes=False, ax=axes[axisCount])
# axes[axisCount].set(ylabel="Seconds", xlabel="Frames")
# axes[axisCount].set_title("Average time")
# axisCount+=1;

sea.lineplot(data=cellAdv_df, dashes=False, ax=axes[axisCount])
axes[axisCount].set(ylabel="Seconds", xlabel="Frames")
axes[axisCount].set_title("Cell advection average\n(substep time)")

plt.show()
