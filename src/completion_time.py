import datetime as dt
import glob
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter, NullFormatter

dir = "/home/lyh970817/Data/GLAD/data_raw/*.csv"
files = glob.glob(dir)
# opt_files = [f for f in files if "Sign-up" not in f ]

opt_files = [f for f in files]

dats_list = [None] * len(opt_files)
names_list = [None] * len(opt_files)
dat = pd.read_csv(opt_files[1])

list(dat.columns)
dat["StartDate"]
for i, file in enumerate(opt_files):
    dats_list[i] = pd.read_csv(file).iloc[2:, :]
    dats_list[i].iloc[:, 0] = pd.to_datetime(dats_list[i].iloc[:, 0])
    dats_list[i].iloc[:, 1] = pd.to_datetime(dats_list[i].iloc[:, 1])
    # Extract the second instance of [A-Z]{2,}, since the first instance could be \
    # to do with variablees not in the questionnaire
    names_list[i] = pd.Series(
        dats_list[i].columns) \
        .str.extract(r"([A-Z]{2,})") \
        .dropna().values[1][0]

time_list = [dat.iloc[:, 1] - dat.iloc[:, 0] for dat in dats_list]
time_dict = dict(zip(names_list, time_list))


def PlotSummary(questionnaire, time_dict):
    """
    Generates a plot for the completion time of `questionnaire`
    and compute descriptives (including the mode)
    """
    hist = (time_dict[questionnaire].dt.seconds / 60).values
    descrptives = pd.Series(hist).describe().astype(int)

    # The mode seems to be the most suitable summary statistic here.
    mode = pd.Series(hist).mode()[0:1]
    mode.index = ["mode"]

    # plt.hist(hist[hist <= np.quantile(hist, q = 0.8)])
    hist_plot = plt.hist(hist)
    plt.yscale("symlog")
    ax = plt.axes()
    ax.yaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
    ax.yaxis.set_minor_formatter(NullFormatter())
    plt.suptitle('Completion Time of' + ' ' + questionnaire)
    # plt.title("Only data smaller than or equal to the 80% quantile were plotted.")
    plt.title("y axis in log10 scale")
    plt.xlabel('Minutes')
    plt.ylabel('Frequency')
    plt.text(plt.xlim()[1] / 2, 100,
             descrptives.append(mode).astype(int).to_string())
    # plt.show()


for q in time_dict:
    PlotSummary(q, time_dict)
    fname = "./" + q
    plt.savefig(fname)
    plt.close()
