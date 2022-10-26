import numpy as np
import corner #pip install corner
import matplotlib.pyplot as plt

Parameters = ["Period", "T0", "RpRs", "aRs", "b", "u1", "u2"]


Data = np.load("Figures/HD189733Sample.npy")
print(np.shape(Data))

X, Y, Z = np.shape(Data)

Data = Data[:,Y//2:,:]
X, Y, Z = np.shape(Data)

Data = Data.reshape(X*Y, Z)

#Maker corner plot
corner.corner(Data,plot_contours="True", labels=Parameters, title_fmt=".2e",quantiles=[0.16, 0.5, 0.84], show_titles=True, title_kwargs={"fontsize": 12})
plt.savefig("Figures/HD189733b_CornerPlot.png")


for i in range(Z):
    Lower, Median, Upper = np.percentile(Data[:, i], [16., 50, 84.])

    ErrorUpper = str(Upper-Median)
    ErrorLower = str(Median-Lower)
    Text = str(Median)+"$^{"+ErrorUpper+"}"+"_{"+ErrorLower+"}$"
    print(Parameters[i], Text)
