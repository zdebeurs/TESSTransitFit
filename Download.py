import numpy as np
import matplotlib.pyplot as plt

#pip install lightkurve
import lightkurve as lk

#Downloading the
#search_result = lk.search_lightcurve('HD189733')
#search_result.download_all(download_dir="data")

search_result = lk.search_lightcurve('HD147506')
search_result.download_all(download_dir="data")
