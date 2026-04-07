import pandas as pd


df = pd.read_csv("ASTRO501/seminar/table_Science-Exposure.csv")
obslist = [2458963.7175579001,2458936.7646181001,2458900.8169328999,2459001.6712731002,2459037.7144559999,2458838.0099769002]
df_selected_images = df[df["obsjd"].isin(obslist)]
df_selected_images["mjd"] = df_selected_images["obsjd"]- 2400000.5
print(df_selected_images)

import matplotlib.pyplot as plt

plt.plot(df_selected_images["mjd"], df_selected_images["seeing"])
plt.show()