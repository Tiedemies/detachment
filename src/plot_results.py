# FILEPATH: Untitled-1

import os
import pandas as pd

#The script should generate a pandas table with the following
# The name of the file is of the form N_C_p..num. The columns of the table should be N, C, p
# and also as each file contains four values, those should be as columns as well. These columns should be
# called: cut_epoi, base_epoi1, greedy_epoi, base_epoi2

# The script should also generate a plot with the following characteristics:
# The x axis should be N and there should a line for each value of P and also three lines,
# one for base epoi which is the average of base_epoi1 and base_epoi2, one for cut_epoi and one for greedy_epoi.

table = pd.DataFrame(columns = ["N", "C", "p", "cut_epoi", "base_epoi", "greedy_epoi"])

directory = "../results/"
for filename in os.listdir(directory):
  if filename.endswith(".num"):
    with open(os.path.join(directory, filename), "r") as file:
      lines = file.readlines()
      try:
        N, C, rest = filename.split("_")
        N = int(N)
        C = int(C)
      except:
        print(lines)
        continue
      p = rest.split("..")[0]
      cut,cut_epoi, base_epoi1, greedy_epoi, base_epoi2 = lines[0].split(",")
      
      table = table.append({"N": int(N), "C": int(C), "p": float(p), "cut_epoi": float(cut_epoi), "base_epoi": (float(base_epoi1)+float(base_epoi2))/2, "greedy_epoi": float(greedy_epoi)}, ignore_index=True)


import matplotlib.pyplot as plt
import numpy as np

print(table)

# With N as the x axis and "epoi" as the y-axis, plot cut_epoi's as red dots, base_epois as blue dots and greedy_epois as green dots. 
# The plot should have a legend and the axes should be labeled.

fig, ax = plt.subplots()
ax.plot(table["N"], table["cut_epoi"], 'ro', label='cut_epoi')
ax.plot(table["N"], table["base_epoi"], 'bo', label='base_epoi')
ax.plot(table["N"], table["greedy_epoi"], 'go', label='greedy_epoi')
# plot a line for each value of p, average over all C
for p in table["p"].unique():
  p_table = table[table["p"] == p]
  p_table = p_table.groupby("N").mean()
  ax.plot(p_table.index, p_table["cut_epoi"], 'r--')
  ax.plot(p_table.index, p_table["base_epoi"], 'b--')
  ax.plot(p_table.index, p_table["greedy_epoi"], 'g--')
ax.legend()
ax.set_xlabel('N')
ax.set_ylabel('epoi')
plt.show()



