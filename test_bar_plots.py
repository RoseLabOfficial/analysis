import matplotlib.pyplot as plt
import numpy as np

# Sample data
N = 5
menMeans = (20, 35, 30, 35, 27)
womenMeans = (25, 32, 34, 20, 25)

ind = np.arange(N)    # the x locations for the groups
men_width = 0.35       # the width of the men's bars
women_width = 0.5     # the width of the women's bars: intentionally larger

fig, ax = plt.subplots()

# Plotting women's bars, starting from the same y-position but with larger width
p2 = ax.bar(ind, womenMeans, women_width, label='Women', align='center', alpha=0.5)

# Plotting men's bars
p1 = ax.bar(ind, menMeans, men_width, label='Men', align='center')



ax.axhline(0, color='grey', linewidth=0.8)
ax.set_ylabel('Scores')
ax.set_title('Scores by group and gender')
ax.set_xticks(ind)
ax.set_xticklabels(('G1', 'G2', 'G3', 'G4', 'G5'))
ax.legend()

plt.show()
