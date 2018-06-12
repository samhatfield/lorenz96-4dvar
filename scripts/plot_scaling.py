from seaborn import plt
import seaborn as sns
from numpy import load
from os import chdir

chdir('../')

sns.set_style('whitegrid', {'font.sans-serif': "Helvetica"})
sns.set_palette(sns.color_palette('Set1'))

colour_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

sbits = range(8,53)

data = load('incr_4dvar_performance_scaling.npz')['arr_0']

cost, n_iters = data[:,0], data[:,1]

fig, ax1 = plt.subplots(figsize=(5,4))
ax1.plot(sbits, cost, color=colour_cycle[0])
ax1.set_xlabel('Number of significand bits')
ax1.set_ylabel('Minimum obtained', color=colour_cycle[0])
ax2 = ax1.twinx()
ax2.plot(sbits, n_iters, color=colour_cycle[1])
ax2.set_ylabel('Number of iterations required', color=colour_cycle[1])
plt.savefig('incr_scaling_plot.pdf', bbox_inches='tight')
plt.show()
