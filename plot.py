from numpy import loadtxt, where
from seaborn import plt
import seaborn as sns

sns.set_style('whitegrid', {'font.sans-serif': "Helvetica"})
sns.set_palette(sns.color_palette('Set1'))

# Number of variables to show
n_plots = 5

f, axarr = plt.subplots(n_plots, sharex=True, figsize=(5,8))

# Load all data
truth = loadtxt('truth.txt')
obs = loadtxt('obs.txt')
first_guess = loadtxt('first_guess.txt')
final_guess = loadtxt('final_guess.txt')

# Plot all data
for i in range(n_plots):
    axarr[i].plot(truth[:,0], truth[:,i+1], label="truth")
    axarr[i].plot(obs[:,0], obs[:,i+1], 'x', label="observations", alpha=1.0, mew=1)
    axarr[i].plot(first_guess[:,0], first_guess[:,i+1], '--', label="first guess", alpha=1.0)
    axarr[i].plot(final_guess[:,0], final_guess[:,i+1], ':', label="final guess")

plt.tight_layout()
leg = plt.figlegend(axarr[n_plots-1].get_lines(), ['truth', 'observations', 'first guess', 'final guess'], loc=(0.025, 0.0), ncol=4)
plt.xlabel('Time (nondimensional units)')
f.subplots_adjust(bottom=0.08)
plt.savefig('vars.pdf', bbox_extra_artists=(leg,), bbox_inches='tight')

# Plot diagnostics
plt.figure(figsize=(6,4))
diagn = loadtxt('diagnostics.txt')
length = where(diagn==0)[0]
length = diagn.shape[0] if len(length) == 0 else length[0]
plt.semilogy(diagn[:length-1,0], diagn[:length-1,1])
plt.xlabel('Iterations')
plt.ylabel('Cost function')

plt.tight_layout()
plt.savefig('cost_function.pdf', bbox_inches='tight')
plt.show()
