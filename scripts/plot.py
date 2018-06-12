from numpy import loadtxt, where
from seaborn import plt
import seaborn as sns

sns.set_style('whitegrid', {'font.sans-serif': "Helvetica"})
sns.set_palette(sns.color_palette('Set1'))

# Variables to show
display_vars = [3, 13, 25, 35]

f, axarr = plt.subplots(len(display_vars), sharex=True, figsize=(4,6))

# Load all data
truth = loadtxt('truth.txt')
obs = loadtxt('obs.txt')
first_guess = loadtxt('first_guess.txt')
final_guess = loadtxt('final_guess.txt')

# Plot all data
for i, var in enumerate(display_vars):
    axarr[i].plot(truth[:,0], truth[:,var+1], label="truth")
    axarr[i].plot(obs[:,0], obs[:,var+1], 'x', label="observations", alpha=1.0, mew=1)
    axarr[i].plot(first_guess[:,0], first_guess[:,var+1], '--', label="first guess", alpha=1.0)
    axarr[i].plot(final_guess[:,0], final_guess[:,var+1], ':', label="final guess")

plt.tight_layout()
leg = plt.figlegend(axarr[len(display_vars)-1].get_lines(), ['truth', 'observations', 'first guess', 'final guess'], loc='lower center', ncol=2)
plt.xlabel('Time (nondimensional units)')
f.subplots_adjust(bottom=0.17)
plt.savefig('vars.pdf', bbox_extra_artists=(leg,), bbox_inches='tight')

# Plot diagnostics
plt.figure(figsize=(6,3))
diagn = loadtxt('diagnostics.txt')
length = where(diagn<=0)[0]
length = diagn.shape[0] if len(length) == 0 else length[0]
plt.semilogy(diagn[:length-1,0], diagn[:length-1,1])
plt.xlabel('Iterations')
plt.ylabel('Cost function')

plt.tight_layout()
plt.savefig('cost_function.pdf', bbox_inches='tight')
plt.show()
