from subprocess import check_output
from fileinput import input
from numpy import std, mean
from re import sub
from numpy import savez
from os import chdir

chdir('../')

# Run for each precision
results = []
for j in range(8,53):
    # Splice parameters
    with open("params.template.f90", "r") as template:
        lines = template.readlines()
    with open("params.f90", "w") as params:
        for line in lines:
            params.write(sub('SBITS', str(j), line))

    # Build code
    check_output(['make', '-s', 'clean'])
    check_output(['make', '-s'])

    cost, iters = tuple(float(f) for f in check_output(['./main']).split())
    results.append([cost, iters])
    print((cost, iters))

savez('scripts/incr_4dvar_performance_scaling.npz', results)
