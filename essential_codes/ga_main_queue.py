#!/home/djshin/venv/bin/python
from random import random
from ase.io import write
from ase.ga.data import DataConnection
from ase.ga.population import Population
from ase.ga.standard_comparators import InteratomicDistanceComparator
from ase.ga.cutandsplicepairing import CutAndSplicePairing
from ase.ga.offspring_creator import OperationSelector
from ase.ga.standardmutations import MirrorMutation
from ase.ga.standardmutations import RattleMutation
from ase.ga.standardmutations import PermutationMutation
from ase.ga.utilities import closest_distances_generator
from ase.ga.utilities import get_all_atom_types
from ase.ga.pbs_queue_run import PBSQueueRun
from os import system
from ase.io.trajectory import Trajectory

# GA script for supported metal structure tailored by D. Shin based on ASE 
# version: 200214
# Dongjae Shin @POSTECH

def jtg(job_name, traj_file):
    s = '#!/bin/bash\n'
    s += '#SBATCH --nodes=1\n'
    s += '#SBATCH --ntasks-per-node=20\n'
    s += '#SBATCH --partition=g3\n'
    s += '##\n'
    s += '#SBATCH --job-name="{0}"\n'.format(job_name)
    s += '#SBATCH --time=7-12:30\n'
    s += '#SBATCH -o STDOUT.%N.%j.out\n'
    s += '#SBATCH -e STDERR.%N.%j.err\n'
    s += '\n'
    s += '. /etc/profile.d/TMI.sh\n'
    s += '\n'
    s += 'python calc.py {0}\n'.format(traj_file)
    return s

# Parameters to be defined by a user -----------------------------------------
directory_name       = 'Pd8_CeO2_111_GA' # for termination of crontab
job_prefix           = 'Pd8-CeO2_111_opt' # for job names
n_converge           = 70
population_size      = 20
mutation_probability = 0.3
n_simul              = 11
# ----------------------------------------------------------------------------

# Initialize the different components of the GA
da = DataConnection('gadb.db')
tmp_folder = 'tmp_folder/'
# The PBS queing interface is created
pbs_run = PBSQueueRun(da,
                      tmp_folder=tmp_folder,
                      job_prefix=job_prefix,
                      n_simul=n_simul,
                      job_template_generator=jtg)

atom_numbers_to_optimize = da.get_atom_numbers_to_optimize()
n_to_optimize = len(atom_numbers_to_optimize)
slab = da.get_slab()
all_atom_types = get_all_atom_types(slab, atom_numbers_to_optimize)
blmin = closest_distances_generator(all_atom_types,
                                    ratio_of_covalent_radii=0.7)

comp = InteratomicDistanceComparator(n_top=n_to_optimize,
                                     pair_cor_cum_diff=0.015,
                                     pair_cor_max=0.7,
                                     dE=0.02,
                                     mic=False)
pairing = CutAndSplicePairing(slab, n_to_optimize, blmin)
mutations = OperationSelector([1., 1., 0.],
                              [MirrorMutation(blmin, n_to_optimize),
                               RattleMutation(blmin, n_to_optimize),
                               PermutationMutation(n_to_optimize)])
                              # no permutation (not alloy)

# Relax all unrelaxed structures (e.g. the starting population)
while (da.get_number_of_unrelaxed_candidates() > 0 and
       not pbs_run.enough_jobs_running()):
    a = da.get_an_unrelaxed_candidate()
    pbs_run.relax(a)

# create the population
population = Population(data_connection=da,
                        population_size=population_size,
                        comparator=comp)

# Submit new candidates until enough are running
while (not pbs_run.enough_jobs_running() and
       len(population.get_current_population()) > 2):
    a1, a2 = population.get_two_candidates()
    a3, desc = pairing.get_new_individual([a1, a2])
    if a3 is None:
        continue
    da.add_unrelaxed_candidate(a3, description=desc)

    if random() < mutation_probability:
        a3_mut, desc = mutations.get_new_individual([a3])
        if a3_mut is not None:
            da.add_unrelaxed_step(a3_mut, desc)
            a3 = a3_mut
    pbs_run.relax(a3)

    # Convergence Check whenever updated (by DJ)
    # Read candidates with energies from the trajectory file
    calculated_population = Trajectory('all_candidates.traj')
    N = len(calculated_population)

    # Make list of dictionary of gaids
    list = [{'gaid': calculated_population[i].info['key_value_pairs']['gaid']} 
              for i in range(N)] # to be sorted; for now sorted for energy
    gaid_Emin = list[0]['gaid']

    # Sort the list in the order of gaid (or time)
    list.sort(key=lambda x : x['gaid'])

    # Show the position of the GM in terms of gaid order
    for i in range(len(list)):
        if list[i]['gaid'] == gaid_Emin:
            i_max = i+1
            break

    print "GA: Current GM @ {0:5d} among total {1:5d} candidates".format(i_max, N)
    print "GA: Putative GM not changed for {0} cycles".format(N-i_max)

write('all_candidates.traj', da.get_all_relaxed_candidates())

# If GM is unchaged for over n_converge, delete GA in crontab jobs
if 'N' in locals() and 'i_max' in locals():  
    if (N-i_max) >= n_converge:
        print "GA: Convergence achieved!"
        print "GA: Terminating GA search ... (deleting crontab job)"
        system("crontab -l | grep -v {0}| crontab -".format(directory_name))
