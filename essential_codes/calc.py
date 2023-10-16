#from ase.optimize import BFGS
from ase.io import read, write
#from ase.calculators.emt import EMT
#from ase.ga.relax_attaches import VariansBreak
import sys

from ase.calculators.vasp import Vasp

fname = sys.argv[1]

print('Now relaxing {0}'.format(fname))
a = read(fname)
# fname = 'tmp_folder//*.traj'

#a.set_calculator(EMT())
calc = Vasp(xc='pbe',
            setups={'Ce': '_h'},
            istart=0,
            icharg=1,
            ispin=2,
            idipol=3,
            encut=300.,
            prec='Normal',
            ediff=1.e-4,
            gga='PE',
            ismear=0,
            sigma=0.1,
            lreal=False,
            isym=1,
            nsw=800,
            ibrion=2,
            potim=0.1,
            isif=2,
            ediffg=-0.05,
            lmaxmix=6,
            ldau=True,
            ldau_luj={'Ce': {'L': 3, 'U': 6, 'J': 1}, 
                      'O' : {'L':-1, 'U': 0, 'J': 0},
                      'Pd': {'L':-1, 'U': 0, 'J': 0}},
            ldauprint=1,
            lwave=False,
            lcharg=False)

a.set_calculator(calc)

#dyn = BFGS(a, trajectory=None, logfile=None)
#vb = VariansBreak(a, dyn)
#dyn.attach(vb.write)
#dyn.run(fmax=0.05)

# Calculation in an independent directory
# not to overwrite VASP input files
import os, datetime
now = datetime.datetime.now().strftime("%Y_%m_%d-%H%M%S_%f")
os.system('mkdir {0}'.format(now))
os.chdir('{0}'.format(now))

## for old versioned run_vasp.py
#if os.path.exists(path):
#    os.system('cp -rf ../run_slurm_vasp.sh .')
#else:
#    print('No run_slurm_vasp.sh in the current directory!')
#    exit()

a.info['key_value_pairs']['raw_score'] = -a.get_potential_energy()

os.chdir('..')

write(fname[:-5] + '_done.traj', a)
# write to 'tmp_folder//*_done.traj'

print('Done relaxing {0}'.format(fname))
