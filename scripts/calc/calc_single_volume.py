# author : artemis lightman
# date created : jan 2, 2023
# last modified: jan 2, 2023

# command line arguments
    # fname - pdb filename (not including extension)

# dependencies:
    # arty.py

# input: pdb files
# output: calculates voronoi tessellation for provided pdb files

##############################################################################
# completes voronoi volume/neighbour calculation
##############################################################################

import os
from datetime import datetime
import argparse
import arty

dt_start = datetime.now()

parser = argparse.ArgumentParser(description = "this script completes the voronoi volume/neighbour calculation")

parser.add_argument("-fname", help = "pdb filename (not including extension)", required = True)

args = parser.parse_args()

os.chdir("pomelo_vor")

##### generate surface ########

os.system("module load MATLAB/2017b")

cmd = '''matlab -nosplash -nodisplay -nojvm -r "generate_surface('generating_surfaces/''' + args.fname + '''_complex')"'''

os.system(cmd)

os.chdir("vor_code/bin")

os.system("chmod 777 *")

os.chdir("../..")

####### run pomelo ###########

cmd = "vor_code/bin/pomelo -mode GENERIC -i generating_surfaces/" + args.fname + "_complex_param.lua -o vor_calc/" + args.fname + "_complex_result"

os.system(cmd)

os.chdir("../pomelo_vor_10")

##### generate surface for bigger box ########

os.system("module load MATLAB/2017b")

cmd = '''matlab -nosplash -nodisplay -nojvm -r "generate_surface('generating_surfaces/''' + args.fname + '''_complex')"'''

os.system(cmd)

os.chdir("vor_code/bin")

os.system("chmod 777 *")

os.chdir("../..")

####### run pomelo (again) ###########

cmd = "vor_code/bin/pomelo -mode GENERIC -i generating_surfaces/" + args.fname + "_complex_param.lua -o vor_calc/" + args.fname + "_complex_result"

os.system(cmd)

dt = datetime.now()

str_dt = dt.strftime("%H:%M:%S on %d %B, %Y")

str_runtime = arty.get_runtime(dt_start, dt)

print ("calc_single_volume.py finished at " + str_dt + " (runtime " + str_runtime + ")")