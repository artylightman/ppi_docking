import glob
import os

a = glob.glob("*result")

for item in a:
    command =  "cp " + item + "/setVoronoiVolumes.dat volume_files/" + item[:-15] + "_volume10.txt"
    os.system(command)
