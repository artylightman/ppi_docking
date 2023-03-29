import glob
import sys

a = glob.glob("*model*complex_H.pdb")

output = open("models.lst", "w")

chains = sys.argv[1]

for item in a:
    output.write(item + " " + chains + "\n")