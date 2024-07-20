import os
import subprocess
from os import sys

start_num = int(sys.argv[1])
end_num = int(sys.argv[2])

indir = 'data/split/'

fis = [indir + fi for fi in sorted(os.listdir(indir))]

for i in range(start_num,end_num):
    subprocess.run(['python3', 'process_primary_tracks.py', fis[i]])


