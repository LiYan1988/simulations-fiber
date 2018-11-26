import os
from subprocess import call
from shutil import copyfile

for file in os.listdir('.'):
    try:
        extension = file.split('.')[1]
        if extension=='sh':
            call(['sbatch', file])
    except:
        pass