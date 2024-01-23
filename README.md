Mamba packages required:
 - enlighten
 - numpy
 - pandas
 - oct2py
 - matplotlib

## Setup octave
apt install `sudo apt-get install`:
 - octave
 - octave-signal
Create `~/.octaverc` file. Add line: `pkg load signal`. This will automatically load the signal package everytime you load octave.

Ensure your path contains the GAMMA bin files, this directory and the directory containing the az_proc_tdbp_gpu binary file.