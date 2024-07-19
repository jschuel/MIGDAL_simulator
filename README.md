# MIGDAL_simulator
Tools for convenient degrad simulation for MIGDAL. Plan to add in a quick digitizer too.

# Getting started

1. Clone this repository `git clone git@github.com:jschuel/degrad_tools.git`
2. In the `degrad_tools/` directory download the latest version of degrad
   ```sh
   wget https://degrad.web.cern.ch/degrad/degrad-3.19.f
   ```
3. Compile degrad with
   ```sh
   gfortran degrad-3.19.f -o degrad
   ````
   this will create an executable called `degrad`.
4. Adjust the number of electron recoil tracks and energy of the tracks you want to simulate in `configuration.yaml`. **Note:** Eventually it would be nice to replace all contents of the card file with this configuration file. Currently I don't have gas properties here, so you'll have to manually edit the `.card` file to adjust things like gas mixture, pressure, electric field strength. As it stands the .card file has 100% CF4. The comments in `degrad-3.19.f` explain the parameters of the `.card` file if you want to manually edit it yourself
5. Run `python3 run_and_process_degrad.py`. This script (1) automatically edits the `.card` file with the parameters of `configuration.yaml`, (2) runs degrad, (3) processes degrad's outputs to a track-indexed pandas dataframe, and (4) removes degrad's raw output.
6. Analyze or further process the output `.feather` file as you please!
