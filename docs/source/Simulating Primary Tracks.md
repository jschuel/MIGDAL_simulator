# Simulating Primary Tracks
This is a guide detailing how to simulate primary tracks with MIGDAL_simulator. At the highest level, the workflow of this package can be broken down into two steps:

1. Simulate primary tracks

2. Process these tracks with `process_primary_tracks.py`

This guide details step (1). Step (1) provides primary electronic recoils (ERs), nuclear recoils (NRs), and Migdal effect tracks (Migdals). These are the inputs to step (2) which the simulation of these primary tracks inthe MIGDAL detector. Details of step (2) are provided in the [MIGDAL Detector Simulation](https://migdal-simulator.readthedocs.io/en/latest/MIGDAL%20Detector%20Simulation.html) section of these documents.

## Simulating ERs with Degrad

If you have not yet compiled Degrad on your machine, please work through the [Getting Started](https://migdal-simulator.readthedocs.io/en/latest/Getting%20Started.html) section of these documents.

Nearly all of MIGDAL_Simulator's functionality is configurable with `configuration.yaml`. Currently, `configuration.yaml` only supports simulating ERs in 100% CF4 gas. We hope to update this soon.

Below are the `configuration.yaml` parameters specifically relevant to primary ER track generation:

1. Under the `Degrad_card` heading:

- **input_file** - The filepath of Degrad's input card file. `input.card` comes default and doesn't need to change

- **n_tracks** - The number of primary tracks to simulate

- **seed** - The pseudorandom seed input into Degrad. If you run degrad parallely, it doesn't matter what you put in here; the seed will be set by `run_degrad.py`

- **energy** - Recoil energy of the primary track ERs **in eV**

- **primary_track_output_dir** - directory you want Degrad primary tracks output file to be saved to

2. Under the `Sim_settings` heading:

- **randomize_primary_track_order** - If True, this will shuffle the order of produced tracks. This does not have any effects on the tracks themselves, just where they show up in the output dataframe.

- **rotate_tracks** - If true, this "isotropizes" the angular distribution of primary tracks. If false, all tracks will point in the +x direction

- **parallel** - Setting this to True splits the Degrad processing into chunks and then concatenates the output. For instance if we're generating 5000 electrons and we set parallel to True and then parallel chunks (below) to 500. This will run 500 instances (each with different pseudorandom seeds) each generating 10 ERs. The 500 output file-chunks will then be concatenated into a single file and the file-chunks will then be deleted.

- **parallel_chunks** - See previous parameter

Once you have specified all of your parameters in `configuration.yaml`, there are three scripts you can run:

### Script 1: run_degrad.py

Run this script with `python3 run_degrad.py`. This script will generate tracks following the specifications you put in `configuration.yaml`. The output file will be placed in the `primary_track_output_dir` as specified in `configuration.yaml`. **This is the recommended script to use when you're producing a monoenergetic sample of ER primary tracks.**

### Script 2: run_degrad_multiple_energies.py

This is a script that runs `run_degrad.py` multiple times as a subprocess. The script still inherits everything from `configuration.yaml`, however it can also modify the `'n_tracks'`, `'energy'` and `'parallel_chunks'` field of `configuration.yaml`. As currently written, the script is meant to loop over a list of energies and generate ERs over that energy range. **The recommended usage of this script is therefore to generate ERs at user-defined energy steps**. Changing the energy range of this script is **not** supported by `configuration.yaml`, so to update this, you'll have to open `run_degrad_multiple_energies.py` in your favorite text editor and modify the information shown in the code block below.

```python
if __name__ == "__main__":
    yaml_file_path = 'configuration.yaml'
    processing_script_path = 'run_degrad.py'

    # Modify these values as needed
    n_tracks = 400
    energies = [200 * i for i in range(10, 61)]  #These are in eV so this range is 1-12 keV
    parallel_chunks = 1

    for e in energies:
        read_and_modify_yaml(yaml_file_path, n_tracks, e, parallel_chunks)
        run_processing_script(processing_script_path)
```

The specific information to update is `n_tracks`, `energy` and `parallel_chunks`, though I've found this script to run fastest when using 1 parallel_chunk (i.e. not splitting Degrad into sub-subprocesses).

Once you've configured everything you can run this script with `python3 run_degrad_multiple_energies.py`.

As currently written, this script doesn't combine the output primary track files. If the primary track files generated from this script are the only files in their directory, you can use something like the code block below to easily combine them

```python
import pandas as pd
import os

#Combined dataframe
combined = pd.concat([pd.read_feather(fi) for fi in sorted(os.listdir())])
combined = combined.sample(frac=1) #shuffle the order of the tracks
combined.index = [i for i in range(0,len(combined))] #reset the indices
#save
combined.to_feather("combined_output.feather") #chage to whatever name you want
```

### Script 3: run_and_process_degrad.py (**not recommended**)

This script creates primary tracks and performs the entire MIGDAL detector simulation for these tracks in a single swoop. While this can be useful, I prefer the modular approach of producing primary tracks and then feeding those into the MIGDAL simulation script. If you do run this script, read the description of `configuration.yaml` parameters in the [MIGDAL Detector Simulation](https://migdal-simulator.readthedocs.io/en/latest/MIGDAL%20Detector%20Simulation.html) section of these documents to understand the relevant parameters to tweak.

## Simulating NRs with SRIM
(coming soon)

## Generating Migdal primary tracks

You can generate simulated Migdals (ER and NR primary tracks stitched at their vertices) by adjusting the subfields of the `Migdal_sim` header in `configuration.yaml` and then running

```sh
python3 create_migdals.py
```
