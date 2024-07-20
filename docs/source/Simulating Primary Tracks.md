#Simulating Primary Tracks
This is a guide detailing how to simulate primary tracks with MIGDAL_simulator. At the highest level, the workflow of this package can be broken down into two steps:

1. Simulate primary tracks
2. Process these tracks with `process_primary_tracks.py`

This guide details step (1). Step (1) provides primary electronic recoils (ERs), nuclear recoils (NRs), and Migdal effect tracks (Migdals). These are the inputs to step (2) which the simulation of these primary tracks inthe MIGDAL detector. Details of step (2) are provided in 