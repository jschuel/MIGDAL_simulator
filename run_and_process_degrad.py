'''JTS 06/26/2024: Wrapper code for running FORTRAN executables. This code adjusts the number of tracks, energy, and seed
values from the FORTRAN input card file. After running the executable, the code processes the output into a track-indexed
pandas dataframe that includes the number of ionization electrons, x,y, and z information of the ER (in um), and the timing 
information in ps. A flag is also included for  1, fluorescence; 2, pair production; 3, bremsstrahlung; 0, otherwise'''

import subprocess
import pandas as pd
import numpy as np
import os

class run_degrad:
    def __init__(self,card_file,num_events,energy,random_seed=0):
        self.card_file = card_file
        self.num_events = num_events
        self.seed = random_seed
        self.energy = format(energy, '.5f')
        self.outfile_name = '%s_events_%skeV_%s_seed.out'%(self.num_events,float(self.energy)/1000,self.seed)
        '''Modify num_events, random_seed, and energy in CARD file'''
        self.modify_file()

        '''Run executable and change the name of the output to outfile_name'''
        self.run_degrad_executable()

        self.tracks = self.convert_output_to_numpy()
        '''Group tracks into a track-indexed dataframe'''
        self.output = self.convert_numpy_to_pandas()

        '''Save output dataframe'''
        outname = os.path.splitext(self.outfile_name)[0]+'.feather'
        self.output.to_feather(outname)

        '''Delete raw DEGRAD output'''
        self.remove_degrad_output()
        
    '''Modifies num_events, random_seed, and energy in CARD file'''
    def modify_file(self):
        # Read the input file
        with open(self.card_file, 'r') as file:
            lines = file.readlines()

        # Modify the specific numbers in the first line
        lines[0] = '        1      %s         2         0         %s    %s    2.0000    0.0000\n'%(self.num_events,self.seed,self.energy)
    
        # Write the modified content to the output file
        with open(self.card_file, 'w') as file:
            file.writelines(lines)

    def run_degrad_executable(self):
        # Run the Fortran executable
        with open(self.card_file, 'r') as infile:
            subprocess.run(['./degrad'],stdin=infile)

        # Rename the output file
        subprocess.run(['mv', 'DEGRAD.OUT', self.outfile_name])

    def convert_output_to_numpy(self):
        '''Modified version of Tim's conversion code'''
        tracks = []
        with open(self.outfile_name,"r") as f:
            for i, line in enumerate(f):
      
                if i%2==0:
                    num_es = int(line[21:42])
                    track = np.zeros((num_es,6))
        
                else:
                    index = i//2
                    track[:,0] = index
        
                    for e_num in range(num_es):
                        j = e_num * 167
                        track[e_num,1] = float(line[j:j+26])
                        track[e_num,2] = float(line[j+26:j+52])
                        track[e_num,3] = float(line[j+52:j+78])
                        track[e_num,4] = float(line[j+78:j+104])
          
                        # code to 0,1,2,3:
                        f = int(line[j+104:j+125])
                        pp = line[j+145:j+146]
                        b = line[j+166:j+167]
          
                        if f != 0:
                            track[e_num,5] = 1.
                        elif pp != "0":
                            track[e_num,5] = 2.
                        elif b != "0":
                            track[e_num,5] = 3.

                    track = track[np.argsort(track[:,4]),:] # sort by time
        
                    if abs(track[0,1]) > 25 or abs(track[0,2]) > 25 or abs(track[0,3]) > 25:
                        track[:,1:4] -= track[0:1,1:4] + np.random.normal(0,5,(1,3))
        
                    tracks.append(track)
  
        tracks = np.concatenate(tracks)
        return tracks

    def convert_numpy_to_pandas(self):

        df = pd.DataFrame()
        df['index'] = self.tracks[:,0]
        df['x']     = self.tracks[:,1]
        df['y']     = self.tracks[:,2]
        df['z']     = self.tracks[:,3]
        df['t']     = self.tracks[:,4]
        df['flag']  = self.tracks[:,5]

        grp = df.groupby('index').agg(list)
        grp.index = [i for i in range(0,len(grp))]
        for col in grp.columns:
            grp[col] = grp[col].apply(np.array)
        grp['flag'] = grp['flag'].apply(lambda x: x.astype('uint8'))
        grp['nHits'] = grp['x'].apply(lambda x: len(x))

        return grp[['nHits','x','y','z','t','flag']]

    def remove_degrad_output(self):
        os.remove(self.outfile_name)
    
if __name__ == '__main__':
    import yaml
    
    '''Load configuration.yaml'''
    with open('configuration.yaml','r') as cfg:
        config = yaml.safe_load(cfg)
    
    card = config['input_file']
    n_events = config['n_tracks']
    seed = config['seed']
    E = config['energy']

    '''Update card, run degrad, process output, save pandas dataframe as feather file'''
    run_degrad(card_file=card,
               num_events=n_events,
               energy = E,
               random_seed = seed)
