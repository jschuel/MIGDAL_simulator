'''JTS 06/26/2024: Wrapper code for running FORTRAN executables. This code adjusts the number of tracks, energy, and seed
values from the FORTRAN input card file. After running the executable, the code processes the output into a track-indexed
pandas dataframe that includes the number of ionization electrons, x,y, and z information of the ER (in um), and the timing 
information in ps. A flag is also included for  1, fluorescence; 2, pair production; 3, bremsstrahlung; 0, otherwise'''

import subprocess
import pandas as pd
import numpy as np
import os
from tqdm import tqdm
import time
from scipy.spatial import KDTree #Determines if charge passes through GEM holes FAST

class run_degrad:
    def __init__(self,card_file,num_events,energy,drift,v_drift,sigmaT,sigmaL,sigmaT_trans,sigmaL_trans,
                 sigmaT_induc,sigmaL_induc,W,GEM_width, GEM_height, hole_diameter, hole_pitch, amplify,gain,
                 truth_dir=np.array([1,0,0]),random_seed=0,rotate=False):
        self.card_file = card_file
        self.num_events = num_events
        self.seed = random_seed
        self.truth_dir = truth_dir
        self.energy = format(energy, '.5f')
        self.gain=gain
        self.W = W
        self.v_drift = v_drift
        self.sigmaT = sigmaT
        self.sigmaL = sigmaL
        self.GEM_width = GEM_width
        self.GEM_height = GEM_height
        self.hole_diameter = hole_diameter
        self.hole_pitch = hole_pitch
        self.outfile_name = '%s_events_%skeV_%s_seed.out'%(self.num_events,float(self.energy)/1000,self.seed)
        '''Modify num_events, random_seed, and energy in CARD file'''
        self.modify_file()

        '''Run executable and change the name of the output to outfile_name'''
        self.run_degrad_executable()

        self.tracks = self.convert_output_to_numpy()
        '''Group tracks into a track-indexed dataframe'''
        self.output = self.convert_numpy_to_pandas()

        '''Isotropize track direction if rotate is set to True'''
        if rotate:
            self.output = self.rotate_tracks(self.output)
        else:
            self.output['truth_dir'] = [self.truth_dir for i in range(0,len(self.output))]

        self.output['truth_energy'] = float(self.energy)
        self.output['ionization_energy'] = self.output['nHits']*self.W / 1000 #keV
        self.output['truth_theta'] = self.output['truth_dir'].apply(lambda x: np.arccos(x[2]))
        self.output['truth_phi'] = self.output['truth_dir'].apply(lambda x: np.arctan2(x[1],x[0]))

        '''Center tracks (make vertex the origin in z)'''
        self.output = self.center_tracks(self.output)
        '''Apply drift if specified'''
        if drift:
            self.output['drift_length'] = np.random.uniform(1,3,len(self.output)) #cm
            self.output = self.apply_diffusion_all_tracks(self.output)
            
        '''Amplify if specified'''
        if amplify:
            self.hole_positions = self.create_GEM_holes()
            #Create KD tree
            self.hole_tree = KDTree(self.hole_positions)
            xamps = []
            yamps = []
            zamps = []
            xgains = []
            ygains = []
            zgains = []
            for i in range(len(self.output)):
                event = self.output.iloc[i]
                track = np.array([event['xdiff'],event['ydiff'],event['zdiff']]).T
                charges_passing_through = np.array([charge for charge in track if self.is_within_GEMhole(charge[0], charge[1])])
                xamp = charges_passing_through[:,0]
                yamp = charges_passing_through[:,1]
                zamp = charges_passing_through[:,2]
                xamps.append(xamp)
                yamps.append(yamp)
                zamps.append(zamp)
                '''Apply first amplification'''
                xgain,ygain,zgain = self.GEM_gain_and_diffusion(xamp,yamp,zamp,gap='transfer')
                xgains.append(xgain)
                ygains.append(ygain)
                zgains.append(zgain)
                
            self.output['xamp'] = xamps
            self.output['yamp'] = yamps
            self.output['zamp'] = zamps

            self.output['xgain'] = xgains
            self.output['ygain'] = ygains
            self.output['zgain'] = zgains

        '''Save output dataframe'''
        if rotate:
            outname = os.path.splitext(self.outfile_name)[0]+'_isotropic.feather'
        else:
            outname = os.path.splitext(self.outfile_name)[0]+'.feather'
        if not os.path.exists('data'):
            os.makedirs('data')
        self.output.to_feather('data/'+outname)

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
            with open(os.devnull,'w') as devnull:
                subprocess.run(['./degrad'],stdin=infile,stdout=devnull,stderr=devnull)

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
        df['x']     = self.tracks[:,1] / 10000 #cm
        df['y']     = self.tracks[:,2] / 10000 #cm
        df['z']     = self.tracks[:,3] / 10000 #cm
        df['t']     = self.tracks[:,4] #ps
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

    '''Isotropizes track angular distributions'''
    def rotate_track(self,track,init_dir):
        
        def random_theta_phi(): #get random theta and phis for rotation
            ctheta = np.random.uniform(-1,1) #draw from uniform cos(theta) distribution
            phi = np.random.uniform(0,2*np.pi)
            theta = np.arccos(ctheta)
            x = np.sin(theta)*np.cos(phi)
            y = np.sin(theta)*np.sin(phi)
            z = np.cos(theta)
            return theta, np.arctan2(y,x)

        def rotate_y(x,y,z,angle): #rotate about y axis
            xp = np.cos(angle)*x+np.sin(angle)*z
            yp = y
            zp = -np.sin(angle)*x+np.cos(angle)*z
            return xp,yp,zp
    
        def rotate_z(x,y,z,angle): #rotate about z axis
            xp = np.cos(angle)*x-np.sin(angle)*y
            yp = np.sin(angle)*x+np.cos(angle)*y
            zp = z
            return xp,yp,zp

        theta,phi = random_theta_phi()

        '''Rotate tracks to make them directionally isotropic'''
        x_r1, y_r1, z_r1 = rotate_y(track['x'], track['y'], track['z'], -(np.pi/2-theta)) #rotate track coordinates about y axis
        x_r2, y_r2, z_r2 = rotate_z(x_r1, y_r1, z_r1, phi) #rotate track coordinates about z axis
        
        '''Do the same for truth directions'''
        dir_rx1, dir_ry1, dir_rz1 = rotate_y(init_dir[0], init_dir[1], init_dir[2], -(np.pi/2-theta)) #rotate reocil direction about y axis
        dir_rx2, dir_ry2, dir_rz2 = rotate_z(dir_rx1, dir_ry1, dir_rz1, phi) #rotate reocil direction about z axis
        
        dir_r = np.array([dir_rx2, dir_ry2, dir_rz2])
        
        return x_r2, y_r2, z_r2, dir_r

    def rotate_tracks(self,df):
        #print("\nRandomly rotating tracks\n")
        xrot   = []
        yrot   = []
        zrot   = []
        dirrot = []
        for i in range(0,len(df)):
            track = df.iloc[i]
            xr, yr, zr, dirr = self.rotate_track(track,init_dir = self.truth_dir) #Rotate track coordinates..also centers tracks
            xrot.append(xr)
            yrot.append(yr)
            zrot.append(zr)
            dirrot.append(dirr)
        df['x'] = xrot
        df['y'] = yrot
        df['z'] = zrot
        df['truth_dir'] = dirrot

        return df

    def center_tracks(self,df):
        df['x'] = df['x'].apply(lambda x: x-x.mean())
        df['y'] = df['y'].apply(lambda x: x-x.mean())
        df['z'] = df['z'].apply(lambda x: x-x[0])
        return df
    
    def apply_diffusion(self, track): #applies diffusion using input diffusion parameters
        xs = track['x']
        ys = track['y']
        zs = track['z']+track['drift_length']
        x_diff = np.sqrt(zs)*self.sigmaT*1e-4*np.random.normal(0,1, len(zs))
        y_diff = np.sqrt(zs)*self.sigmaT*1e-4*np.random.normal(0,1, len(zs))
        z_diff = np.sqrt(zs)*self.sigmaL*1e-4*np.random.normal(0,1, len(zs))
        xs = xs+x_diff
        ys = ys+y_diff
        zs = zs+z_diff
        
        del x_diff, y_diff, z_diff
        return xs,ys,zs

    def apply_diffusion_all_tracks(self,df):
        xdiff = []
        ydiff = []
        zdiff = []
        for i in range(0,len(df)):
            track = df.iloc[i]
            xd,yd,zd = self.apply_diffusion(track)
            xdiff.append(xd)
            ydiff.append(yd)
            zdiff.append(zd)
        df['xdiff'] = xdiff
        df['ydiff'] = ydiff
        df['zdiff'] = zdiff
        return df

    '''Creates honeycomb grid of GEM holes with user-input diameter and pitch'''
    def create_GEM_holes(self):
        # Convert micrometers to centimeters for plotting
        hole_diameter_cm = self.hole_diameter / 10000  # cm
        self.hole_radius_cm = hole_diameter_cm / 2
        hole_pitch_cm = self.hole_pitch / 10000  # cm

        # Calculate the number of holes in the x and y directions
        num_holes_x = int(self.GEM_width / hole_pitch_cm)
        num_holes_y = int(self.GEM_height / (hole_pitch_cm * np.sqrt(3) / 2))

        # Generate the positions of the holes in a honeycomb pattern
        hole_positions = []

        for i in range(0,num_holes_x):
            for j in range(0,num_holes_y):
                x = i * hole_pitch_cm
                y = j * (hole_pitch_cm * np.sqrt(3) / 2)
                if j % 2 == 1:
                    x += hole_pitch_cm / 2  # Offset every other row
                hole_positions.append((x, y))

        shift = self.GEM_width / 2
        pos = pd.Series(hole_positions).apply(np.array)-shift
        hole_positions = pos[pos.apply(lambda x: (x[0]>-1*shift) & (x[1]>-1*shift) & (x[0] < shift) & (x[1] < shift))].apply(tuple).to_list()
        return hole_positions

    def is_within_GEMhole(self,x, y):
        distance, index = self.hole_tree.query([x, y])
        return distance <= self.hole_radius_cm

    def generate_gain_points(self, x, x_post, gain_electrons, gap_length, diff_coeff): #Generates x, y, and z coordiantes after gain
        for enum, val in np.ndenumerate(gain_electrons):
            start_ind = np.sum(gain_electrons[:enum[0]])
            end_ind = np.sum(gain_electrons[:enum[0]+1])
            x_post[start_ind:end_ind] = x[enum] + np.sqrt(gap_length)*diff_coeff*1E-4*np.random.normal(0,1,val)

    def GEM_gain_and_diffusion(self,x,y,z,gap): #Applies gain and readout resolution smearing
        gain_electrons = np.random.exponential(np.sqrt(self.gain), len(x))
        gain_electrons = np.asarray(gain_electrons, dtype=int)
        
        x_post = np.ascontiguousarray(np.zeros(np.sum(gain_electrons)),dtype=np.float32)
        y_post = np.ascontiguousarray(np.zeros(np.sum(gain_electrons)),dtype=np.float32)
        z_post = np.ascontiguousarray(np.zeros(np.sum(gain_electrons)),dtype=np.float32)

        if gap.lower() == 'transfer':
            self.generate_gain_points(x, x_post, gain_electrons, gap_length = 0.2, diff_coeff = sigmaT_trans)
            self.generate_gain_points(y, y_post, gain_electrons, gap_length = 0.2, diff_coeff = sigmaT_trans)
            self.generate_gain_points(z, z_post, gain_electrons, gap_length = 0.2, diff_coeff = sigmaT_trans)
            
        elif gap.lower() == 'induction':
            self.generate_gain_points(x, x_post, gain_electrons, gap_length = 0.2, diff_coeff = sigmaT_induc)
            self.generate_gain_points(y, y_post, gain_electrons, gap_length = 0.2, diff_coeff = sigmaT_induc)
            self.generate_gain_points(z, z_post, gain_electrons, gap_length = 0.2, diff_coeff = sigmaT_induc)

        else:
            raise ValueError("Must specify 'transfer' or 'induction' gap for simulating gain across the GEM")

        return x_post, y_post, z_post
    
if __name__ == '__main__':
    import yaml
    
    '''Load configuration.yaml'''
    with open('configuration.yaml','r') as cfg:
        config = yaml.safe_load(cfg)
    
    card = config['input_file']
    n_events = config['n_tracks']
    seed = config['seed']
    E = config['energy']
    rot = config['rotate_tracks']
    parallel = config['parallel']
    nchunks = config['parallel_chunks']
    W = config['W']
    drift = config['apply_drift']
    v_drift = config['vd']
    sigmaT = config['sigmaT']
    sigmaL = config['sigmaL']
    sigmaT_trans = config['sigmaT_trans']
    sigmaL_trans = config['sigmaL_trans']
    sigmaT_induc = config['sigmaT_induc']
    sigmaL_induc = config['sigmaL_induc']
    GEM_width = config['GEM_width']
    GEM_height = config['GEM_height']
    hole_diameter = config['hole_diameter']
    hole_pitch = config['hole_pitch']
    amplify = config['apply_amplification']
    gain = config['gain']

    if not parallel:
        '''Update card, run degrad, process output, save pandas dataframe as feather file'''
        run_degrad(card_file=card,
                   num_events=n_events,
                   energy = E,
                   random_seed = seed,
                   rotate = rot,
                   W = W,
                   drift = drift,
                   v_drift = v_drift,
                   sigmaT = sigmaT,
                   sigmaL = sigmaL,
                   sigmaT_trans = sigmaT_trans,
                   sigmaL_trans = sigmaL_trans,
                   sigmaT_induc = sigmaT_induc,
                   sigmaL_induc = sigmaL_induc,
                   GEM_width=GEM_width,
                   GEM_height=GEM_height,
                   hole_diameter=hole_diameter,
                   hole_pitch=hole_pitch,
                   amplify=amplify,
                   gain=gain)
        '''If we run in parallel chunks, we just set random_seed to the chunk number'''
    else:
        n_events = n_events//nchunks
        print('Running DEGRAD %s times'%(nchunks))
        start = time.time() #log start time
        for chunk in tqdm(range(0,nchunks)):
            run_degrad(card_file=card,
                       num_events=n_events,
                       energy = E,
                       random_seed = chunk,
                       rotate = rot,
                       W = W,
                       drift = drift,
                       v_drift = v_drift,
                       sigmaT = sigmaT,
                       sigmaL = sigmaL,
                       sigmaT_trans = sigmaT_trans,
                       sigmaL_trans = sigmaL_trans,
                       sigmaT_induc = sigmaT_induc,
                       sigmaL_induc = sigmaL_induc,
                       GEM_width=GEM_width,
                       GEM_height=GEM_height,
                       hole_diameter=hole_diameter,
                       hole_pitch=hole_pitch,
                       amplify=amplify,
                       gain=gain)

        '''Check all files that were created after start'''
        files = []
        outpath = 'data/'
        for filename in os.listdir(outpath):
            filepath = os.path.join(outpath, filename)
            if os.path.isfile(filepath):
                creation_time = os.path.getctime(filepath)
                if creation_time > start:
                    files.append(filename)

        print('Concatenating files created this run\n')
        df = pd.concat([pd.read_feather(outpath+fi) for fi in sorted(files)])
        df.index = [i for i in range(0,len(df))]
        df.to_feather(outpath+"%skeV_%sEvents_all.feather"%(float(E)/1000,n_events*nchunks))
        print('DONE\n')
        print('Cleaning intermediate files\n')
        for fi in files:
            os.remove(outpath+fi)

        print('DONE')
