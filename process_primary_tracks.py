'''Loads in an input file processed with Jeff's fork of RCTRIM and then simulates, drift, diffusion, amplification and digitization'''

import pandas as pd
import numpy as np
import os
from tqdm import tqdm
import time
from scipy.spatial import KDTree #Determines if charge passes through GEM holes FAST

class process_tracks:
    def __init__(self,infile,outpath,nGEM,drift,min_drift_length,max_drift_length,v_drift,sigmaT,
                 sigmaL, sigmaT_trans,sigmaL_trans, sigmaT_induc,sigmaL_induc, drift_gap_length,
                 GEM_width, GEM_height, GEM_thickness, hole_diameter, hole_pitch, extra_GEM_diffusion, amplify,gain,
                 transfer_gap_length, induction_gap_length, GEM_offsetsx, GEM_offsetsy,
                 randomize_position, write_gain = False, overwrite = False, use_gpu = False):

        self.gpu = use_gpu
        self.nGEM = int(nGEM)
        self.gain=gain**(1/self.nGEM) #scale of random exponential for amplification.
        self.drift_gap_length = drift_gap_length
        self.randomize_position = randomize_position
        self.v_drift = v_drift
        self.sigmaT = sigmaT
        self.sigmaL = sigmaL
        self.GEM_width = GEM_width
        self.GEM_height = GEM_height
        self.hole_diameter = hole_diameter
        self.hole_pitch = hole_pitch
        self.GEM_offsetsx = GEM_offsetsx
        self.GEM_offsetsy = GEM_offsetsy
        self.extra_GEM_diffusion = extra_GEM_diffusion
        
        print(infile)

        self.data = pd.read_feather(infile)
        if 'ID' in self.data.columns: #autocheck if the event is a Migdal
            self.migdal = True
        else:
            self.migdal = False

        '''Randomize position if specified'''
        if self.randomize_position:
            self.data['xshift'] = np.random.uniform(-4,4,len(self.data))
            self.data['yshift'] = np.random.uniform(-2.25,2.25,len(self.data))
        else:
            self.data['xshift'] = 0
            self.data['yshift'] = 0

        self.data['x'] = self.data['x'] + self.data['xshift']
        self.data['y'] = self.data['y'] + self.data['yshift']

        '''Fiducialize in z. We'll put the z-vertex at 0 and then add a drift length to this to figure out the primary track z position'''
        self.data['z'] = self.data['z'].apply(lambda x: x-x[0])
        self.data['drift_length'] = np.random.uniform(min_drift_length,max_drift_length,len(self.data)) #cm
        self.data['z'] = self.data['drift_length']+self.data['z']

        #Now let's fiducialize by removing all coordinates with z < 0 and z > drift_gap. We will add flags for whether the track hits the GEM or the cathode
        self.data['GEM_clipped'] = self.data['z'].apply(lambda x: 1 if len(np.where(x<0)[0]) > 0 else 0)
        self.data['Cathode_clipped'] = self.data['z'].apply(lambda x: 1 if len(np.where(x>self.drift_gap_length)[0]) > 0 else 0)
        self.data['fiducial_indices'] = self.data['z'].apply(lambda x: np.where((x >= 0) & (x <= self.drift_gap_length))[0])
        #remove nonfiducial coordinates
        for col in ['x','y','z']:
            self.data[col] = self.data.apply(lambda x: x[col][x['fiducial_indices']], axis=1)
        
        '''Apply drift+diffusion if specified'''
        if drift:            
            xdiff = []
            ydiff = []
            zdiff = []
            for i in range(0,len(self.data)):
                track = self.data.iloc[i]
                xd,yd,zd = self.apply_diffusion(track['x'],track['y'],track['z'],diffusion_length=track['z'])
                xdiff.append(xd)
                ydiff.append(yd)
                zdiff.append(zd)
            self.data['xdiff'] = xdiff
            self.data['ydiff'] = ydiff
            self.data['zdiff'] = zdiff
        else:
            self.data['xdiff'] = self.data['x']
            self.data['ydiff'] = self.data['y']
            self.data['zdiff'] = self.data['z']
            
        '''Amplify if specified'''
        if amplify:
            #Populate lists of hole positions for each GEM
            if type(self.GEM_offsetsx) is list or type(self.GEM_offsetsy) is list:
                if len(self.GEM_offsetsx) != len(self.GEM_offsetsy):
                    raise ValueError("If GEM_offsetsx or GEM_offsetsy are declared as a list, both must have the same length. Modify these fields in configuration.yaml appropriately.")
                elif self.nGEM > len(self.GEM_offsetsx):
                    raise ValueError("Number of GEMs must be <= number of entries in GEM_offset{x,y}. Modify these fields in configuration.yaml accordingly.")
                else:
                    #Create a list of GEM hole positions with appropriate offset. We only fill up to nGEM, so if nGEM < len(GEM_offsets), we ignore the rest of GEM_offsets
                    self.hole_positions = [self.create_GEM_holes(offsetx=self.GEM_offsetsx[i],offsety=self.GEM_offsetsy[i]) for i in range(self.nGEM)]
                    
            else: #Make GEMs aligned if offsets are not a list
                self.hole_positions = [self.create_GEM_holes(offsetx=0,offsety=0) for i in range(self.nGEM)]

            #Create a KDTree of hole positions to quickly check if charge passes through GEM
            self.hole_trees = [KDTree(self.hole_positions[i]) for i in range(self.nGEM)]

            #Output branches
            xcam, ycam, qcam, fraccam = [], [], [], []
            xITO, zITO, qITO, fracITO = [], [], [], []

            #Apply gain and digitize on a track-by-track basis
            for i in tqdm(range(len(self.data))):
                event = self.data.iloc[i]
                #Initial track condition
                if self.migdal:
                    track = np.array([event['xdiff'],event['ydiff'],event['zdiff'],event['ID']]).T
                else:
                    track = np.array([event['xdiff'],event['ydiff'],event['zdiff']]).T

                #loop and apply amplification through each GEM
                for i in range(self.nGEM):
                    # Filter the charges that fall within the GEM holes for this stage.
                    track = np.array([charge for charge in track if self.is_within_GEMhole(charge[0], charge[1], self.hole_trees[i])])
                    if len(track) == 0:
                        # If no charge makes it through, exit the loop.
                        break
                    gap_length = transfer_gap_length + self.extra_GEM_diffusion if i != self.nGEM-1 else self.extra_GEM_diffusion
                    track = self.apply_amplification(track,gap_length)

                    
                # After processing all GEM stages, digitize the camera
                if self.migdal:
                    xc, yc, qc, fracc = self.digitize_camera_migdal(track[:, 0], track[:, 1], track[:, 3])
                    fraccam.append(fracc)
                else:
                    xc, yc, qc = self.digitize_camera(track[:, 0], track[:, 1])
                xcam.append(xc)
                ycam.append(yc)
                qcam.append(qc)

                # Apply extra induction gap diffusion before digitizing ITO
                track[:, 0] += np.sqrt(induction_gap_length) * sigmaT_induc * 1e-4 * np.random.normal(0, 1, len(track))
                track[:, 1] += np.sqrt(induction_gap_length) * sigmaT_induc * 1e-4 * np.random.normal(0, 1, len(track))
                track[:, 2] += np.sqrt(induction_gap_length) * sigmaL_induc * 1e-4 * np.random.normal(0, 1, len(track))
            
                # Digitize ITO
                if self.migdal:
                    xi, zi, qi, fraci = self.digitize_ITO_migdal(track[:, 0], track[:, 2], track[:, 3])
                    fracITO.append(fraci)
                else:
                    xi, zi, qi = self.digitize_ITO(track[:, 0], track[:, 2])
                    xITO.append(xi)
                    zITO.append(zi)
                    qITO.append(qi)
                    
                if self.gpu:
                    torch.cuda.empty_cache()

            # Save the digitized outputs to self.data.
            self.data['xcam'] = xcam
            self.data['ycam'] = ycam
            # Convert the camera charge to an unsigned 16-bit integer array.
            self.data['qcam'] = [np.array(q).astype('uint16') for q in qcam]

            self.data['xITO'] = xITO
            self.data['zITO'] = zITO
            self.data['qITO'] = qITO
            if self.migdal:
                self.data['ER_frac_cam'] = fraccam
                self.data['ER_frac_ITO'] = fracITO

        '''Save output dataframe'''
        if overwrite:
            self.data.to_feather(infile)
        else:
            if not os.path.exists(outpath):
                os.makedirs(outpath)
            outname = os.path.split(os.path.splitext(infile)[0])[1]+'_%sGEMs_%sgain_digitized.feather'%(self.nGEM,gain)
            self.data.to_feather(os.path.join(outpath,outname))

    #Method to apply GEM amplification
    def apply_amplification(self,track,gap_length):
        if self.migdal:
            # Separate out NR and ER charges.
            NRidx = np.where(track[:, 3] == 1)[0]
            ERidx = np.where(track[:, 3] == 0)[0]

            # Process the NR charges if they exist.
            if len(NRidx) > 0:
                xNR, yNR, zNR = self.GEM_gain_and_diffusion(
                    track[NRidx, 0], track[NRidx, 1], track[NRidx, 2],
                    gap_length=gap_length,
                    diff_coeff_trans=sigmaT_trans,
                    diff_coeff_long=sigmaL_trans
                )
            else:
                xNR, yNR, zNR = np.array([]), np.array([]), np.array([])

            # Process the ER charges if they exist.
            if len(ERidx) > 0:
                xER, yER, zER = self.GEM_gain_and_diffusion(
                    track[ERidx, 0], track[ERidx, 1], track[ERidx, 2],
                    gap_length=gap_length,
                    diff_coeff_trans=sigmaT_trans,
                    diff_coeff_long=sigmaL_trans
                )
            else:
                xER, yER, zER = np.array([]), np.array([]), np.array([])

            # Recombine the results and rebuild the track.
            new_x = np.concatenate((xNR, xER))
            new_y = np.concatenate((yNR, yER))
            new_z = np.concatenate((zNR, zER))
            # Recreate the ID array (order matters; here NR charges get a value of 1, ER of 0)
            new_ID = np.concatenate((np.ones(len(xNR)), np.zeros(len(xER))))
            track = np.column_stack((new_x, new_y, new_z, new_ID))
            
        else:
            # Non-migdal case: simply update x, y, z.
            xgain, ygain, zgain = self.GEM_gain_and_diffusion(
                track[:, 0], track[:, 1], track[:, 2],
                gap_length=gap_length,
                diff_coeff_trans=sigmaT_trans,
                diff_coeff_long=sigmaL_trans
            )
            track = np.column_stack((xgain, ygain, zgain))
            
            return track
    
    def apply_diffusion(self, xs, ys, zs ,diffusion_length): #applies diffusion using input diffusion parameters
        x_diff = np.sqrt(diffusion_length)*self.sigmaT*1e-4*np.random.normal(0,1, len(zs))
        y_diff = np.sqrt(diffusion_length)*self.sigmaT*1e-4*np.random.normal(0,1, len(zs))
        z_diff = np.sqrt(diffusion_length)*self.sigmaL*1e-4*np.random.normal(0,1, len(zs))
        xs = xs+x_diff
        ys = ys+y_diff
        zs = zs+z_diff
        
        del x_diff, y_diff, z_diff
        return xs,ys,zs
    '''
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
    '''
    '''Creates honeycomb grid of GEM holes with user-input diameter and pitch'''
    def create_GEM_holes(self,offsetx=0,offsety=0): #offsets of second GEM relative to first
        # Convert micrometers to centimeters
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
                x += offsetx / 10000 #convert to cm
                y += offsety / 10000
                hole_positions.append((x, y))

        shift = self.GEM_width / 2
        pos = pd.Series(hole_positions).apply(np.array)-shift
        hole_positions = pos[pos.apply(lambda x: (x[0]>-1*shift) & (x[1]>-1*shift) & (x[0] < shift) & (x[1] < shift))].apply(tuple).to_list()
        return hole_positions

    def is_within_GEMhole(self,x, y, tree):
        distance, index = tree.query([x, y])
        return distance <= self.hole_radius_cm

    def generate_gain_points(self, x, x_post, gain_electrons, gap_length, diff_coeff): #Generates x, y, and z coordiantes after gain
        for enum, val in np.ndenumerate(gain_electrons):
            start_ind = np.sum(gain_electrons[:enum[0]])
            end_ind = np.sum(gain_electrons[:enum[0]+1])
            x_post[start_ind:end_ind] = x[enum] + np.sqrt(gap_length)*diff_coeff*1E-4*np.random.normal(0,1,val)

    def generate_gain_points_GPU(self, x, gain_electrons, gap_length, diff_coeff):
        gain_electrons = gain_electrons.to('cuda')
        x = torch.tensor(x).to('cuda')
        gain_cumsum = torch.cumsum(gain_electrons, dim=0)
        start_indices = torch.zeros_like(gain_electrons, dtype=torch.int, device='cuda')
        start_indices[1:] = gain_cumsum[:-1]
        end_indices = gain_cumsum
        noise = np.sqrt(gap_length) * diff_coeff * 1E-4 * torch.randn(int(torch.sum(gain_electrons).item()), device='cuda')
        x_repeated = x.repeat_interleave(gain_electrons)
        return x_repeated+noise

    def GEM_gain_and_diffusion(self,x,y,z,gap_length, diff_coeff_trans, diff_coeff_long): #Applies gain and readout resolution smearing

        if not self.gpu:
            gain_electrons = np.random.exponential(self.gain, len(x))
            gain_electrons = np.asarray(gain_electrons, dtype=int)
            x_post = np.ascontiguousarray(np.zeros(np.sum(gain_electrons)),dtype=np.float32)
            y_post = np.ascontiguousarray(np.zeros(np.sum(gain_electrons)),dtype=np.float32)
            z_post = np.ascontiguousarray(np.zeros(np.sum(gain_electrons)),dtype=np.float32)

            self.generate_gain_points(x, x_post, gain_electrons, gap_length = gap_length, diff_coeff = diff_coeff_trans)
            self.generate_gain_points(y, y_post, gain_electrons, gap_length = gap_length, diff_coeff = diff_coeff_trans)
            self.generate_gain_points(z, z_post, gain_electrons, gap_length = gap_length, diff_coeff = diff_coeff_long)

            return x_post, y_post, z_post

        else:
            exponential_dist = torch.distributions.Exponential(1.0 / self.gain)
            gain_electrons = exponential_dist.sample((len(x),)).to(dtype=torch.int)
            x_post = self.generate_gain_points_GPU(x, gain_electrons, gap_length = gap_length, diff_coeff = diff_coeff_trans)
            y_post = self.generate_gain_points_GPU(y, gain_electrons, gap_length = gap_length, diff_coeff = diff_coeff_trans)
            z_post = self.generate_gain_points_GPU(z, gain_electrons, gap_length = gap_length, diff_coeff = diff_coeff_long)

            return x_post.cpu().numpy(), y_post.cpu().numpy(), z_post.cpu().numpy()

    def digitize_camera_migdal(self,x,y,ID):
        NRidx = np.where(ID == 1)[0]
        ERidx = np.where(ID == 0)[0]
        NRhist = np.histogram2d(x[NRidx],y[NRidx],bins=(2048,1152),range=((-4,4),(-2.25,2.25)))[0].T
        ERhist = np.histogram2d(x[ERidx],y[ERidx],bins=(2048,1152),range=((-4,4),(-2.25,2.25)))[0].T
        totalhist = NRhist + ERhist
        fraction = np.divide(ERhist, totalhist, out=np.zeros_like(ERhist, dtype=float), where=totalhist != 0)
        sparse_hist = np.where(totalhist > 0)
        y,x = sparse_hist
        q = totalhist[sparse_hist]
        frac = fraction[sparse_hist]
        return x,y,q,frac
    
    def digitize_camera(self,x,y):
        a = np.histogram2d(x,y,bins=(2048,1152),range=((-4,4),(-2.25,2.25)))[0].T
        sparse_hist = np.where(a > 0)
        y,x = sparse_hist
        q = a[sparse_hist]
        return x,y,q

    def digitize_ITO_migdal(self,x,z,ID):
        NRidx = np.where(ID == 1)[0]
        ERidx = np.where(ID == 0)[0]
        NRhist = np.histogram2d(x[NRidx],z[NRidx],bins=(120,150),range=((-5,5),(0,3.9)))[0].T
        ERhist = np.histogram2d(x[ERidx],z[ERidx],bins=(120,150),range=((-5,5),(0,3.9)))[0].T
        totalhist = NRhist + ERhist
        fraction = np.divide(ERhist, totalhist, out=np.zeros_like(ERhist, dtype=float), where=totalhist != 0)
        sparse_hist = np.where(totalhist > 0)
        z,x = sparse_hist
        q = totalhist[sparse_hist]
        frac = fraction[sparse_hist]
        return x,z,q,frac

    def digitize_ITO(self,x,z):
        a = np.histogram2d(x,z,bins=(120,150),range=((-5,5),(0,3.9)))[0].T
        sparse_hist = np.where(a > 0)
        z,x = sparse_hist
        q = a[sparse_hist]
        return x,z,q
    
if __name__ == '__main__':
    import yaml
    from os import sys
    
    '''Load configuration.yaml'''
    with open('configuration.yaml','r') as cfg:
        config = yaml.safe_load(cfg)
        degrad_cfg = config['Degrad_card']
        settings = config['Sim_settings']
        gas_cfg = config['Gas_props']
        tpc_cfg = config['TPC_sim']
    if len(sys.argv) == 1:
        infile = settings['digitization_input_file']
    else:
        infile = sys.argv[1]
    drift = settings['apply_drift']
    nGEM = tpc_cfg['nGEM']
    min_drift_length = tpc_cfg['min_drift_length']
    max_drift_length = tpc_cfg['max_drift_length']
    v_drift = gas_cfg['vd']
    sigmaT = gas_cfg['sigmaT']
    sigmaL = gas_cfg['sigmaL']
    sigmaT_trans = gas_cfg['sigmaT_trans']
    sigmaL_trans = gas_cfg['sigmaL_trans']
    sigmaT_induc = gas_cfg['sigmaT_induc']
    sigmaL_induc = gas_cfg['sigmaL_induc']
    drift_gap_length = tpc_cfg['drift_gap_length']
    GEM_width = tpc_cfg['GEM_width']
    GEM_height = tpc_cfg['GEM_height']
    GEM_thickness = tpc_cfg['GEM_thickness']
    hole_diameter = tpc_cfg['hole_diameter']
    hole_pitch = tpc_cfg['hole_pitch']
    trans_gap = tpc_cfg['transfer_gap_length']
    induc_gap = tpc_cfg['induction_gap_length']
    GEM_offsetsx = tpc_cfg['GEM_offsetsx']
    GEM_offsetsy = tpc_cfg['GEM_offsetsy']
    extra_GEM_diffusion = tpc_cfg['extra_GEM_diffusion']
    amplify = settings['apply_amplification']
    gain = tpc_cfg['gain']
    randomize_position = settings['randomize_position']
    write_gain = settings['write_gain']
    overwrite = settings['overwrite_output']
    outpath = settings['output_dir']
    gpu = settings['gpu']

    if gpu:
        import torch

    process_tracks(infile = infile,
                   outpath = outpath,
                   nGEM = nGEM,
                   drift = drift,
                   min_drift_length = min_drift_length,
                   max_drift_length = max_drift_length,
                   v_drift = v_drift,
                   sigmaT = sigmaT,
                   sigmaL = sigmaL,
                   sigmaT_trans = sigmaT_trans,
                   sigmaL_trans = sigmaL_trans,
                   sigmaT_induc = sigmaT_induc,
                   sigmaL_induc = sigmaL_induc,
                   drift_gap_length = drift_gap_length,
                   GEM_width=GEM_width,
                   GEM_height=GEM_height,
                   GEM_thickness = GEM_thickness,
                   hole_diameter=hole_diameter,
                   hole_pitch=hole_pitch,
                   transfer_gap_length = trans_gap,
                   induction_gap_length = induc_gap,
                   GEM_offsetsx = GEM_offsetsx,
                   GEM_offsetsy = GEM_offsetsy,
                   extra_GEM_diffusion = extra_GEM_diffusion,
                   randomize_position = randomize_position,
                   amplify=amplify,
                   write_gain=write_gain,
                   gain=gain,
                   overwrite=overwrite,
                   use_gpu = gpu)
