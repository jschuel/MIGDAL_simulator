'''JTS 06/26/2024: Wrapper code for running FORTRAN executables. This code adjusts the number of tracks, energy, and seed
values from the FORTRAN input card file. After running the executable, the code processes the output into a track-indexed
pandas dataframe that includes the number of ionization electrons, x,y, and z information of the ER (in um), and the timing 
information in ps. A flag is also included for  1, fluorescence; 2, pair production; 3, bremsstrahlung; 0, otherwise'''

import pandas as pd
import numpy as np
import os
from tqdm import tqdm
import time
from scipy.spatial import KDTree #Determines if charge passes through GEM holes FAST

class process_tracks:
    def __init__(self,infile,outpath,drift,min_drift_length,max_drift_length,v_drift,sigmaT,
                 sigmaL, sigmaT_trans,sigmaL_trans, sigmaT_induc,sigmaL_induc, W, GEM_width,
                 GEM_height, GEM_thickness,hole_diameter, hole_pitch, amplify,gain,
                 transfer_gap_length, induction_gap_length, GEM2_offsetx, GEM2_offsety,
                 truth_dir = np.array([1,0,0]), write_gain = False, overwrite = False):

        self.gain=gain
        self.W = W
        self.v_drift = v_drift
        self.sigmaT = sigmaT
        self.sigmaL = sigmaL
        self.GEM_width = GEM_width
        self.GEM_height = GEM_height
        self.hole_diameter = hole_diameter
        self.hole_pitch = hole_pitch
        self.GEM2_offsetx = GEM2_offsetx
        self.GEM2_offsety = GEM2_offsety
        self.truth_dir = truth_dir
        print(infile)

        self.data = pd.read_feather(infile)
        if 'ID' in self.data.columns: #autocheck if the event is a Migdal
            migdal = True
        else:
            migdal = False

        '''Apply drift if specified
        TODO: make drift length ranges a user input parameter in configuration.yaml'''
        if drift:
            self.data['drift_length'] = np.random.uniform(min_drift_length,max_drift_length,len(self.data)) #cm
            self.data['z'] = self.data['z'].apply(lambda x: x-x.mean())
            self.data['z'] = self.data['drift_length']+self.data['z']
            '''Apply diffusion'''
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
            
        '''Amplify if specified'''
        if amplify:
            #GEM 1 with KDtree
            self.hole_positions1 = self.create_GEM_holes() #no offset in GEM holes
            self.hole_tree1 = KDTree(self.hole_positions1)
            #GEM 2 with KDtree
            self.hole_positions2 = self.create_GEM_holes(offsetx=self.GEM2_offsetx,offsety=self.GEM2_offsety)
            self.hole_tree2 = KDTree(self.hole_positions2)
            xamps = []
            yamps = []
            zamps = []
            IDamps = [] #only relevant for Migdals...this keeps track of ER vs NR
            xgains = []
            ygains = []
            zgains = []
            IDgains = [] #only relevant for Migdals...this keeps track of ER vs NR
            xgains2 = [] #Position after second GEM gain
            ygains2 = []
            zgains2 = []
            IDgains2 = [] #only relevant for Migdals...this keeps track of ER vs NR
            xcam = [] #digitized camera
            ycam = []
            qcam = []
            fraccam = [] #only relevant for Migdals...this gives fractional ER content in pixel
            xITO = [] #digitized ITO
            zITO = []
            qITO = []
            fracITO = []
            for i in tqdm(range(len(self.data))):
                event = self.data.iloc[i]
                if migdal:
                    track = np.array([event['xdiff'],event['ydiff'],event['zdiff'],event['ID']]).T
                else:
                    track = np.array([event['xdiff'],event['ydiff'],event['zdiff']]).T
                charges_passing_through = np.array([charge for charge in track if self.is_within_GEMhole(charge[0], charge[1], self.hole_tree1)])
                xamp = charges_passing_through[:,0]
                yamp = charges_passing_through[:,1]
                zamp = charges_passing_through[:,2]
                if migdal:
                    IDamp = charges_passing_through[:,3]
                    IDamps.append(IDamp)
                xamps.append(xamp)
                yamps.append(yamp)
                zamps.append(zamp)
                '''Apply first amplification'''
                if migdal:
                    NRidx = np.where(IDamp == 1)[0]
                    ERidx = np.where(IDamp == 0)[0]
                    xgNR, ygNR, zgNR = self.GEM_gain_and_diffusion(xamp[NRidx],yamp[NRidx],zamp[NRidx],gap_length = transfer_gap_length + GEM_thickness / 3, diff_coeff_trans = sigmaT_trans, diff_coeff_long = sigmaL_trans)
                    xgER, ygER, zgER = self.GEM_gain_and_diffusion(xamp[ERidx],yamp[ERidx],zamp[ERidx],gap_length = transfer_gap_length + GEM_thickness / 3, diff_coeff_trans = sigmaT_trans, diff_coeff_long = sigmaL_trans)
                    IDgainER = [0 for i in range(0,len(xgER))]
                    IDgainNR = [1 for i in range(0,len(xgNR))]
                    xgain = np.concatenate((xgNR,xgER))
                    ygain = np.concatenate((ygNR,ygER))
                    zgain = np.concatenate((zgNR,zgER))
                    IDgain = np.concatenate((IDgainNR,IDgainER))
                else:
                    xgain,ygain,zgain = self.GEM_gain_and_diffusion(xamp,yamp,zamp,gap_length = transfer_gap_length + GEM_thickness / 3, diff_coeff_trans = sigmaT_trans, diff_coeff_long = sigmaL_trans) #transfer gap is 0.2cm GEM_thickness is in cm
                xgains.append(xgain)
                ygains.append(ygain)
                zgains.append(zgain)
                if migdal:
                    IDgains.append(IDgain)
                    
                '''Figure out which charges make it thru second GEM hole'''
                if migdal:
                    gaintrack = np.array([xgain,ygain,zgain,IDgain]).T
                else:
                    gaintrack = np.array([xgain,ygain,zgain]).T
                charges_passing_through_gain = np.array([charge for charge in gaintrack if self.is_within_GEMhole(charge[0], charge[1], self.hole_tree2)])
                xgainamp = charges_passing_through_gain[:,0]
                ygainamp = charges_passing_through_gain[:,1]
                zgainamp = charges_passing_through_gain[:,2]
                if migdal:
                    IDgainamp = charges_passing_through_gain[:,3]
                '''Apply second amplification for camera'''
                if migdal:
                    NRgainidx = np.where(IDgainamp == 1)[0]
                    ERgainidx = np.where(IDgainamp == 0)[0]
                    xgNR2, ygNR2, zgNR2 = self.GEM_gain_and_diffusion(xgainamp[NRgainidx],ygainamp[NRgainidx],zgainamp[NRgainidx],gap_length = 0 + GEM_thickness / 3, diff_coeff_trans = sigmaT_trans, diff_coeff_long = sigmaL_trans)
                    xgER2, ygER2, zgER2 = self.GEM_gain_and_diffusion(xgainamp[ERgainidx],ygainamp[ERgainidx],zgainamp[ERgainidx],gap_length = 0 + GEM_thickness / 3, diff_coeff_trans = sigmaT_trans, diff_coeff_long = sigmaL_trans)
                    IDgainER2 = [0 for i in range(0,len(xgER2))]
                    IDgainNR2 = [1 for i in range(0,len(xgNR2))]
                    xgain2 = np.concatenate((xgNR2,xgER2))
                    ygain2 = np.concatenate((ygNR2,ygER2))
                    zgain2 = np.concatenate((zgNR2,zgER2))
                    IDgain2 = np.concatenate((IDgainNR2,IDgainER2))
                else:
                    xgain2,ygain2,zgain2 = self.GEM_gain_and_diffusion(xgainamp,ygainamp,zgainamp,gap_length = 0 + GEM_thickness / 3, diff_coeff_trans = sigmaT_trans, diff_coeff_long = sigmaL_trans) #transfer gap is 0.2cm; GEM_thickness is in cm
                xgains2.append(xgain2)
                ygains2.append(ygain2)
                zgains2.append(zgain2)
                if migdal:
                    IDgains2.append(IDgain2)
                '''Digitize camera readout'''
                if migdal:
                    xc,yc,qc,fracc = self.digitize_camera_migdal(xgain2,ygain2,IDgain2)
                else:
                    xc,yc,qc = self.digitize_camera(xgain2,ygain2)
                xcam.append(xc)
                ycam.append(yc)
                qcam.append(qc)
                if migdal:
                    fraccam.append(fracc)
            
            self.data['xdiff'] = xdiff
            self.data['ydiff'] = ydiff
            self.data['zdiff'] = zdiff
                
            self.data['xamp'] = xamps
            self.data['yamp'] = yamps
            self.data['zamp'] = zamps
            if migdal:
                self.data['IDamp'] = IDamps

            self.data['xgain'] = xgains
            self.data['ygain'] = ygains
            self.data['zgain'] = zgains
            if migdal:
                self.data['IDgain'] = IDgains

            self.data['xgain2'] = xgains2
            self.data['ygain2'] = ygains2
            self.data['zgain2'] = zgains2
            if migdal:
                self.data['IDgain2'] = IDgains2

            self.data['xcam'] = xcam
            self.data['ycam'] = ycam
            self.data['qcam'] = qcam
            self.data['qcam'] = self.data['qcam'].apply(lambda x: x.astype('int16'))
            if migdal:
                self.data['ER_frac_cam'] = fraccam

            '''Apply extra induction gap diffusion for ITO sim'''
            self.data['xgain2'] = self.data['xgain2'].apply(lambda x: x+np.sqrt(induction_gap_length)*sigmaT_induc*1e-4*np.random.normal(0,1,len(x)))
            self.data['ygain2'] = self.data['ygain2'].apply(lambda x: x+np.sqrt(induction_gap_length)*sigmaT_induc*1e-4*np.random.normal(0,1,len(x)))
            self.data['zgain2'] = self.data['zgain2'].apply(lambda x: x+np.sqrt(induction_gap_length)*sigmaL_induc*1e-4*np.random.normal(0,1,len(x)))

            '''Digitize ITO'''
            for i in range(0,len(self.data)):
                track = self.data.iloc[i]
                if migdal:
                    xi,zi,qi,fraci = self.digitize_ITO_migdal(track['xgain2'],track['zgain2'],track['IDgain2'])
                    fracITO.append(fraci)
                else:
                    xi,zi,qi = self.digitize_ITO(track['xgain2'],track['zgain2'])
                xITO.append(xi)
                zITO.append(zi)
                qITO.append(qi)

            self.data['xITO'] = xITO
            self.data['zITO'] = zITO
            self.data['qITO'] = qITO
            if migdal:
                self.data['ER_frac_ITO'] = fracITO
            
            if not write_gain:
                if migdal:
                    del(self.data['xgain'],self.data['xgain2'],self.data['ygain'],
                        self.data['ygain2'], self.data['zgain'], self.data['zgain2'],
                        self.data['IDgain']), self.data['IDgain2'])
                else:
                    del(self.data['xgain'],self.data['xgain2'],self.data['ygain'],
                    self.data['ygain2'], self.data['zgain'], self.data['zgain2'])
                
        '''Save output dataframe'''
        if overwrite:
            self.data.to_feather(infile)
        else:
            if not os.path.exists(outpath):
                os.makedirs(outpath)
            outname = os.path.split(os.path.splitext(infile)[0])[1]+'_digitized.feather'
            self.data.to_feather(outpath+'/'+outname)
    
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

    def GEM_gain_and_diffusion(self,x,y,z,gap_length, diff_coeff_trans, diff_coeff_long): #Applies gain and readout resolution smearing
        gain_electrons = np.random.exponential(np.sqrt(self.gain), len(x))
        gain_electrons = np.asarray(gain_electrons, dtype=int)
        
        x_post = np.ascontiguousarray(np.zeros(np.sum(gain_electrons)),dtype=np.float32)
        y_post = np.ascontiguousarray(np.zeros(np.sum(gain_electrons)),dtype=np.float32)
        z_post = np.ascontiguousarray(np.zeros(np.sum(gain_electrons)),dtype=np.float32)

        self.generate_gain_points(x, x_post, gain_electrons, gap_length = gap_length, diff_coeff = diff_coeff_trans)
        self.generate_gain_points(y, y_post, gain_electrons, gap_length = gap_length, diff_coeff = diff_coeff_trans)
        self.generate_gain_points(z, z_post, gain_electrons, gap_length = gap_length, diff_coeff = diff_coeff_long)

        return x_post, y_post, z_post

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
        
    infile = settings['digitization_input_file']
    W = gas_cfg['W']
    drift = settings['apply_drift']
    min_drift_length = tpc_cfg['min_drift_length']
    max_drift_length = tpc_cfg['max_drift_length']
    v_drift = gas_cfg['vd']
    sigmaT = gas_cfg['sigmaT']
    sigmaL = gas_cfg['sigmaL']
    sigmaT_trans = gas_cfg['sigmaT_trans']
    sigmaL_trans = gas_cfg['sigmaL_trans']
    sigmaT_induc = gas_cfg['sigmaT_induc']
    sigmaL_induc = gas_cfg['sigmaL_induc']
    GEM_width = tpc_cfg['GEM_width']
    GEM_height = tpc_cfg['GEM_height']
    GEM_thickness = tpc_cfg['GEM_thickness']
    hole_diameter = tpc_cfg['hole_diameter']
    hole_pitch = tpc_cfg['hole_pitch']
    trans_gap = tpc_cfg['transfer_gap']
    induc_gap = tpc_cfg['induction_gap']
    GEM2_offsetx = tpc_cfg['GEM2_offsetx']
    GEM2_offsety = tpc_cfg['GEM2_offsety']
    amplify = settings['apply_amplification']
    gain = tpc_cfg['gain']
    write_gain = settings['write_gain']
    overwrite = settings['overwrite_output']
    outpath = settings['output_dir']

    process_tracks(infile = infile,
                   outpath = outpath,
                   W = W,
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
                   GEM_width=GEM_width,
                   GEM_height=GEM_height,
                   GEM_thickness = GEM_thickness,
                   hole_diameter=hole_diameter,
                   hole_pitch=hole_pitch,
                   transfer_gap_length = trans_gap,
                   induction_gap_length = induc_gap,
                   GEM2_offsetx = GEM2_offsetx,
                   GEM2_offsety = GEM2_offsety,
                   amplify=amplify,
                   write_gain=write_gain,
                   gain=gain,
                   overwrite=overwrite)
