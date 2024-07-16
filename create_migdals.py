import pandas as pd
from tqdm import tqdm
import numpy as np
import os
pd.options.mode.chained_assignment = None

class form_migdals:
    def __init__(self,ERinfile,NRinfile,outdir,randomize = False):
        self.ERinfile = ERinfile
        self.NRinfile = NRinfile
        self.outdir = outdir
        
        '''Load tracks'''
        self.ERs = self.load_ERs()
        self.NRs = self.load_NRs()

        if randomize: #randomize track ordering
            self.NRs = self.NRs.sample(frac = 1)
            self.ERs = self.ERs.sample(frac = 1)
            self.ERs.index = [i for i in range(0,len(self.ERs))]
            self.NRs.index = [i for i in range(0,len(self.NRs))]
            
        '''Form Migdals'''
        self.migdals = pd.concat([self.make_migdal(i) for i in tqdm(range(0,min(len(self.ERs),len(self.NRs))))])
        self.migdals.index = [i for i in range(0,len(self.migdals))]

        '''Save'''
        self.save()

    def load_ERs(self):
        return pd.read_feather(self.ERinfile)

    def load_NRs(self):
        return pd.read_feather(self.NRinfile)

    def make_migdal(self,i):
        mig = pd.DataFrame()
        nr = self.NRs.iloc[i]
        nr['nHits'] = len(nr['x'])
        er = self.ERs.iloc[i]
        for col in ['nHits','truthE','ionizationE','truth_dir','truth_theta','truth_phi']:
            try:
                mig['ER_%s'%(col)] = [er[col]]
                mig['NR_%s'%(col)] = [nr[col]]
            except:
                pass
        er['ID'] = [0 for i in range(0,len(er['x']))]
        nr['ID'] = [1 for i in range(0,len(nr['x']))]
        for col in ['x','y','z']:
            er['shift%s'%(col)] = er[col][0]-nr[col][0]
            er[col] = er[col] - er['shift%s'%(col)]
            mig[col] = [list(nr[col])+list(er[col])]
        mig['ID'] = [nr['ID']+er['ID']]
        for col in [col for col in mig.columns if mig[col].dtype == 'O']:
            mig[col] = np.array(mig[col])
        return mig

    def save(self):
        outfile = 'ER_%s-%skeV_%s_NRs_%s_tracks_migdals.feather'%(self.ERs['truthE'].min(),self.ERs['truthE'].max(),os.path.split(self.NRinfile)[1][0],len(self.migdals))
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
        self.migdals.to_feather(self.outdir+outfile)

if __name__ == '__main__':
    import yaml

    '''Load configuration.yaml'''
    with open('configuration.yaml','r') as cfg:
        config = yaml.safe_load(cfg)
        mig_cfg = config['Migdal_sim']

    form_migdals(ERinfile = mig_cfg['ER_input_filepath'],
                     NRinfile = mig_cfg['NR_input_filepath'],
                     outdir = mig_cfg['output_directory'],
                     randomize = mig_cfg['randomize'])
