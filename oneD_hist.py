import numpy as np
import pandas as pd
import os
from scipy import signal

current_directory = os.path.dirname(os.path.abspath(__file__))
csv_path = os.path.join(current_directory, 'BioCSV/')

class Histogram:
    def __init__(self, PDB, UNI, path = csv_path, bin= [-180,-170,-150,-130,-110,-90,-70,-50,-30,-10,10,30,50,70,90,110,130,150,170,180]):
        '''
        :param file: ficheio .csv para converter
        '''
        self.PDB = PDB
        self.csv = path + PDB +'.csv'
        self.Uni = ' '.join(UNI)
        self.path = path
        self.bin = bin
        self.get_angle()
        self.kernel= [[1/9,1/9,1/9],[1/9,1/9,1/9],[1/9,1/9,1/9]]

    def get_angle(self):
        '''
        :return: Lista ALL com todos os ID's Uniprot no self.csv
        '''
        df = pd.read_csv(self.csv)
        df.dropna()
        self.psi = df['PSI']
        self.phi = df['PHI']

    def build_hist(self):
        ops=np.histogram2d(self.psi,self.phi, bins= self.bin, range=[[-180, +180], [-180, +180]])
        hist = ops[0]
        self.size_bins= ops
        kernel= signal.convolve2d(hist,self.kernel,'same')
        if np.sum(kernel) != 0:
            self.kernel = kernel/sum(sum(kernel))
        else:
            self.kernel = kernel

    def create_vector(self):
        self.build_hist()
        hist_v = self.kernel.flatten()
        return list(hist_v)


    def li_dataset(self):
        Info = [self.Uni,self.PDB, self.create_vector()]
        return Info


