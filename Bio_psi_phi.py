from Bio.PDB import MMCIFParser
import os
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import warnings
import Bio.PDB
import csv
import math

warnings.filterwarnings("ignore", category=PDBConstructionWarning)
current_directory = os.path.dirname(os.path.abspath(__file__))
PDB_path = os.path.join(current_directory, 'BioCSV/')
out ='BioCSV/'

class Psi_phi:
    def __init__(self,PDB,path = PDB_path, CSV_out = out):
        self.PDB = PDB #nome do dic
        self.path = path
        self.filepath = self.path + PDB + '.cif'
        self.out = CSV_out

    def Bio_psi_phi(self):
        chains = []
        for model in Bio.PDB.MMCIFParser().get_structure(self.PDB, self.filepath):
            for chain in model:
                seq = []
                poly = Bio.PDB.Polypeptide.Polypeptide(chain)
                ta = poly.get_phi_psi_list()
                sequence = poly.get_sequence()
                for i in range(len(ta)):
                    seq.append((chain.id, sequence[i]))
                chains.append([seq,ta])
        self.csv_out(chains)


    def csv_out(self, chains):
        '''
        :return: cria ficheiro .csv com os outputs da conversão
        '''
        file = open(self.out + self.PDB + '.csv', 'w', newline='')
        header = ['Chain','AA', "PSI", 'PHI']
        writer = csv.DictWriter(file, fieldnames=header)
        writer.writeheader()
        if len(chains) != 1:
            for i in chains:
                    seq,angulos = i[0], i[1]
                    for index in range(len(angulos)):
                        chain, AA = seq[index]
                        PHI, PSI  = angulos[index]
                        if PSI != None: PSI = math.degrees(float(PSI))
                        if PHI != None: PHI = math.degrees(float(PHI))
                        if PSI != None and PHI != None:
                            writer.writerow({'Chain': chain,'AA': AA,"PSI": PSI,'PHI': PHI})
        else:
            seq, angulos = chains[0][0], chains[0][1]
            for index in range(len(angulos)):
                chain, AA = seq[index]
                PHI,PSI  = angulos[index]
                if PSI != None: PSI = math.degrees(float(PSI))
                if PHI != None: PHI = math.degrees(float(PHI))
                if PSI != None and PHI != None:
                    writer.writerow({'Chain': chain,
                                     'AA': AA,
                                     "PSI": PSI,
                                     'PHI': PHI})
