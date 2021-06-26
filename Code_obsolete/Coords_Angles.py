import Bio.PDB
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.PDB.Polypeptide import is_aa
from Calc_angle import Calc_angles
import os
import csv

warnings.filterwarnings("ignore", category=PDBConstructionWarning)
backbone_ids = ['C', 'N', 'CA']
current_directory = os.path.dirname(os.path.abspath(__file__))
PDB_path = os.path.join(current_directory, 'PDBfiles/')

class Get_PsiandPhi:

    def __init__(self, PDB, path = PDB_path, CSV_out = 'CSVfiles/' ,auto = True):
        '''
        :param PDB: ficheiro do PDB
        :param path: path do ficheiro
        :param auto: Processo automatico
        '''
        self.path = path
        self.file = self.path + PDB + '.cif'
        self.name = PDB
        self.structure = Bio.PDB.MMCIFParser().get_structure(self.name,self.file)
        self.out= CSV_out
        if auto == True:
            self.get_coords()
            self.angulos_todos()

    def get_seq_chain(self):
        res = []
        parser = Bio.PDB.MMCIFParser()
        struct = parser.get_structure(self.name,self.file)
        for model in struct:
            for chain in model:
                seq = ''
                for residue in chain:
                    if is_aa(residue, standard=True):
                        seq += residue.resname+','
                res.append([chain.id, seq[:len(seq)-1]])
        return res

    def get_coords(self):
        backbone_at = []
        backbone_coord = []
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    atoma = []
                    atm_coord= []
                    if is_aa(residue, standard=True):  # Check if amino acid
                        for atom in residue:
                            if atom.get_name() in backbone_ids:
                                atoma.append(atom.get_name())
                                atm_coord.append(atom.get_coord().tolist())
                        backbone_coord = backbone_coord + [atm_coord]
                        backbone_at = backbone_at + [[residue.resname, atoma]]
        self.backbone_coord, self.backbone_at = backbone_coord, backbone_at

    def all_phi(self): #considera faulty files
        phi = []
        for i in range(1,len(self.backbone_at)):
            if len(self.backbone_coord[i]) ==3 and 'C' in self.backbone_at[i-1][1]: #verifica se tem os 3 atomos para o residuo e se existe um C no residuo anterior
                index = self.backbone_at[i - 1][1].index('C') #indexa o C do residuo anterior caso não tenha os atomos todos
                C1= self.backbone_coord[i-1][index]
                N= self.backbone_coord[i][0]
                CA= self.backbone_coord[i][1]
                C2= self.backbone_coord[i][2]
                Coords_phi = [C1,N,CA,C2]
                ang = Calc_angles(Coords_phi)
                a_phi = ang.DihedralAngle()
                a_phi = round(a_phi, 5)
                phi.append(a_phi)
            else:
                phi.append('NA')
        return phi

    def all_psi(self):
        psi=[]
        for i in range(len(self.backbone_at)-1):
            if len(self.backbone_coord[i]) == 3 and 'N' in self.backbone_at[i+1][1]:  #verifica se tem os 3 atomos para o residuo e se existe um N no residuo anterior
                index = self.backbone_at[i+1][1].index('N') #indexa o N do residuo anterior caso não tenha os atomos todos
                N1= self.backbone_coord[i][0]
                CA= self.backbone_coord[i][1]
                C= self.backbone_coord[i][2]
                N2 = self.backbone_coord[i+1][index]
                Coords_psi = [N1,C,CA,N2]
                ang = Calc_angles(Coords_psi)
                a_psi = ang.DihedralAngle()
                psi.append(-a_psi)
            else:
                psi.append('NA')
        return psi

    def angulos_todos(self):
        psi, phi = self.all_psi(), self.all_phi()
        matrix_res = []
        matrix_res.append([psi[0],'NA'])
        for i in range(1, len(self.backbone_at)-1):
            matrix_res.append([psi[i], phi[i-1]])
        last = len(self.backbone_at)-1
        matrix_res.append(['NA',phi[last-1]])
        return matrix_res

    def angulos_todos_pseq(self):
        psi, phi = self.all_psi(), self.all_phi()
        matrix_res = []
        matrix_res.append([self.backbone_at[0][0],[psi[0],'NA']])
        for i in range(1, len(self.backbone_at)-1):
            matrix_res.append([self.backbone_at[i][0], [psi[i], phi[i-1]]])
        last = len(self.backbone_at)-1
        matrix_res.append([self.backbone_at[last][0], ['NA',phi[last-1]]])
        self.csv_out(matrix_res)
        return matrix_res

    def csv_out(self, lista):
        '''
        :return: cria ficheiro .csv com os outputs da conversão
        '''
        file = open('CSVfiles/' + self.name + '.csv', 'w', newline='')
        header = ['AA', "PSI", 'PHI']
        writer = csv.DictWriter(file, fieldnames=header)
        writer.writeheader()
        for i in lista:
            AA, angulos = i[0], i[1]
            PSI, PHI = angulos
            if PSI != None and PHI != None:
                writer.writerow({'AA': AA,
                                 "PSI": PSI,
                                 'PHI': PHI})

