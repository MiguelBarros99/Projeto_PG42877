from Bio.PDB import MMCIFParser
import os
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import warnings
import random

current_directory = os.path.dirname(os.path.abspath(__file__))
PDB_path = os.path.join(current_directory, 'PDBfiles/')
warnings.filterwarnings("ignore", category=PDBConstructionWarning)

class PDB_batch:
    def __init__(self, dic, out):
        '''
        :param dic: dicionario com os dados da conversão
        :param out: Output para ficheiro .txt para batch download
        '''
        self.dic = dic
        self.out = out

    def test_txt(self, out):
        '''
        :return: Imprime os resultados num .txt
        '''
        dic = self.dic
        f = open(out+'.txt', "w+")
        for i in dic.keys():
            f.write('_uni '+i+' _pdb ')
            for j in dic[i]:
                f.write(j + ' ')
            f.write('\n')
        f.close()

    def PDB_downloadbatch(self):
        All = []
        for i in self.dic.keys():
            for j in self.dic[i]:
                if j not in All:
                    All.append(j)
        All = ','.join(All)
        f = open(self.out + '.txt', "w+")
        f.write(All[:len(All)-1])
        f.close()

class PDB_info:
    def __init__(self, lista, uni, path = PDB_path, auto = True):
        '''
        :param dic: dicionario com os dados da conversão
        :param path: path para a pasta onde estão os ficheiros .cif
        :param auto: realizar o processo completo
        '''
        self.lis = lista
        self.uni = uni
        self.path = path
        if auto == True:
            self.lis_optimo = self.PDB_selection()

    def PDB_selection(self):
        res = float(99999)
        pdb_opt = 'ERRO'
        parser = MMCIFParser()
        if len(self.lis) == 1:
            if os.path.isfile(self.path+self.lis[0]+'.cif'):
                pdb_opt = self.lis[0]
                print(self.lis[0] + ' para ' + self.uni)
        else:
            ids = []
            for j in self.lis:
                if os.path.isfile(self.path+j+'.cif'):
                    ids.append(j)
                    structure = parser.get_structure(j, self.path +j+'.cif')
                    resolution = structure.header["resolution"]
                    if resolution != None and resolution < res:
                        res = resolution
                        pdb_opt = j
                else:
                    pass
            if pdb_opt == 'ERRO' and len(ids)>0:
                pdb_opt  = random.choices(ids, k=1)[0]
            if pdb_opt == 'ERRO' and len(ids) == 0: print('No files!')
            print(pdb_opt  + ' para ' + self.uni)
        return pdb_opt.upper()
    

