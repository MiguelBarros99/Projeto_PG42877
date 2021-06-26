import pandas as pd
import urllib.parse
import urllib.request

class Convert:
    def __init__(self, file,coluna, orig,fin, auto = True):
        '''
        :param file: ficheio .csv para converter
        :param coluna: Coluna onde está o ID Uniprot
        :param out: Nome para ficheiro gerado
        :param auto: fazer o processo inteiro devolvendo um dicionario
        '''
        self.csv = file
        self.col = coluna
        self.orig = orig
        self.fin = fin
        if auto ==True:
            self.organize_dic()

    def isolar_IDs(self):
        '''
        :return: Lista ALL com todos os ID's Uniprot no self.csv
        '''
        df = pd.read_csv(self.csv)
        IDs= df[self.col] #isolar IDS do dataset
        ALL = []
        for i in IDs: #isolar os IDs para a lista ALL
            ALL.append(i)
        return ALL

    def create_query(self):
        '''
        :return: Cria a query(string) para a conversão de ID's
        '''
        lista = self.isolar_IDs()
        str = ' '.join(lista)
        return str

    def converter(self):
        '''
        :return: convertor de IS's da Uniprot para ID's do PDB
        '''
        query = self.create_query()
        url = 'https://www.uniprot.org/uploadlists/'
        res = []
        params = {'from': self.orig, 'to': self.fin, 'format': 'tab', 'sort': 'score', 'query': query}
        data = urllib.parse.urlencode(params)
        data = data.encode('utf-8')
        req = urllib.request.Request(url, data)
        with urllib.request.urlopen(req) as f:
            response = f.read()
        x = response.decode('utf-8').splitlines()
        for i in x[1:]:
            uni, pdb = i.split('\t')
            res.append((uni, pdb))
        return res

    def organize_dic(self):
        '''
        :return: converte o output para um dicionario
        '''
        lista = self.converter()
        d={}
        for i in lista:
            k,l = i
            if k not in d:
                d[k] = [l]
            else:
                d[k] = d[k]+ [l]
        return d


