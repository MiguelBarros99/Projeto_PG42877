import pandas as pd
from oneD_hist import Histogram
from Auxiliary import csv_out
from Bio_psi_phi import Psi_phi
import os

##### Conjugar datasets #######

print('#### juntar datasets ####')
PDBreps = []
Unireps = []

f = open('textfiles/AMIDI.txt', "r")
x = f.readlines()
dic = {}
for i in x:
    pdb=[]
    xi = i.split()
    uni, ids = xi[1], len(xi)-3
    Unireps.append(uni)
    for i in range(ids):
        pdb.append(xi[i+3])
        if xi[i+3] not in PDBreps:
            PDBreps.append(xi[i+3])
    dic[uni]= ' '.join(pdb)

k = open('textfiles/ECPRED.txt', "r")
kx = k.readlines()
for i in kx:
    pdb = []
    xi = i.split()
    uni, ids = xi[1], len(xi) - 3
    if uni not in Unireps: Unireps.append(uni)
    for i in range(ids):
        pdb.append(xi[i + 3])
        if xi[i+3] not in PDBreps:
            PDBreps.append(xi[i+3])
    if uni not in dic:
        dic[uni] = ' '.join(pdb)
    else:
        xi = dic[uni].split()
        for i in pdb:
            if i not in xi:
                xi.append(i)
        dic[uni] = ' '.join(xi)

print(' UNIPROTS: '+ str(len(Unireps)) + ' PDBs: ' + str(len(PDBreps)))

DIC_S = {}
PDBf = []
for i in dic.keys():
    xi = dic[i].split(' ')
    ids =[]
    for k in xi:
        if os.path.isfile('PDBfiles/' + k + '.cif'):
            ids.append(k)
            if k not in PDBf: PDBf.append(k)
    if ids != []: DIC_S[i] = ' '.join(ids)

print(' UNIPROTS: '+ str(len(DIC_S.keys())) + ' PDBs: ' + str(len(PDBf)))

f = open('textfiles/fulldata.txt', "w+")
for i in DIC_S.keys():
    f.write('_uni ' +i+' _pdb ' + DIC_S[i] + '\n')
f.close()

####### gerar dicionario final e csv com a feature ############
print('#### Feature extraction e transformação do dataset ####')

k = open('textfiles/fulldata.txt', "r")
kx = k.readlines()
dic = {}
for i in kx:
    xi = i.split()
    pdb=[]
    uni, ids = xi[1], len(xi) - 3
    for i in range(ids):
        pdb.append(xi[i + 3])
    for i in pdb:
        if i not in dic:
            dic[i]= uni
        else:
            xi = dic[i].split(' ')
            if uni not in xi:
                xi.append(uni)
                dic[i] = ' '.join(xi)

UNI = []
DIC_F ={}
total, cont = 0,0
for i in dic.keys():
    total += 1
    if os.path.isfile('PDBfiles/' + i + '.cif'):
        cont +=1
        # print(i+ ' in progress...')
        DIC_F[i] = dic[i]
        xi = dic[i].split(' ')
        for x in xi:
            if x not in UNI:
                UNI.append(x)
            if os.path.isfile('BioCSV/' + i + '.csv') == False:
                res2 =Psi_phi(i, path='PDBfiles/', CSV_out= 'BioCSV/')
                res2.Bio_psi_phi()

print('PDB finais: '+ str(total) + '  UniProt finais: ' + str(len(UNI)))

f = open('textfiles/finalres.txt', "w+")
for i in DIC_F.keys():
    f.write('_uni ' +DIC_F[i]+' _pdb ' + i + '\n')
f.close()

# ##### GERAR DATASET COMPLETO #######
print('### Gerar o dataset ###')

#dicionario dos ECS -> Uniprots
ECS = {}
df = pd.read_csv('textfiles/ecpred_uniprot_uniref_90.csv')
IDs_start= df['Entry']
ECs = df['ec_number']
for i in range(len(ECs)):
    x, y = IDs_start[i], ECs[i]
    y = y.replace(' ','').split(';')
    if y != 'nan':
        if x not in ECS:
            ECS[x] = y
        else:
            for EC in y:
                if EC not in ECS[x]:
                    ECS[x] = ECS[x] + [EC]

df = pd.read_csv('ECnumber_39251.csv')
IDs_39251= df['Entry']
ECs = df['EC number']
for i in range(len(ECs)):
    x, y = IDs_39251[i], str(ECs[i])
    y = y.replace(' ','').split(';')
    if y[0] != 'nan' and y[0] != ['nan']:
        if x not in ECS:
            ECS[x] = y
        else:
            for EC in y:
                if EC not in ECS[x]:
                    ECS[x] = ECS[x]+ [EC]

#PDB -> UNIPROT
f = open('textfiles/finalres.txt', "r")
x = f.readlines()
FILE = {}
contar = []
for i in x:
    xi = i.split()
    unis, pdb = xi[1:len(xi)-2],xi[len(xi)-1].strip()
    FILE[pdb]= unis
    for k in unis:
        if k not in contar: contar.append(k)

print('Uniprot inicial ' + str(len(contar)))

def att_uni(lista, add):
    for x in add:
        if x not in lista:
            lista.append(x)
    return lista

#PDB -> EC
PDB_EC = {}
PDB_UNI = {}
all_uni = []
lost=[]
for i in FILE.keys():
    ec = []
    realuni = []
    for uni in FILE[i]:
        if uni in ECS:
            realuni.append(uni)
            ec = ec + ECS[uni]
        else:
            lost = att_uni(lost, [uni])
    if ec != []:
        PDB_EC[i] = ';'.join(ec)
        PDB_UNI[i] = realuni
        all_uni= att_uni(all_uni, realuni)

print('Uniprots: perdidos - ' + str(len(lost)), ' mantidos- ' + str(len(all_uni)))
print('PDB: perdidos - ' + str(len(FILE.keys())-len(PDB_UNI)), ' mantidos- ' + str(len(PDB_UNI)))

def dif_uni(lista, lista2):
    Miss = 0
    for uni in lista:
        if uni not in lista2: Miss += 1
    return Miss

def dataset(dic_s, PDB_EC):
    totalv,totalc,totald,contagem,comec = 0,0,0,0,0
    Univ,Unic,Unicome, Unidel, Unict = [],[],[],[], []
    Fulldataset = []
    for i in dic_s.keys():
        contagem +=1
        Unict = att_uni(Unict, dic_s[i])
        if i in PDB_EC:
            comec +=1
            Unicome = att_uni(Unicome, dic_s[i])
            # print(i + ' in progress...')
            res = Histogram(i, dic_s[i])
            vetor = res.li_dataset()
            if sum(vetor[2]) != 0: #remoção dos NA's
                totalv += 1
                Univ = att_uni(Univ, dic_s[i])
                if PDB_EC[i][0] != '0':
                    Unic = att_uni(Unic, dic_s[i])
                    totalc +=1
                    vetor.append(PDB_EC[i])
                    Fulldataset.append(vetor)
                    if len(PDB_EC[i].split(';')) > 1:
                        x = PDB_EC[i].split(';')
                        last = x[0][0]
                        for j in range(1,len(x)):
                            y = x[j].strip()
                            if y[0] == last:
                                last = y[0]
                            else:
                                totald +=1
                                Unidel = att_uni(Unidel, dic_s[i])
                                del Fulldataset[len(Fulldataset)-1]
    all_uni=[]
    for k in Fulldataset:
        x = k[0].split(' ')
        for y in x:
            if y not in all_uni: all_uni.append(y)
    for i in Unicome:
        if i not in Unict: print(i)
    dele = dif_uni(Unidel, all_uni)
    nozero = dif_uni(Unic, all_uni) -dele
    gvet = dif_uni(Univ, all_uni) -dele -nozero
    semEC = dif_uni(Unicome, all_uni) -dele -nozero - gvet
    print('Sem EC: PDB- ' + str(contagem - comec) + ' UNIPROT- ' + str(semEC))
    print('Falta de dados : PDB- ' + str(comec- totalv) + ' UNIPROT- ' + str(gvet))
    print('EC = 0 : PDB- ' + str(totalv - totalc) + ' UNIPROT- ' + str(nozero))
    print('Dismatch de ECs: PBD-' + str(totald) + ' UNIPROT- ' + str(dele))
    print('Total: PDB- '+str(len(Fulldataset)) + ' Uniprot: '+ str(len(all_uni)))
    return Fulldataset

FDT = dataset(PDB_UNI, PDB_EC)
csv_out(FDT, 'Fulldata')