import pandas as pd

ECS = {}
df = pd.read_csv('ecpred_uniprot_uniref_90.csv')
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

print('#### FULLDATA ####')

df = pd.read_csv('Fulldata.csv')
df.dropna()
PDB = df['PDB']
UNI = df['UniProt']
Numb = df['EC']
IsoUNI = []
for i in UNI:
    xi = i.split(' ')
    for k in xi:
        if k not in IsoUNI: IsoUNI.append(k)

print('PDB finais: '+ str(len(PDB)) + '   UniProt finais: '+ str(len(IsoUNI)))

contagemUNI = {}
for j in IsoUNI:
    if j in ECS:
        Ec = ECS[j][0][0]
        if Ec not in contagemUNI:
            contagemUNI[Ec] = 1
        else: contagemUNI[Ec] += 1
print('UniProt')
print(contagemUNI )

contagemPDB = {}
for i in range(len(PDB)):
    Ec = Numb[i][0]
    if Ec not in contagemPDB:
        contagemPDB[Ec] = 1
    else:
        contagemPDB[Ec] += 1
print('PDBs')
print(contagemPDB)
#
print('#### Reduddata ####')

df = pd.read_csv('Reduddata.csv')
df.dropna()
PDB = df['PDB']
UNI = df['UniProt']
print(len(UNI))
Numb = df['EC']
IsoUNI = []
for i in UNI:
    xi = i.split(' ')
    for k in xi:
        if k not in IsoUNI: IsoUNI.append(k)

print('PDB finais: '+ str(len(PDB)) + '   UniProt finais: '+ str(len(IsoUNI)))

contagemUNI = {}
for j in IsoUNI:
    if j in ECS:
        Ec = ECS[j][0][0]
        if Ec not in contagemUNI:
            contagemUNI[Ec] = 1
        else: contagemUNI[Ec] += 1
print('UniProt')
print(contagemUNI)

contagemPDB = {}
for i in range(len(PDB)):
    Ec = Numb[i][0]
    if Ec not in contagemPDB:
        contagemPDB[Ec] = 1
    else:
        contagemPDB[Ec] += 1
print('PDBs')
print(contagemPDB)

print('#### Red_one ####')

df = pd.read_csv('Red_one_data.csv')
df.dropna()
PDB = df['PDB']
UNI = df['UniProt']
Numb = df['EC']
IsoUNI = []
for i in UNI:
    xi = i.split(' ')
    for k in xi:
        if k not in IsoUNI: IsoUNI.append(k)

print('PDB finais: '+ str(len(PDB)) + '   UniProt finais: '+ str(len(IsoUNI)))

contagemUNI = {}
for j in IsoUNI:
    if j in ECS:
        Ec = ECS[j][0][0]
        if Ec not in contagemUNI:
            contagemUNI[Ec] = 1
        else: contagemUNI[Ec] += 1
print('UniProt')
print(contagemUNI)

contagemPDB = {}
for i in range(len(PDB)):
    Ec = Numb[i][0]
    if Ec not in contagemPDB:
        contagemPDB[Ec] = 1
    else:
        contagemPDB[Ec] += 1
print('PDBs')
print(contagemPDB)


