import csv

def csv_out(lista, out):
    '''
    :return: cria ficheiro .csv com os outputs da conversão
    '''
    file = open('textfiles/'+out+'.csv', 'w', newline='')
    header = ['UniProt', "PDB",'Vetor','EC']
    writer = csv.DictWriter(file, fieldnames=header)
    writer.writeheader()
    for i in lista:
        Uni,PDB, vetor,EC = i
        writer.writerow({'UniProt': Uni, "PDB": PDB,'EC': EC, 'Vetor': vetor })
