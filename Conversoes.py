from Database_converter import Convert
from PDB_info import PDB_batch

########## DATASET ECPRED ###########

ECPRED = Convert('textfiles/ecpred_uniprot_uniref_90.csv','uniref_90','ACC+ID','PDB_ID', auto= False)
dicECPRED = ECPRED.organize_dic()

batch = PDB_batch(dicECPRED,'textfiles/All_ECPRED')
file = batch.test_txt('textfiles/ECPRED')
filebatch = batch.PDB_downloadbatch()

########## DATASET AMIDI ###########

Amidi = Convert('textfiles/Dataset39251.csv','PDB','PDB_ID','ACC', auto= False)
dicAmidi = Amidi.organize_dic()

dicAmidi2 = {}
for i in dicAmidi:
    for j in dicAmidi[i]:
        if j not in dicAmidi2:
            dicAmidi2[j] = [i]
        else:
            dicAmidi2[j] = dicAmidi2[j] + [i]
print(dicAmidi2)
batch = PDB_batch(dicAmidi2,'textfiles/All_AMIDI')
file = batch.test_txt('textfiles/AMIDI')
filebatch = batch.PDB_downloadbatch()
