import pandas as pd
import tqdm
import pickle
import requests
from rdkit.Chem import MACCSkeys
from rdkit.Chem import AllChem
from rdkit import DataStructs, Chem


def generate_fingerprint(data, smile_name):
    '''
    :param data: the original table with smiles columns
    :param smile_name: give the column name that represent molecular smiles
    :return: the extend fingerprint dataframe and maccs fingerprint dataframe
    '''

    # remove records rows whose smiles is NA or blank
    data = data[(data[smile_name] != 'NA') & (data[smile_name].notnull()) & (data[smile_name].notna())]

    # set fingerrpint parameter
    radius = 2
    nBits = 1024

    # generate fingerprint on by one
    ext = []
    mac = []

    for s in tqdm.tqdm(list(data['smiles'])):

        # generate extend fingerrpint by smiles
        mol = Chem.MolFromSmiles(s)
        finger_ext = AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=nBits).ToBitString()
        finger_ext = list(finger_ext)
        ext.append([s]+finger_ext)

        # generate maccs fingerrpint by smiles
        finger_maccs = MACCSkeys.GenMACCSKeys(mol).ToBitString()
        finger_maccs = list(finger_maccs)
        mac.append([s] + finger_maccs)

    # set result dataframe column names
    col_ext = [smile_name] + ['EXTENDED_' + str(i) for i in range(1, 1025)]
    col_mac = [smile_name] + ['MACCS_' + str(i) for i in range(1, 168)]

    # prepare result to dataframe and merge with original data
    ext_pd = pd.DataFrame(ext, columns=col_ext)
    ext_pd_merged = pd.merge(data, ext_pd, left_on=smile_name, right_on=smile_name, how='left')
    maccs_pd = pd.DataFrame(mac, columns=col_mac)
    maccs_pd_merged = pd.merge(data, maccs_pd, left_on=smile_name, right_on=smile_name, how='left')

    return ext_pd_merged, maccs_pd_merged


def main():
    ctrp = pd.read_csv('drug_id_list/ctrp_idlist.csv')
    ctrp.name = 'ctrp'
    gdsc = pd.read_csv('drug_id_list/gdsc2_idlist.csv')
    gdsc.name = 'gdsc'
    prism = pd.read_csv('drug_id_list/prism_idlist.csv')
    prism.name = 'prism'

    # Loop each dataset and save out as csv
    for data in [ctrp, gdsc, prism]:
        ext_pd, maccs_pd = generate_fingerprint(data, 'smiles')
        ext_pd.to_csv('result/extended_finger_{}.csv'.format(data.name), index = None)
        maccs_pd.to_csv('result/maccs_finger_{}.csv'.format(data.name), index = None)