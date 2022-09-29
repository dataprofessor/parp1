import streamlit as st
import os
import pickle
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from stmol import speck_plot
from padelpy import padeldescriptor

st.set_page_config(
  page_title='QSAR',
  page_icon='💊',
  layout='wide',
  initial_sidebar_state='collapsed')

if 'smiles_input' not in st.session_state:
  st.session_state.smiles_input = ''

if os.path.isfile('molecule.smi'):
  os.remove('molecule.smi') 
  
st.title('💊 QSAR app')

with st.expander('About this app'):
  st.write('''
    This QSAR app allow users to predict the biological activity of a query molecule against the target protein being investigated.
    
    This app is based on the following Python libraries:
    - `streamlit`
    - `pandas`
    - `rdkit`
    - `padelpy`
    ''')

# Input SMILES
st.sidebar.subheader('Input molecule')

def insert_example_smiles():
    st.session_state.smiles_input = 'CC(=O)OC1=CC=CC=C1C(=O)O'
def clear_smiles():
    st.session_state.smiles_input = ''

smiles_txt = st.sidebar.text_input('Enter SMILES notation', st.session_state.smiles_input)
st.sidebar.button('Example input', on_click=insert_example_smiles)
st.sidebar.button('Clear input', on_click=clear_smiles)


if st.session_state.smiles_input == '':
  st.subheader('Welcome to the app!')
  st.info('Enter SMILES notation in the sidebar to proceed', icon='👈')
else:
  st.subheader('⚛️ Input molecule:')
  st.info(smiles_txt)

f = open('molecule.smi', 'w')
f.write(f'{smiles_txt}\tmol_001')
f.close()


# Compute PADEL descriptors
if st.session_state.smiles_input != '':
  st.subheader('🔢 Descriptors')
  if os.path.isfile('molecule.smi'):
    padeldescriptor(mol_dir='molecule.smi', 
                    d_file='descriptors.csv',
                    descriptortypes='data/PubchemFingerprinter.xml', 
                    detectaromaticity=True,
                    standardizenitro=True,
                    standardizetautomers=True,
                    threads=2,
                    removesalt=True,
                    log=True,
                    fingerprints=True)

  descriptors = pd.read_csv('descriptors.csv')
  descriptors.drop('Name', axis=1, inplace=True)
  st.write(descriptors)
  

# Load descriptor subset used in trained model
pubchem = pd.read_csv('data/PubChem_removed_low_variance_0.1_2class_chemblID.csv')
pubchem = pubchem.drop(['molecule_chembl_id', 'class'], axis=1)
st.write(pubchem.columns)
  
# Read in saved classification model
if st.session_state.smiles_input != '':
  model = pickle.load(open('data/oversampling_PubChem_RandomForestClassifier.pkl', 'rb'))
  importances = pd.Series(model.feature_importances_)
  st.write(importances)
