import streamlit as st
import os
import pickle
import pandas as pd
from PIL import Image
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from padelpy import padeldescriptor

# Page configuration
st.set_page_config(
  page_title='PARP1pred',
  page_icon='💊',
  initial_sidebar_state='expanded')

# Session state
if 'smiles_input' not in st.session_state:
  st.session_state.smiles_input = ''

if os.path.isfile('molecule.smi'):
  os.remove('molecule.smi') 
  
st.sidebar.title('💊 PARP1pred')


# Input SMILES
st.sidebar.subheader('Input molecule')

def insert_example_smiles():
    st.session_state.smiles_input = 'O=C(c1cc(Cc2n[nH]c(=O)c3ccccc23)ccc1F)N1CCN(C(=O)C2CC2)CC1'
def clear_smiles():
    st.session_state.smiles_input = ''

smiles_txt = st.sidebar.text_input('Enter SMILES notation', st.session_state.smiles_input)
st.sidebar.button('Example input', on_click=insert_example_smiles)
st.sidebar.button('Clear input', on_click=clear_smiles)


# Default page (loading for the first time)
if st.session_state.smiles_input == '':
  st.subheader('Welcome to the app!')
  
  with st.expander('About this app'):
    st.write('''
      PARP1pred allow users to predict whether a query molecule is active/inactive towards the PARP1 target protein.

      This app is based on the following Python libraries:
      - `streamlit`
      - `pandas`
      - `rdkit`
      - `padelpy`
      ''')
  
  st.info('Enter SMILES notation in the sidebar to proceed', icon='👈')
  
else:
  st.subheader('⚛️ Input molecule:')
  st.write('**SMILES**')
  st.text(smiles_txt)
  
  st.write('**Chemical structure**')
  smi = Chem.MolFromSmiles(smiles_txt)
  Chem.Draw.MolToFile(smi, 'molecule.png', width=900)
  mol_image = Image.open('molecule.png')
  st.image(mol_image)

# Input SMILES saved to file
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
  st.write('**Full set of descriptors (calculated for query molecule)**')
  st.write(descriptors)
  st.write(descriptors.shape)
  

# Load descriptor subset used in trained model
if st.session_state.smiles_input != '':
  model = pickle.load(open('data/oversampling_PubChem_RandomForestClassifier.pkl', 'rb'))
  pubchem_subset = model.feature_names_in_

  query_desc_1 = descriptors.columns.difference(pubchem_subset)
  query_desc_2 = descriptors.drop(query_desc_1, axis=1)

  st.write('**Subset of descriptor (used in trained model)**')
  st.write(query_desc_2)
  st.write(query_desc_2.shape)

  
# Read in saved classification model
if st.session_state.smiles_input != '':
  st.subheader('🤖 Predictions')
  pred = int(model.predict(query_desc_2))
  if pred == 0:
    st.error('Inactive')
  if pred == 1:
    st.success('Active')
  
