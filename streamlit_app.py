import streamlit as st
import os
import pickle
import pandas as pd
from PIL import Image
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from padelpy import padeldescriptor

st.set_page_config(
  page_title='PARP1',
  page_icon='💊',
  layout='wide',
  initial_sidebar_state='expanded')

if 'smiles_input' not in st.session_state:
  st.session_state.smiles_input = ''

if os.path.isfile('molecule.smi'):
  os.remove('molecule.smi') 
  
st.sidebar.title('💊 PARP1 bioactivity prediction')


# Input SMILES
st.sidebar.subheader('Input molecule')

def insert_example_smiles():
    st.session_state.smiles_input = 'CC(=O)OC1=CC=CC=C1C(=O)O'
def clear_smiles():
    st.session_state.smiles_input = ''

smiles_txt = st.sidebar.text_input('Enter SMILES notation', st.session_state.smiles_input)
st.sidebar.button('Example input', on_click=insert_example_smiles)
st.sidebar.button('Clear input', on_click=clear_smiles)


# Default page (loading for the first time)
if st.session_state.smiles_input == '':
  
  with st.expander('About this app'):
    st.write('''
      This QSAR app allow users to predict the biological activity of a query molecule against the target protein being investigated.

      This app is based on the following Python libraries:
      - `streamlit`
      - `pandas`
      - `rdkit`
      - `padelpy`
      ''')
  
  st.subheader('Welcome to the app!')
  st.info('Enter SMILES notation in the sidebar to proceed', icon='👈')
else:
  st.subheader('⚛️ Input molecule:')
  st.info(smiles_txt)
  
  smi = Chem.MolFromSmiles(smiles_txt)
  Chem.Draw.MolToFile(smi, 'molecule.png')
  mol_image = Image.open('molecule.png')
  
  col1, col2, col3 = st.columns(3)
  with col1:
    st.write(' ')
  with col2:
      st.image(mol_image)
  with col3:
      st.write(' ')


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
  
