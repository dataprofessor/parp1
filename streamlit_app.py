import streamlit as st
import os
import pickle
import pandas as pd
from PIL import Image
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from padelpy import padeldescriptor
from pathlib import Path

APP_DIR = Path(__file__).parent

st.set_page_config(
  page_title='PARP1pred',
  page_icon='💊',
  initial_sidebar_state='expanded')

if 'smiles_input' not in st.session_state:
  st.session_state.smiles_input = ''
if 'submitted' not in st.session_state:
  st.session_state.submitted = False

if os.path.isfile(APP_DIR / 'molecule.smi'):
  os.remove(APP_DIR / 'molecule.smi') 

st.title('💊 PARP1pred app')
st.info('PARP1pred allow users to predict whether a query molecule is active/inactive towards the PARP1 target protein.')

tab1,tab2,tab3,tab4,tab5,tab6,tab7 = st.tabs(['Main', 'About', 'What is PARP1?', 'Dataset', 'Model performance', 'Python libraries', 'Citing us'])

with tab1:
  with st.form('my_form'):
    st.subheader('Predict PARP1 inhibitory activity')
    smiles_txt = st.text_input('Enter SMILES notation', st.session_state.smiles_input)
    with st.expander('Example SMILES'):
      st.code('O=C(c1cc(Cc2n[nH]c(=O)c3ccccc23)ccc1F)N1CCN(C(=O)C2CC2)CC1')
    submit_button = st.form_submit_button('Submit')

  if submit_button and smiles_txt:
    st.session_state.smiles_input = smiles_txt
    st.session_state.submitted = True

  if st.session_state.submitted and st.session_state.smiles_input:
    smiles = st.session_state.smiles_input
    
    st.subheader('⚛️ Input molecule:')
    with st.expander('Show SMILES', expanded=True):
      st.text(smiles)

    with st.expander('Show chemical structures', expanded=True):
      smi = Chem.MolFromSmiles(smiles)
      Chem.Draw.MolToFile(smi, str(APP_DIR / 'molecule.png'), width=900)
      mol_image = Image.open(APP_DIR / 'molecule.png')
      st.image(mol_image)

    f = open(APP_DIR / 'molecule.smi', 'w')
    f.write(f'{smiles}\tmol_001')
    f.close()

    st.subheader('🔢 Descriptors')
    with st.spinner('Calculating molecular descriptors...'):
      if os.path.isfile(APP_DIR / 'molecule.smi'):
        padeldescriptor(mol_dir=str(APP_DIR / 'molecule.smi'), 
                          d_file=str(APP_DIR / 'descriptors.csv'),
                          descriptortypes=str(APP_DIR / 'data/PubchemFingerprinter.xml'), 
                          detectaromaticity=True,
                          standardizenitro=True,
                          standardizetautomers=True,
                          threads=2,
                          removesalt=True,
                          log=True,
                          fingerprints=True)

    descriptors = pd.read_csv(APP_DIR / 'descriptors.csv')
    descriptors.drop('Name', axis=1, inplace=True)

    with st.expander('Show full set of descriptors as calculated for query molecule'):
      st.write(descriptors)
      st.write(descriptors.shape)

    import warnings
    with warnings.catch_warnings():
      warnings.simplefilter("ignore")
      model = pickle.load(open(APP_DIR / 'data/original.pkl', 'rb'))
    pubchem_subset = model.feature_names_in_

    query_desc_1 = descriptors.columns.difference(pubchem_subset)
    query_desc_2 = descriptors.drop(query_desc_1, axis=1)

    with st.expander('Show subset of descriptors as used in trained model'):
      st.write(query_desc_2)
      st.write(query_desc_2.shape)

    st.subheader('🤖 Predictions')
    pred = int(model.predict(query_desc_2))
    if pred == 0:
      st.error('Inactive')
    if pred == 1:
      st.success('Active')

with tab2:
  coverimage = Image.open(APP_DIR / 'PARP1pred.jpg')
  st.image(coverimage)
with tab3:
  st.header('What is PARP1?')
  st.write('Poly (ADP-ribose) polymerase-1 (PARP-1) is an enzyme that catalyzes the ADP-ribosylation of a specific protein and plays a vital role in DNA repair. It has become an attractive target as inhibition of PARP-1 causes a toxic accumulation of DNA double strand breaks in cancer cells, particularly with BRCA1/2 deficiency, found in breast, ovarian, prostate, and pancreatic cancers.')
with tab4:
  st.header('Dataset')
  st.write('''
    In our work, we retrieved a human PARP-1 biological dataset from the ChEMBL database. The data was curated and resulted in a non-redundant set of 2,018 PARP-1 inhibitors, which can be divided into:
    - 1,720 active compounds
    - 298 inactive compounds
    ''')
with tab5:
  st.header('Model performance')
  st.write('We selected PubChem as a molecular fingerprint and used a random forest with an oversampling approach to construct the best model. The Matthews correlation coefficients for training, cross-validation, and test sets were 1.00, 0.96, and 0.74, respectively.')
with tab6:
  st.header('Python libraries')
  st.markdown('''
    This app is based on the following Python libraries:
    - `streamlit`
    - `pandas`
    - `rdkit`
    - `padelpy`
  ''')
with tab7:
  st.markdown('T. Lerksuthirat, S. Chitphuk, W. Stitchantrakul, D. Dejsuphong, A.A. Malik, C. Nantasenamat, PARP1PRED: A web server for screening the bioactivity of inhibitors against DNA repair enzyme PARP-1, ***EXCLI Journal*** (2023) DOI: https://doi.org/10.17179/excli2022-5602.')
