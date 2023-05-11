import streamlit as st
import openai
from openai import api_key
import chromadb
import pandas as pd

from search_update_db import *

st.set_page_config(layout="wide", page_title="SmartInsight")

hide_default_format = """
       <style>
       footer {visibility: hidden;}
       </style>
       """
	   #MainMenu {visibility: hidden; }
st.markdown(hide_default_format, unsafe_allow_html=True)

if 'key' not in st.session_state:
    st.session_state.key = None

if 'analysis' not in st.session_state:
       st.session_state.analysis = []

if 'current_analysis' not in st.session_state:
       st.session_state.current_analysis = None

st.title("Welcome to SmartInsight")
st.subheader("Define settings")

st.session_state.key = st.text_input("Your Open API Key", "sk...")

if st.session_state.key is not None:
       try:
              st.session_state.analysis = find_analysis()
       except:
              st.error("Enter a valid key")


option = st.radio('Pick one', ['Select an existing analysis', 'Create new analysis'])

if option == "Create new analysis":
       new_analysis = st.text_input("Enter name for new analysis")
       if new_analysis is not None and ';' not in new_analysis:
              # Valid input: update the current analysis state
              st.session_state.current_analysis = new_analysis
       else:
              # Invalid input: show an error message
              st.error("Invalid analysis name. Please enter a name that does not contain semicolon (;) character.")

elif option == "Select an existing analysis":
       st.session_state.current_analysis = st.selectbox('Pick an existing analysis', st.session_state.analysis)  



submit_button = st.button("Submit")

if submit_button:
       if option == "Create new analysis" and st.session_state.current_analysis in st.session_state.analysis:
              st.warning("Analysis already exists")
       elif len(st.session_state.current_analysis) > 0:
              #st.session_state.analysis.append(st.session_state.current_analysis)
              st.write(f'Ready to use analysis "{st.session_state.current_analysis}"')
       else:
              st.warning("Enter a valid name for the analysis")

       #Agregar validación de la key 



		




