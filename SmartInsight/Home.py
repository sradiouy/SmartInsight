import streamlit as st
from streamlit_extras.stateful_button import button

import openai
from openai import api_key
import chromadb
import pandas as pd
from PIL import Image

from updateDb import *
from utils import validate_key

#im = Image.open('/AppIcon.png')
#st.set_page_config(layout="wide", page_title="SmartInsight",page_icon = im)
st.set_page_config(layout="wide", page_title="SmartInsight")

hide_default_format = """
       <style>
       footer {visibility: hidden;}
       </style>
       """
	   #MainMenu {visibility: hidden; }
st.markdown(hide_default_format, unsafe_allow_html=True)


# Initialize session state variables
if 'key' not in st.session_state:
    st.session_state.key = None

if 'project' not in st.session_state:
       st.session_state.project = None

if 'current_project' not in st.session_state:
       st.session_state.current_project = None



st.title("Welcome to SmartInsight")

st.write("SmartInsight is a comprehensive research platform designed to streamline scientific literature analysis. Powered by AI and machine learning algorithms, SmartInsight provides scientists and researchers with efficient tools to explore and extract valuable insights from diverse resources.")


st.markdown("""---""")


# Sidebar 
st.sidebar.title("Settings")

# Radio button in the sidebar to navigate between pages.
selected_page = st.sidebar.radio("", ["Define OpenAI API Key","Select Project"])

st.subheader("Settings")
st.write("\n")

if selected_page == "Define OpenAI API Key":
       st.write("**Define OpenAI API Key**")
       if st.session_state.key is not None and validate_key():
              st.success("OpenAI API Key already in use. To change it enter a new one")
       key = st.text_input("Your OpenAI API Key",label_visibility="collapsed")
       if st.button("Set key"):
              st.session_state.key = key
              if validate_key():
                     st.success("Ready to use OpenAI API Key")
              else:
                     st.error("Invalid OpenAI API Key")
       

elif selected_page == "Select Project":

       st.write("**Select project to work with SmartInsight**")
       if validate_key():
              option = st.radio('', ["Choose existing project", "Create new project"],label_visibility="collapsed")

              if option == "Choose existing project":
                     project_lst = find_project()
                     if project_lst is not None:
                            selection = st.selectbox('Select project', project_lst, label_visibility="collapsed")
                     elif not validate_key():
                            st.error("Invalid OpenAI API Key")
                     else:
                            st.error("No projects available. Create a new project")
                     
                     done = st.button("Select")

                     if selection is not None and done:
                            st.session_state.current_project = selection
                            st.success(f"Ready to use {st.session_state.current_project}")
                     elif done:
                            st.error("Choose a project to continue")

              elif option == "Create new project":
                     new_project = st.text_input("Enter name for new project")
                     if new_project is not None and ';' not in new_project:
                            # Valid input: update the current project state
                            st.session_state.current_project = new_project
                     else:
                            # Invalid input: show an error message
                            st.error("Invalid project name. Please enter a name that does not contain semicolon (;) character.")
                     
                     st.session_state.project = find_project()
                     
                     if st.button("Create"):
                            if st.session_state.current_project in st.session_state.project:
                                   st.warning("Project already exists")
                            elif len(st.session_state.current_project) > 0:
                                   #st.session_state.project.append(st.session_state.current_project)
                                   st.success(f"Ready to use project '{st.session_state.current_project}'")
                            else:
                                   st.warning("Enter a valid name for the project")

       else:
              st.error("Enter a valid OpenAI API Key")
		




