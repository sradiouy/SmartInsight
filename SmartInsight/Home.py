import streamlit as st

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

if 'collections' not in st.session_state:
       st.session_state.collections = []

if 'current_collection' not in st.session_state:
       st.session_state.current_collection = None

#collections = st.session_state.get('collections', [])

st.title("Welcome to SmartInsight")
st.subheader("Define settings")

option = st.radio('Pick one', ['Select an existing collection', 'Create new collection'])

if option == "Create new collection":
       st.session_state.current_collection = st.text_input("Enter name for new collection")
       if len(st.session_state.current_collection) > 0:
              if st.session_state.current_collection not in st.session_state.collections:
                     st.session_state.collections.append(st.session_state.current_collection)
              else:
                     st.warning("Collection already created")
elif option == "Select an existing collection":
       st.session_state.current_collection = st.selectbox('Pick an existing collection', st.session_state.collections)  

st.session_state.key = st.text_input("Your Open API Key", "sk...")

submit_button = st.button("Submit")

if submit_button:
       if len(st.session_state.current_collection) > 0:
              st.write(f'Ready to use collection "{st.session_state.current_collection}"')
       else:
              st.warning("Enter a valid name for the collection")

		




