import streamlit as st
import openai
from streamlit_option_menu import option_menu

from chatBot import *
from chatBotContext import *
from updateDb import find_project

def chat_page():
    with st.sidebar:
        st.write("**ChatBot Settings**")
        st.selectbox(label= "Select AI model", key="model", options=["gpt-3.5-turbo","davinci", "curie", "babbage", "ada"])
        
        option = st.radio("Interact with a project or a cluster:",options=["Project","Cluster"])
        if option == "Cluster":
            all_analysis = find_analysis()
            name = st.selectbox(label="Select analysis",options=all_analysis)
            topic_dict = find_cluster(name)
            topics = list(topic_dict.keys())
            selected_topic = st.selectbox("Select a topic", topics)
            id_lst = topic_dict[selected_topic]

        elif option == "Project":
            project_lst = find_project()
            if project_lst is not None:
                selection = st.selectbox('Select project', project_lst, label_visibility="collapsed")
                id_lst = query_database_project(selection)

    if option == "Cluster":
        st.write(f"Using topic '{selected_topic}' inside analysis '{name}'")
    elif option == "Project":
        st.write(f"Using project '{selection}'")

    if st.session_state.user_text:
        #st.write(st.session_state.messages)
        show_conversation(id_lst)
        st.session_state.user_text = ""
    show_text_input()
    show_chat_buttons()
    