import streamlit as st
import openai
from openai.error import InvalidRequestError
from openai.error import OpenAIError

from typing import List
from streamlit_chat import message
from chatBotContext import *

if "messages" not in st.session_state:
    st.session_state.messages = []

# Clear session state
def clear_chat() -> None:
    st.session_state.generated = []
    st.session_state.past = []
    st.session_state.messages = []
    st.session_state.user_text = ""
    st.success("Chat cleared")

# Display user input area
def show_text_input() -> None:
    st.text_area(label="Write your message:", value=st.session_state.user_text, key="user_text")

# Display chat buttons
def show_chat_buttons() -> None:
    b0, b1, b2 = st.columns(3)
    with b0, b1, b2:
        b0.button(label="Ask")
        b1.button(label="Clear", on_click=clear_chat)
        b2.download_button(
            label="Save",
            data="\\n".join([str(d) for d in st.session_state.messages[1:]]),
            file_name="ai-talks-chat.json",
            mime="application/json",
        )

def create_gpt_completion(ai_model: str, messages: List[dict],text) -> dict:
    """
    Arguments:
    ai_model -- the GPT model
    messages -- a list of previous chat messages

    """
    #st.write("st.session state messages")
    #st.write(messages)
    openai.api_key = st.session_state.key
    completion = openai.ChatCompletion.create(
        model=ai_model,
        messages=[
            {"role": "system", "content": "You are a helpful assistant. That try to use concise answers"},
            {"role": "user", "content": "I need your help to ask question about scientific abstracts"},
            {"role": "assistant", "content": "As a science expert, I am here to help you. Please provide me the list of abstracts so I can assist you:\n"},
            {"role": "user", "content": "I will provide the abstract in the following format \n ID: [ID] Document:[Abstract] Title:[Title] Citation: [Citation] doi: [doi] Year: [year] url: [url]. For example if I give you the following: ID: 148 Document: the sky is black . Url: url1 \n ID: 256 Document: the sky is black and the sun is red. Url: url2\n ID: 31 Document: the sky is blue. Url: url3\n Your answer must be: Based on [148, url1] and [256, url2], the sky is black, however [31, url3] indicates that the sky is blue. So which is the color of the sky?"},
            {"role": "assistant", "content": "Several sources indicates that the color of the sky is black [PubMed ID 148,256], however [Pubmed ID 31] says that is blue."},
            {"role": "user", "content": "Great! Now based in the following context: \n" + ". ".join(messages) + text + ". Please try to cite and the end of a sentence! And if possible obtain information from multiple abstracts. Only answer if information is explicitly in the abstract"}
        ]
        )
    return completion

def show_chat(ai_content: str, user_text: str) -> None:
    """
    This function displays the conversation messages between the AI and the user.

    Arguments:
    ai_content -- the response from the AI
    user_text -- the input text from the user
    """

    if ai_content not in st.session_state.generated:
        # store the ai content
        st.session_state.past.append(user_text)
        st.session_state.generated.append(ai_content)
    if st.session_state.generated:
        for i in range(len(st.session_state.generated)):
            message(st.session_state.past[i], is_user=True, key=str(i) + "_user", avatar_style="pixel-art-neutral")
            #message("", key=str(i))
            message(st.session_state.generated[i], key=str(i), avatar_style="bottts")
            #st.markdown(st.session_state.generated[i])

def show_gpt_conversation() -> None:
    """
    Generate AI response and displaying conversation to the user. 
    """
    try:
        # create the AI response using the GPT model and the user input
        completion = create_gpt_completion(st.session_state.model, st.session_state.messages,st.session_state.user_text)
        #AI response is added to the messages list
        ai_content = completion.get("choices")[0].get("message").get("content")
        st.session_state.messages.append({"role": "assistant", "content": ai_content})
        #If ai_content is not empty, call the show_chat function to display the conversation messages
        if ai_content:
            show_chat(ai_content, st.session_state.user_text)
    except InvalidRequestError as err:
        if err.code == "context_length_exceeded":
            st.session_state.messages.pop(1)
            if len(st.session_state.messages) == 1:
                st.session_state.user_text = ""
            show_conversation()
        else:
            st.error(err)
    except (OpenAIError, UnboundLocalError) as err:
        st.error(err)

def show_conversation(id_lst) -> None:
    #Manage conversation state and update message list
    if st.session_state.messages:
        if id_lst:
            data = query_database(id_lst)
            create_temp_database(data)
            #st.write("result")
            result = query_database_question(st.session_state.user_text)
            #st.write(result)
            context_messages = create_context_messages(result)
            #st.write("context_messages")
            #st.write(context_messages)
            st.session_state.messages = context_messages
        #st.session_state.messages.append({"role": "user", "content": st.session_state.user_text})
        #st.write(st.session_state.messages)
    else:
        ai_role = f"You are a friendly and helpful assistant." 
        st.session_state.messages = [
            {"role": "system", "content": ai_role},
            {"role": "user", "content": st.session_state.user_text},
        ]
        #st.write(ai_role)
    show_gpt_conversation()