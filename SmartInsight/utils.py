import random
import string
import streamlit as st

def generate_random_emails():
    """Generates random emails.
    
    Returns:
        [string] -- Random email.
    """

    return ''.join(random.choice(string.ascii_lowercase[:12]) for i in range(7)) + '@' + random.choice([ "hotmail.com", "gmail.com", "aol.com", "mail.com" , "mail.kz", "yahoo.com"])

@st.cache_data
def convert_df(df):
    return df.to_csv(index=False).encode('utf-8')