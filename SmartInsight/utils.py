import random
import string

def generate_random_emails():
    """Generates random emails.
    
    Returns:
        [string] -- Random email.
    """

    return ''.join(random.choice(string.ascii_lowercase[:12]) for i in range(7)) + '@' + random.choice([ "hotmail.com", "gmail.com", "aol.com", "mail.com" , "mail.kz", "yahoo.com"])


def convert_df(df):
    # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return df.to_csv(index=False).encode('utf-8')