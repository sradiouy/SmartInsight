from utils import get_keepterms, generate_stop_grams, nltk
from nltk import FreqDist
from wordfreq import zipf_frequency
import string
from collections import Counter
import pandas as pd
import numpy as np

def get_words(abstractdict):
    """ Tokenization of articles abstracts. Breaking up sequences of strings into pieces (words/terms). For each abstract remove punctuation and transform to lower case. 
    
    Arguments:
        abstractdict {[dict]} -- Dictionary of articles, where the key is the gene and the value is the combined abstract 
    
    Returns:
        [type] -- Dictionary of articles, where the key is the gene and the value is a list of tokens.  
    """

    worddict = {} #empty dict
    removepunct = string.punctuation.replace("-","").replace("_","") #keep - and _ in text other punctuation removed
    table = str.maketrans('', '', removepunct) # create table to make the replacements
    for article in abstractdict: 
        worddict[article] = ([*map(str.lower, nltk.word_tokenize(abstractdict[article].translate(table)))]) #Tokenization and pre-processing
    return worddict

def cleanwords(abstractdict,keep,remove,zipf):
    """ Transform the list of word of each articles into a cleaned list (words without common terms, and trying to reduce the amount of words via lemmatization and stemmatization)
    
    Arguments:
        abstractdict {[dict]} -- Dictionary of articles, where the key is the pubmed id and the value is the abstract
    
    Returns:
        [dict] --Dictionary of articles, where the key is the pubmed id and the value is a list of cleaned (remove common words, stemm, lemma, etc.) tokens 
    """
    stopgrams = generate_stop_grams()
    stopgrams += remove
    keepgrams = get_keepterms()
    keepgrams += keep
    worddict = get_words(abstractdict)
    lemmatizer = nltk.stem.WordNetLemmatizer() #lemmatizer
    for article in worddict.keys():
        cleaned_words = [w for w in worddict[article] if (len(w) > 2 and not w[0].isdigit() and not w[-1] == "-")] # remove tokens if len < 3, start with a number or end with -.
        cleaned_words = [w for w in cleaned_words if w not in stopgrams] # remove common words 
        cleaned_words = [lemmatizer.lemmatize(l,pos="v") for l in cleaned_words] #transform to verb form
        cleaned_words = [lemmatizer.lemmatize(l) if l[-1] == "s" else l for l in cleaned_words] #try to remove plurals from not s ending words.
        cleaned_words = [w for w in cleaned_words if (zipf_frequency(w, 'en') <= zipf or w in keepgrams)] #keep word if zipf score is low (not frequent in english) or if it is in keep terms (words common in biology))
        cleaned_words = [w for w in cleaned_words if (len(w) > 2 and not w in stopgrams)]
        worddict[article] = cleaned_words
    return worddict

def countwords(cleaned_tokens):
    gene_term_counts = {}

    # Iterate over the genes and their tokens
    for gene, tokens in cleaned_tokens.items():
        # Count the occurrences of each term
        term_counts = Counter(tokens)
        
        # Filter out terms that appear only once
        filtered_counts = {term: count for term, count in term_counts.items() if count > 1}
        
        # Store the gene and its term counts in the resulting dictionary
        gene_term_counts[gene] = filtered_counts

    rows = []
    # Print the gene-term count dictionaries
    for gene, term_counts in gene_term_counts.items():
        for term, weight in term_counts.items():
            rows.append({'gene': gene, 'term': term, 'actual_weight': weight})

    df = pd.DataFrame(rows)

    total_weights = df.groupby('gene')['actual_weight'].sum()

    # Normalize weights for each gene
    df['weight'] = df.apply(lambda row: np.round((row['actual_weight'] / total_weights[row['gene']]) * 100, 2), axis=1)


    # Filter keeping terms that appear at least 2 times
    filtered_df = df[df['term'].map(df['term'].value_counts()) >= 2]
    
    return filtered_df
