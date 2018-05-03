import config
import xml.etree.ElementTree as ET
import urllib.request as url
import time
import numpy as np
import pandas as pd

# now build a function to do all of this, given a user id
def pullBooks(user_id):
    
    developer_key = config.DEVELOPER_KEY
    url_root = 'https://www.goodreads.com/review/list/'
    shelf = 'read'
    per_page = '10'
    page = '1'
    webpage = url_root + user_id + '.xml?key=' + developer_key + '&v=2&shelf=' + shelf + '&per_page=' + per_page + '&page=' + page
    xml = ET.parse(url.urlopen(webpage)).getroot()
    
    books = []
    ratings = []
    
    n_books = int(xml.find('reviews').get('total'))
    n_books_per_page = 200
    n_pages = int(np.ceil(n_books/n_books_per_page))
    
    # step through each page and pull data
    for count in range(0,n_pages):
        page = str(count + 1) # determine page number
        per_page = str(n_books_per_page) # 200 books per page
        webpage = url_root + user_id + '.xml?key=' + developer_key + '&v=2&shelf=' + shelf + '&per_page=' + per_page + '&page=' + page
        xml = ET.parse(url.urlopen(webpage)).getroot()
        time.sleep(2) # always wait
    
        # now go through XML and pull out desired information
        for book in xml.findall('reviews/review'):
            books.append(book.find('book/title_without_series').text)
            ratings.append(book.find('rating').text)
            
    d = {'user' : user_id,
     'title' : books,
     'rating' : ratings,
     }

    #import pdb; pdb.set_trace()
    #df = pd.DataFrame(data=d).drop_duplicates(subset='title').pivot(index='title', columns='user', values='rating')
    df = pd.DataFrame(data=d)
    print(df.shape)        
    return df