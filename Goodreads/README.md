## Goodreads

Hi! This folder contains Python code to scrape book ratings from Goodreads and use those ratings to predict what books a user may like (collaborative filtering). This is very similar to the "Netflix Prize", where users were given (movie, user, rating) triplets and asked to predict the ratings for unobserved movie-user pairings: but with books! This project was mostly borne out of a desire to find new books to read :)

### *goodreads.ipynb*
This Jupyter notebook contains the main code used for the project. The notebook is pretty well-commented, so I won't dive down into it too much here. In the code, I use non-metric multidimensional scaling (MDS) to show why nearest-neighbors doesn't work well when the dataset is very sparse, and then I use incremental singular-value decomposition (SVD) to perform matrix factorization and predict ratings for novel book-user pairings.

### *goodreadsFns.py*
This file contains the function used scrape all of the books and ratings for a given user from the Goodreads website (XML).

### Data: *book_ratings.csv*, *distance.pkl*, and *data.pkl*
I've included some data files so that the more time-intensive parts of the code don't need to be run everytime the notebook is restarted. *book_ratings.csv* contains 2715 (book, user, rating) triplets from 13 different users. This dataset is tiny compared to the 100 million ratings that Netflix provided, but it should be large enough for a proof-of-concept. Currently, I'm scraping all of this ratings from the web. In the future, I plan on using the Goodreads API (which is non-trivial to work with in Python) to pull more ratings faster. *distance.pkl* contains the distances between pairs of books in rating-space, which is used in MDS. *data.pkl* contains various outputs from incremental SVD.
