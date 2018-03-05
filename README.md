# CodeDemos
Demos of Code! Just starting this on 2018/02/28 to begin turning my in-classroom experience with Python into something that the external world can see. All of my MATLAB code is in my main 2AFC-Analysis repository currently, although I will be gradually adding some of it to this folder over time.

## Table of Contents

### rtidatascience-exercise1 (Python, SQL)
(_Work in progress_) This exercise involves using Python and SQL to turn census data (stored in a normalized relation SQLite database) into a pandas dataframe, then analyze those data and build models to predict whether individuals make over $50,000/year. The README.md file (which describes the task quite comprehensively) and the SQLite database were both provided (courtesy of [the RTI Center for Data Science](https://www.rti.org/service-capability/data-science) in the exercise's Github folder [here](https://github.com/rtidatascience/data-scientist-exercise01).

### rtidatascience-exercise2 (Python)
(*Not yet begun*)

### Probabilistic-Modeling
This goal of this project was to predict human choices as a function of two parallel streams of information: the outcomes of past choices, and visual cues that carry information about how outcome probabilities are changing over time. I built 34 candidate behavioral models, fit data from each subject to each model, and used [Bayesian Model Selection](https://www.ncbi.nlm.nih.gov/pubmed/19306932) to learn which facets of the behavioral models were important for accurate choice prediction. A manuscript from this work is currently under review.

### Data-Processing-Pipeline (MATLAB, Spike2 scripting language)
During my time in graduate school, I led an internal effort within my lab to adopt new hardware and software technologies into our existing data acquisition and processing pipeline. As part of this transition, I had to write multiple pieces of codes in order to seamlessly integrate these new technologies into our existing pipeline. In particular, I had to quickly learn a whole new script language to create the file _intanBatchImport.s2s_, which was written in the [native script language](http://ced.co.uk/products/spkdpsl) used by our recording software (Spike2, Cambridge Electronic Design Limited).

### Classification-GD-ES (MATLAB)
These are MATLAB scripts and functions that I wrote in 2012 to flexibly implement logistic regression (classification) using gradient descent and early stopping to avoid overfitting. The script can accommodate either _leave-one-out_ or K = 3-fold cross-validation.
