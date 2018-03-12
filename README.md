# pcaPA

The R package **pcaPA** is a set of functions to perform parallel analysis for
principal components analysis intended mainly for large data sets. It performs a 
parallel analysis of continuous, ordered (including dichotomous/binary as a special case) 
or mixed type of data associated with a principal components analysis.  
Polychoric correlations among ordered variables, Pearson correlations among continuous 
variables and polyserial correlation between mixed type variables (one ordered and one continuous)
are used. Whenever the use of polyserial or polychoric correlations yields a non positive
definite correlation matrix, the resulting matrix is transformed into the nearest positive 
definite matrix. This is a continued work based on a previous version developed at the Colombian 
Institute for the evaluation of education - ICFES

