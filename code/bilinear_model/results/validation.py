# -*- coding: utf-8 -*-
"""
Created on Sat Mar 18 01:25:55 2017

@author: Achin Jain, mLAB, UPenn
"""
from __future__ import division
import numpy as np
import pandas as pd
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
from sklearn.metrics import explained_variance_score

y = pd.read_csv('validation1.csv', header=None)
ytrue = y.iloc[:,0]
ypredrt = y.iloc[:,1]
ypredrf = y.iloc[:,2]

print('----------------------------------------------------------------------')

MSE = mean_squared_error(ytrue, ypredrt, multioutput='raw_values')
R2Score = r2_score(ytrue, ypredrt, multioutput='raw_values')
EV = explained_variance_score(ytrue, ypredrt, multioutput='raw_values')
print('root mean square error: %s' %(np.sqrt(MSE)))
print('normalized mean square error: %s' %(np.sqrt(MSE)/np.array(np.mean(ytrue))))
print('R2 score: %s' %(R2Score))
print('explained variance: %s' %(EV))

print('----------------------------------------------------------------------')

MSE = mean_squared_error(ytrue, ypredrf, multioutput='raw_values')
R2Score = r2_score(ytrue, ypredrf, multioutput='raw_values')
EV = explained_variance_score(ytrue, ypredrf, multioutput='raw_values')
print('root mean square error: %s' %(np.sqrt(MSE)))
print('normalized mean square error: %s' %(np.sqrt(MSE)/np.array(np.mean(ytrue))))
print('R2 score: %s' %(R2Score))
print('explained variance: %s' %(EV))

print('----------------------------------------------------------------------')

y = pd.read_csv('validation6.csv', header=None)
ytrue = y.iloc[:,0]
ypredrt = y.iloc[:,1]
ypredrf = y.iloc[:,2]

print('----------------------------------------------------------------------')

MSE = mean_squared_error(ytrue, ypredrt, multioutput='raw_values')
R2Score = r2_score(ytrue, ypredrt, multioutput='raw_values')
EV = explained_variance_score(ytrue, ypredrt, multioutput='raw_values')
print('root mean square error: %s' %(np.sqrt(MSE)))
print('normalized mean square error: %s' %(np.sqrt(MSE)/np.array(np.mean(ytrue))))
print('R2 score: %s' %(R2Score))
print('explained variance: %s' %(EV))

print('----------------------------------------------------------------------')

MSE = mean_squared_error(ytrue, ypredrf, multioutput='raw_values')
R2Score = r2_score(ytrue, ypredrf, multioutput='raw_values')
EV = explained_variance_score(ytrue, ypredrf, multioutput='raw_values')
print('root mean square error: %s' %(np.sqrt(MSE)))
print('normalized mean square error: %s' %(np.sqrt(MSE)/np.array(np.mean(ytrue))))
print('R2 score: %s' %(R2Score))
print('explained variance: %s' %(EV))

print('----------------------------------------------------------------------')