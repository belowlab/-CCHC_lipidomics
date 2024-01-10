# This code is to train boosting models on CCHC lipidomics data.
# TODO
# Implement boosting models
# 1. AdaBoost
# 2. Gradient boost
# 3. XGBoost
# Options of base estimator: Linear regression Multilayer perceptron (MLP), Decision tree regressor

# Refer to code ML_03_model_train.py and ML_07_boosting_models_test_run.ipynb

from sklearn.ensemble import AdaBoostRegressor
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.datasets import make_regression
import logging
import matplotlib.pyplot as plt
import subprocess
from scipy import stats

import pandas as pd
import numpy as np
import os
import datetime
msg = '', datetime.datetime.now().strftime('%Y-%m-%d')
print('', datetime.datetime.now().strftime('%Y-%m-%d'))

logging.getLogger().setLevel(logging.INFO)
logging.basicConfig(level=logging.INFO)