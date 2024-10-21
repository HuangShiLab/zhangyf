# Import necessary modules
import numpy as np
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import auc

def AUPR(expected_abd, observed_abd):
    expected_abd = np.array(expected_abd)
    observed_abd = np.array(observed_abd)

    # Convert expected abundance values to a binary array
    bi_expected_abd = np.array(expected_abd != 0, dtype=int)

    # Calculate precision and recall values at different thresholds using precision_recall_curve
    precision, recall,_ = precision_recall_curve(bi_expected_abd, observed_abd)

    # Calculate the Area Under the Precision-Recall Curve (AUPR) using the auc function
    aupr = auc(recall, precision)

    return (aupr)
