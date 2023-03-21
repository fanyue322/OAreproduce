import sklearn
from sklearn.linear_model import LogisticRegression

def LogisticClassification(x,y):
    clf=LogisticRegression(penalty="elasticnet",solver="saga",l1_ratio=0.5).fit(x,y)
    return(clf)
    



