from read_file import *
import numpy as np
import iminuit as mi
from numpy.linalg import svd
from  sklearn.decomposition import PCA

def img2alpha(X,Y) :
    def errXY(alpha) :
        return np.sum(dot_mul(gate(X,alpha)-Y,gate(X,alpha)-Y))
    M=mi.Minuit(errXY,alpha=0.5,errordef=1,error_alpha=0.001,limit_alpha=(0,0.95))
    return M.values['alpha']


def modify_data(X) :
    XX=[]
    for i in X :
        XX.append(np.concatenate(i,axis=None))
    XY=np.array(XX)
    return XY

def dimension_red(X) :
    pca=PCA(n_components=20)
    Xp=pca.fit_transform(X)
    return Xp

def training_model(X,Y) :
    Xtrain=X[:int(0.8*len(X))]
    Ytrain = Y[:int(0.8 * len(X))]
    Xtest = X[int(0.8 * len(X)):]
    Ytest = Y[int(0.8 * len(X)):]
    YYtrain=np.zeros(len(Ytrain))
    XXtrain=np.zeros(len(Xtrain))
    for i in range(len(Ytrain)) :
        YYtrain[i]=img2alpha(Xtrain[i],Ytrain[i])
        XXtrain[i]=dimension_red(modify_data(Xtrain[i]))
    theta=np.dot(np.invert(np.dot(XXtrain.T,XXtrain)),np.dot(XXtrain.T,YYtrain))
    Yp=np.dot(Xtest,theta)
    print("score : {}".format(np.sum(dot_mul(Yp-Ytest,Yp-Ytest))))
    return theta

def predict(data,Bscan) :
    X,Y=data
    theta=training_model(X,Y)
    alpha=np.dot(dimension_red(Bscan),theta)
    return gate(Bscan,alpha)

