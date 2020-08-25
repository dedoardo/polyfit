import pickle
import numpy as np


def generate_array(dims, data):
    assert np.prod(dims)  == len(data)
    return np.array(data).reshape(*dims)

def flatten_array(array):
    return list(array.shape), list(array.flatten())

def fprint_object(obj, file):
    with open(fname, "w") as ff:    
        ff.write("%s \n"%obj)

def print_object(obj):
    print( obj )

def pickle_object(obj, fname):
    with open(fname, "wb") as ff:
        return pickle.dump( obj, file=ff, protocol=3 )

def unpickle_object(fname):
    with open(fname, "rb") as ff:
        return pickle.load( file=ff )
    
    
