import numpy as np 
from multiprocessing import Pool
from functools import partial



def multiply(val1, val2): 
    output = val1*val2
    return output


storeVal = dict()

for ind in np.arange(50): 
    out = multiply(ind, ind)
    storeVal[ind]=out

print(storeVal)

print("pooling...")

p = Pool() 
partial_getIndex = partial(multiply)

results = p.map(partial_getIndex, (np.arange(50), 10) )

print(results)
