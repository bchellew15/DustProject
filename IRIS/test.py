import numpy as np

tempo = np.ones(100)
weight = np.zeros(100)
result = np.zeros(100)
indw = np.where(tempo != -32768)[0]

for i in range(2):
    weight[indw] = weight[indw]+1
    result[indw] = result[indw] + tempo[indw]

print("weight:", weight)
print("result:", result)
    
indw = np.where(weight > 0)[0] #nbindw, complement...
mask = np.zeros(100)
mask[indw] = True
indempty = mask==0
nbindw = indw.shape[0]

result[indw] = result[indw] / weight[indw]

print("result:", result)
