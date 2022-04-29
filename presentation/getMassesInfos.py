import numpy as np

massesF=np.load("finalMasses.npy")
masseInit=np.load("initMasses.npy")

#print(masseInit)
#i = input()
print("Masse de la Terre:" , "finale: " ,massesF[0],"initiale: ", masseInit[0])

print("Masse du plus grand corps orbitant: ", "finale: ", np.max(massesF[1:]), "initiale: ", np.max(masseInit[1:]))