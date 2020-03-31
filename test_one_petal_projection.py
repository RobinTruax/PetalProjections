import spp_library

petal_permutation = [1, 3, 5, 2, 7, 4, 6]
petal_number = len(petal_permutation)
unsigned = spp_library.unsignedGaussCode(petal_number)
signed = spp_library.signGaussCode(unsigned, petal_permutation)
arcs = spp_library.splitSignedGaussCode(signed)
crossing_matrix = spp_library.createCrossingMatrix(signed, petal_permutation)
knotDeterminant = spp_library.evaluateKnotDeterminant(crossing_matrix)

print("")
print("You asked for analysis of the knot with petal permutation " + str(petal_permutation))
spp_library.presentDeterminant(knotDeterminant)
print("")
