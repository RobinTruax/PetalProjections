#Importing the math library, SymPy (for prime factorization of natural numbers), and Numpy (for determinant calculation)
import math
import numpy
from sympy.ntheory import factorint

#Auxiliary Function 1: r_n(k)
def r(k,n): 
    return k%n

#Auxiliary Function 2: d(k), which depends on the petal number |P|
def d(k, petal_number): 
    if 2*r(k, petal_number-3) < (petal_number-3): 
        return petal_number - 3 - 2*r(k, petal_number-3)
    return petal_number - 5 - 2*r(k, petal_number-3)

#Auxilary Function 3: b(k), which depends on the petal number |P|
def b(k, petal_number): 
    return math.floor(k/(petal_number - 3))

#Auxiliary Function 4: b'(k), which depends on the petal number |P|
def bPrime(k, petal_number): 
    return r(b(k, petal_number)+d(k, petal_number),petal_number)

#Building Block Function 1: Generates the unsigned Gauss Code, depends only on the petal number |P|
def unsignedGaussCode(petal_number): 
    L_array = [i for i in range(int((petal_number-3)/2))] + [int((petal_number-5)/2) - i for i in range(int((petal_number-3)/2))]
    a_array = [0 for i in range(int((petal_number-3)/2))] + [i+1 for i in range(int((petal_number-3)/2))]
    n_array = [r(int((petal_number-1)/2)*math.floor(k/(petal_number-3)) + a_array[r(k,petal_number-3)],petal_number)+1 for k in range(petal_number*(petal_number-3))]

    unsigned_gauss_code = [petal_number*L_array[r(k, petal_number -3)] + n_array[k] for k in range(petal_number*(petal_number-3))]
    return unsigned_gauss_code

#Building Block Function 2: Generates the signed Gauss Code, requires unsigned Gauss code and petal permutation
def signGaussCode(Gauss_code, petal_permutation): 
    petal_number = len(petal_permutation)
    for k in range(petal_number*(petal_number - 3)): 
        if(petal_permutation[b(k, petal_number)] < petal_permutation[bPrime(k, petal_number)]): 
            pass
        else: 
            Gauss_code[k] *= -1
    return Gauss_code

#Building Block Function 3: Splits up a signed Gauss code into an aray of arcs
def splitSignedGaussCode(signed_Gauss_code): 
    #Preparing to split up the signed Gauss code into arcs
    arcs = []
    current_arc = []

    #Creating a list to scan through that starts with the first negative number
    trawl = list(range(len(signed_Gauss_code)))
    for i in range(len(signed_Gauss_code)): 
        if signed_Gauss_code[i] < 0: 
            trawl = trawl[i:] + trawl[:i]
            break

    #Splitting the signed Gauss code up into arcs
    current_arc.append(signed_Gauss_code[trawl[0]])
    for i in trawl[1:]: 
        current_arc.append(signed_Gauss_code[i])
        if signed_Gauss_code[i] < 0: 
            arcs.append(current_arc)
            current_arc = []
            current_arc.append(signed_Gauss_code[i])
    current_arc.append(signed_Gauss_code[trawl[0]])
    arcs.append(current_arc)

    return(arcs)

#Building Block Function 4: Generates the crossing matrix, requires the signed Gauss code and petal permutation
def createCrossingMatrix(signed_Gauss_code, petal_permutation): 
    #Initializing everything
    petal_number = len(petal_permutation)
    matrix_size = math.floor(petal_number*(petal_number-3)/2)
    crossing_matrix = [[0 for i in range(matrix_size)] for i in range(matrix_size)]

    arcs = splitSignedGaussCode(signed_Gauss_code)

    for crossing in range(1, len(crossing_matrix)+1): 
        for arc in range(len(arcs)): 
            if crossing in arcs[arc]: 
                crossing_matrix[crossing-1][arc] = 2
            if -1*crossing in arcs[arc]: 
                crossing_matrix[crossing-1][arc] = -1

    return crossing_matrix

#Building Block Function 5: Evaluates the first minor of the crossing matrix
def evaluateKnotDeterminant(crossing_matrix): 
    #This part slices off the first crossing (first step of getting the minor)
    crossing_matrix = crossing_matrix[1:]

    #This part slices off the first arc in each crossing (second step of getting the minor)
    for crossing in range(len(crossing_matrix)): 
        crossing_matrix[crossing] = crossing_matrix[crossing][1:]
    
    #This part converts our 2D array into a numpy array and then evaluates the determinant
    crossing_numpy_matrix = numpy.array(crossing_matrix)
    return abs(int(round(numpy.linalg.det(crossing_numpy_matrix))))

#This function directly generates the knot determinant of a petal projection with a given petal permutation
def knotDeterminant(petal_permutation): 
    unsigned = unsignedGaussCode(len(petal_permutation))
    signed = signGaussCode(unsigned, petal_permutation)
    crossing_matrix = createCrossingMatrix(signed, petal_permutation)
    return evaluateKnotDeterminant(crossing_matrix)

#This function presents what the knot determinant implies about the possible colorings of the knot
def presentDeterminant(determinant): 
    print("The knot's determinant is " + str(determinant) + ".")
    factorization = factorint(determinant)
    for i in factorization: 
        print("Since " + str(i) + " appears " + str(factorization[i]) + " time(s) in the prime factorization of " + str(determinant) + ", there are " + str(i) + "^(" + str(factorization[i]) + "+1) - " + str(i) + " = " + str(pow(i,factorization[i]+1)-i) + " nontrivial " + str(i) + "-colorings of the knot.")
    print("These are all nontrivial colorings of the knot, though there are p trivial colorings for every prime p.")