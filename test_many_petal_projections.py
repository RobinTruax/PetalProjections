import spp_library
import timeit

petal_projections_to_test = [[1, 3, 5, 2, 4],[1, 3, 5, 2, 7, 4, 6],[1, 3, 5, 2, 8, 4, 6, 9, 7]]
number_of_runs_each = 100

def mean(runtimes): 
    return sum(runtimes) / len(runtimes) 

for petal_projection in petal_projections_to_test: 
    #Compute the knot determinant
    petal_projection_determinant = spp_library.knotDeterminant(petal_projection)

    alg_runtimes = []
    #Run the petal projection calculation number_of_runs_each times and get a runtime list
    for i in range(number_of_runs_each): 
        start = timeit.default_timer()
        spp_library.knotDeterminant(petal_projection)
        alg_runtimes.append(timeit.default_timer() - start)

    #Print everything
    print("Petal Projection: " + str(petal_projection) + " | Determinant: " + str(petal_projection_determinant) + " | Minimum Runtime: " + str(min(alg_runtimes)) + " | Average Runtime: " + str(mean(alg_runtimes)))