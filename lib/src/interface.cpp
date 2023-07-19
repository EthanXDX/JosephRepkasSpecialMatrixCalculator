#include "headers/interface.h"

const char* lastString = "";
lie_algebra* storedLieAlgebra;
lie_algebra* workingLieAlgebra;

void setLieAlgebra(const char* s) {
    if (strcmp(s, lastString) == 0) return;
    storedLieAlgebra = toLieAlgebra(s);
    workingLieAlgebra = storedLieAlgebra;
}

int getLieAlgebraDim() {
    return workingLieAlgebra->get_dim();
}

const char* getLieAlgebraBasis() {
    vector<g::matrix> basis = workingLieAlgebra->get_basis();
    return toCharArray(basis);
}

const char* getLieAlgebraCentralizer() {
    lie_algebra* centralizer = workingLieAlgebra->compute_centralizer();
    return toCharArray(centralizer);
}

const char* getLieAlgebraNormalizer() {
    lie_algebra* normalizer = workingLieAlgebra->compute_normalizer();
    return toCharArray(normalizer);
}

int getLieAlgebraDim(const char* s) {
    workingLieAlgebra = toLieAlgebra(s);
    int dim = getLieAlgebraDim();
    workingLieAlgebra = storedLieAlgebra;
    return dim;
}

const char* getLieAlgebraBasis(const char* s) {
    workingLieAlgebra = toLieAlgebra(s);
    const char* result = getLieAlgebraBasis();
    workingLieAlgebra = storedLieAlgebra;
    return result;
}

const char* getLieAlgebraCentralizer(const char* s) {
    workingLieAlgebra = toLieAlgebra(s);
    const char* result = getLieAlgebraCentralizer();
    workingLieAlgebra = storedLieAlgebra;
    return result;
}

const char* getLieAlgebraNormalizer(const char* s) {
    workingLieAlgebra = toLieAlgebra(s);
    const char* result = getLieAlgebraNormalizer();
    workingLieAlgebra = storedLieAlgebra;
    return result;
}

// int getMatrixRank(const char* s) {
//     g::matrix x = toMatrix(s);
//     int r = lin_alg::rank(x);
//     return r;
// }
// 
// int getZero() {
//     return 0;
// }
