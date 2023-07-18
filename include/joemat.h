#ifndef JOEMAT_H
#define JOEMAT_H

void setLieAlgebra(const char* s);

int getLieAlgebraDim();
const char* getLieAlgebraBasis();
const char* getLieAlgebraCentralizer();
const char* getLieAlgebraNormalizer();

// Two options: Either pass in the Lie algebra for each computation
// or set the Lie algebra before the computation
int getLieAlgebraDim(const char* s);
const char* getLieAlgebraBasis(const char* s);
const char* getLieAlgebraCentralizer(const char* s);
const char* getLieAlgebraNormalizer(const char* s);

// int getMatrixRank(const char* s);
// int getZero();

#endif // JOEMAT_H
