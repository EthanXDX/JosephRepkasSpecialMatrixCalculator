#include <stdio.h>
#include "../src/interface.h"

using std::cout;
using std::endl;

void setUp(const char* generatorsString) {
    setLieAlgebra(generatorsString);
}

void testDim() {
    int dim = getLieAlgebraDim();
    cout << "Dimension: " << dim << endl;
}

void testBasis() {
    const char* basisString = getLieAlgebraBasis();
    cout << "Basis: " << endl << basisString << endl;
}

void testNormalizer() {
    const char* normalizerString = getLieAlgebraNormalizer();
    cout << "Normalizer:" << endl << normalizerString << endl;
}

void testCentralizer() {
    const char* centralizerString = getLieAlgebraCentralizer();
    cout << "Centralizer:" << endl << centralizerString << endl;
}

void testLieAlgebra2LieAlgebra(lie_algebra* l) {
    string strl = toString(l);
    lie_algebra* k = toLieAlgebra(strl);
    string strk = toString(k);
    
    if (!k->equals(l)) {
        cout << "L generators:" << endl << strl << endl;
        cout << "K generators:" << endl << strk << endl;
        throw std::runtime_error("testLieAlgebra2LieAlgebra failed");
    }
}

void testStrToMatrixSameSymbols(string str) {
    g::matrix a = toMatrix(str);
    g::matrix b = toMatrix(str);
    if (!(a - b).is_zero_matrix()) throw std::runtime_error("testStrToMatrixSameSymbols failed");
}

int main(int argc, char* argv[]) {
    // const char* s = "[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 0]";
    // int rank = getMatrixRank(s);
    // cout << "Matrix: " << s << endl;
    // cout << "Rank: " << rank << endl;

    // const char* generatorsString = "[1 0;0 0]\n[-1 0;0 0]";
    // int dim = getLieAlgebraDim(generatorsString);
    // const char* basisString = getLieAlgebraBasis(generatorsString);

    // Working
    const char* L5_1 = "[0.25 1 0 0;0 0.25 0 0;0 0 0.25 0;0 0 0 -0.75]\n[0 1 0 0;0 0 1 0;0 0 0 0;0 0 0 0]\n[0 0 1 0;0 0 0 0;0 0 0 0;0 0 0 0]";
    const char* L5_2 = "[0.25 0 0 0;0 0.25 0 1;0 0 -0.75 0;0 0 0 0.25]\n[0 1 0 0;0 0 0 0;0 0 0 0;0 0 0 0]\n[0 0 0 1;0 0 0 0;0 0 0 0;0 0 0 0]";
    const char* L5_4 = "[0 1 0 0;0 0 1 0;0 0 0 1;0 0 0 0]\n[0 0 a 0;0 0 0 a+1;0 0 0 0;0 0 0 0]\n[0 0 0 1;0 0 0 0;0 0 0 0;0 0 0 0]";
    const char* L5_5 = "[0 1 0 0;0 0 0 0;0 0 0 1;0 0 0 0]\n[0 0 a 0;0 0 0 a+1;0 0 0 0;0 0 0 0]\n[0 0 0 1;0 0 0 0;0 0 0 0;0 0 0 0]";
    // Working but slow
    const char* L3_0_6 = "[0 0 0 0 0 0 2 0 0 0;0 2 0 0 0 0 0 0 0 0;0 0 -2 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 1 0 0 0 1 0;0 0 0 0 0 -1 0 0 0 1;0 0 0 1 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 1 0;0 0 0 0 0 0 0 0 0 -1]\n[-0.5 0 0 0 0 0 0 0 0 0;0 -0.5 0 0 0 0 0 0 0 0;0 0 1.5 0 0 0 0 0 0 0;0 0 0 -0.5 0 0 0 0 0 0;0 0 0 0 -0.5 0 0 0 0 0;0 0 0 0 0 0.5 0 0 0 0;0 0 0 0 0 0 -0.5 0 0 0;0 0 0 0 0 0 0 0.5 0 0;0 0 0 0 0 0 0 0 -0.5 0;0 0 0 0 0 0 0 0 0 0.5]\n[0 0 0 0 2 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0;0 1 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 1 0 0;0 0 0 0 0 0 0 0 1 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0]";
    const char* L3_0_9a = "[0 0 0 0 0 0 2 0 0 0;0 2 0 0 0 0 0 0 0 0;0 0 -2 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 1 0 0 0 1 0;0 0 0 0 0 -1 0 0 0 1;0 0 0 1 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 1 0;0 0 0 0 0 0 0 0 0 -1]\n[2*a 0 0 0 0 0 0 0 0 0;0 2*a+2 0 0 0 0 0 0 0 0;0 0 -6*a-2 0 0 0 0 0 0 0;0 0 0 2*a 0 0 0 0 0 0;0 0 0 0 2*a+1 0 0 0 0 0;0 0 0 0 0 -2*a-1 0 0 0 0;0 0 0 0 0 0 2*a 0 0 0;0 0 0 0 0 0 0 -2*a 0 0;0 0 0 0 0 0 0 0 2*a+1 0;0 0 0 0 0 0 0 0 0 -2*a-1]\n[0 0 0 0 2 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0;0 1 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 1 0 0;0 0 0 0 0 0 0 0 1 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0]";
    const char* L_3__0_25_5 = "[-1 0 0 0 0 2 0 0 0 0;0 1 0 0 0 0 0 0 0 0;0 0 -1 0 0 0 0 0 0 0;0 0 0 1 0 0 0 0 0 0;0 0 0 0 0 0 0 1 0 0;0 0 1 0 0 -1 0 0 0 0;0 0 0 0 0 0 0 0 0 1;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 1 0;0 0 0 0 0 0 0 0 0 0]\n[0 0 0 0 2 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 2 0 0 0;0 0 0 0 0 0 0 0 0 0;0 1 0 0 0 0 0 0 0 0;0 0 0 0 0 0 1 1 0 0;0 0 0 0 0 0 0 0 1 0;0 0 0 0 0 0 0 0 1 0;0 0 0 0 0 0 0 0 0 0;0 0 0 1 0 0 0 0 0 0]\n[0 0 0 0 0 0 1 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 1 0;0 0 0 0 0 0 0 0 0 1;0 0 0 1 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0]";
    // Not working
    const char* L_5_6a = "[0 0 0 0 0 0;0 0 1 1 0 0;0 0 0 0 1 0;0 0 0 0 1 0;0 0 0 0 0 0;0 0 0 0 0 0]\n[0 0 a+1 -a 0 0;0 0 0 1 0 0;0 0 0 0 1 a;0 0 0 0 0 -a-1;0 0 0 0 0 0;0 0 0 0 0 0]\n[0 0 0 0 -1 0;0 0 0 0 0 -1;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0]";
    const char* L_3_0_10a = "[0 0 0 0 0 0 2 0 0 0;0 2 0 0 0 0 0 0 0 0;0 0 -2 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 1 0 0 0 1 0;0 0 0 0 0 -1 0 0 0 1;0 0 0 1 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 1 0;0 0 0 0 0 0 0 0 0 -1]\n[2*a 0 0 0 0 0 0 0 0 0;0 -6*a-2 0 0 0 0 0 0 0 0;0 0 2*a+2 0 0 0 0 0 0 0;0 0 0 2*a 0 0 0 0 0 0;0 0 0 0 -2*a-1 0 0 0 0 0;0 0 0 0 0 2*a+1 0 0 0 0;0 0 0 0 0 0 2*a 0 0 0;0 0 0 0 0 0 0 -2*a 0 0;0 0 0 0 0 0 0 0 -2*a-1 0;0 0 0 0 0 0 0 0 0 2*a+1]\n[0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 2;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 1 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 1 0;0 0 0 0 0 0 0 0 0 0;0 0 0 1 0 0 0 0 0 0]";
    const char* L_3_0_11a = "[0 0 0 0 0 0 2 0 0 0;0 2 0 0 0 0 0 0 0 0;0 0 -2 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 1 0 0 0 1 0;0 0 0 0 0 -1 0 0 0 1;0 0 0 1 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 1 0;0 0 0 0 0 0 0 0 0 -1]\n[2*a 0 0 0 0 0 0 0 0 0;0 -2*a-2 0 0 0 0 0 0 0 0;0 0 -2*a+2 0 0 0 0 0 0 0;0 0 0 2*a 0 0 0 0 0 0;0 0 0 0 -1 0 0 0 0 0;0 0 0 0 0 1 0 0 0 0;0 0 0 0 0 0 2*a 0 0 0;0 0 0 0 0 0 0 -2*a 0 0;0 0 0 0 0 0 0 0 -1 0;0 0 0 0 0 0 0 0 0 1]\n[0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 2 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 1 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0;0 1 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 1 0]";
    const char* L_3__1_4_6 = "[-1 0 0 0 0 0 0 0 0 0;0 1 0 0 0 0 0 0 2 0;0 0 -1 0 0 0 0 0 0 0;0 0 0 1 0 0 0 0 0 0;0 0 0 0 0 0 1 0 0 0;0 0 0 0 0 -1 0 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 1;0 0 0 1 0 0 0 0 1 0;0 0 0 0 0 0 0 0 0 0]\n[0 0 0 0 2 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 2 0 0 0;0 0 0 0 0 0 0 0 0 0;0 1 0 0 0 0 0 0 0 0;0 0 0 0 0 0 1 1 0 0;0 0 0 0 0 0 0 0 1 0;0 0 0 0 0 0 0 0 1 0;0 0 0 0 0 0 0 0 0 0;0 0 0 1 0 0 0 0 0 0]\n[0 0 0 0 0 0 1 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 1 0;0 0 0 0 0 0 0 0 0 1;0 0 0 1 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0]";

    const char* generatorsString = L_3__0_25_5;
    cout << "Generators:" << endl << generatorsString << endl;
    setUp(generatorsString);

    testDim();
    testBasis();
    testCentralizer();
    testNormalizer();

    lie_algebra* l = toLieAlgebra(generatorsString);
    testLieAlgebra2LieAlgebra(l);

    string matrixString = "[x y z;a b c;0 1 0]";
    testStrToMatrixSameSymbols(matrixString);

    string L5_1_string(L5_1);
    string L5_2_string(L5_2);
    string matrixInput = L5_1_string + "\n\n\n@\n" + L5_2_string;
    vector<lie_algebra*> algebraSeq = toLieAlgebraSequence(matrixInput);
    cout << toString(algebraSeq[0]) << endl;
    cout << toString(algebraSeq[1]) << endl;

    return 0;
}
