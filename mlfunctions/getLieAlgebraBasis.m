function basis = getLieAlgebraDim(generators)

generatorsString = matrixSeqToString(generators);
clib.joemat.setLieAlgebra(generatorsString);
basisString = clib.joemat.getLieAlgebraBasis();
basis = stringToMatrixSeq(basisString);
