function norm = getLieAlgebraNormalizer(generators)

generatorsString = matrixSeqToString(generators);
clib.joemat.setLieAlgebra(generatorsString);
normString = clib.joemat.getLieAlgebraNormalizer();
norm = stringToMatrixSeq(normString);
