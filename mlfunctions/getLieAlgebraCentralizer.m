function cent = getLieAlgebraCentralizer(generators)

generatorsString = matrixSeqToString(generators);
clib.joemat.setLieAlgebra(generatorsString);
centString = clib.joemat.getLieAlgebraCentralizer();
cent = stringToMatrixSeq(centString);
