function dim = getLieAlgebraDim(generators)

generatorsString = matrixSeqToString(generators);
clib.joemat.setLieAlgebra(generatorsString);
dim = clib.joemat.getLieAlgebraDim();
