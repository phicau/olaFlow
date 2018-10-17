bool correctPhi
(
    pimple.dict().lookupOrDefault<Switch>("correctPhi", true)
);

bool checkMeshCourantNo
(
    pimple.dict().lookupOrDefault<Switch>("checkMeshCourantNo", false)
);

bool moveMeshOuterCorrectors
(
    pimple.dict().lookupOrDefault<Switch>("moveMeshOuterCorrectors", false)
);


bool massFluxInterpolation
(
    pimple.dict().lookupOrDefault("massFluxInterpolation", false)
);

bool ddtCorr
(
    pimple.dict().lookupOrDefault("ddtCorr", true)
);
