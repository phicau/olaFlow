bool correctPhi
(
    pimple.dict().lookupOrDefault("correctPhi", true)
);

bool checkMeshCourantNo
(
    pimple.dict().lookupOrDefault("checkMeshCourantNo", false)
);

bool moveMeshOuterCorrectors
(
    pimple.dict().lookupOrDefault("moveMeshOuterCorrectors", false)
);


bool massFluxInterpolation
(
    pimple.dict().lookupOrDefault("massFluxInterpolation", false)
);

bool adjustFringe
(
    pimple.dict().lookupOrDefault("oversetAdjustPhi", false)
);

bool ddtCorr
(
    pimple.dict().lookupOrDefault("ddtCorr", true)
);
