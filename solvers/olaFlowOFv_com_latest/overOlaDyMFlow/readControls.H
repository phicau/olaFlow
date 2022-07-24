#include "readTimeControls.H"

correctPhi = pimple.dict().lookupOrDefault("correctPhi", false);

checkMeshCourantNo =
    pimple.dict().lookupOrDefault("checkMeshCourantNo", false);

moveMeshOuterCorrectors =
    pimple.dict().lookupOrDefault("moveMeshOuterCorrectors", false);

massFluxInterpolation =
    pimple.dict().lookupOrDefault("massFluxInterpolation", false);

ddtCorr = pimple.dict().lookupOrDefault("ddtCorr", true);

adjustFringe =  pimple.dict().lookupOrDefault("oversetAdjustPhi", false);
