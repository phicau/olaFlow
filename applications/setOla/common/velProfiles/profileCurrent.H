    // Dummy field, current introduced later

    patchU[cellIndex] = 0 + (1.0 - alpha[cellIndex])*uWind_.component(0);;
    patchV[cellIndex] = 0 + (1.0 - alpha[cellIndex])*uWind_.component(1);;
    patchW[cellIndex] = 0 + (1.0 - alpha[cellIndex])*uWind_.component(2);; 
