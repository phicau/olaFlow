dictionary tmpDict;

tmpDict.add("theory", theory);

tmpDict.add("nPaddles", nPaddles);
tmpDict.add("initialWaterDepths", initialWaterDepths);
tmpDict.add("absDir", absDir);
if (wavePeriod > 0)
{
    tmpDict.add("wavePeriod", wavePeriod);
}
tmpDict.add("uCurrent", uCurrent);
tmpDict.add("nEdgeMin", nEdgeMin);
tmpDict.add("nEdgeMax", nEdgeMax);

tmpDict.add("meanAngles", meanAngles);
tmpDict.add("zSpanL", zSpanL);

AWADict_ = tmpDict;
