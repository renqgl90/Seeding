TCut Signal      =  TCut("DeltaXBestCombinationT1>-999 && (ISLONG || isDown) ");
TCut over5       =  TCut("DeltaXBestCombinationT1>-999 && Track_P>5000");
TCut  over2less5  =  TCut("DeltaXBestCombinationT1>-999 && Track_P>2000 && Track_P<5000");
TCut  less5       =  TCut("DeltaXBestCombinationT1>-999 && Track_P<5000");
TCut  etain25     =  TCut("DeltaXBestCombinationT1>-999 && Track_eta>2 && Track_eta<5");
