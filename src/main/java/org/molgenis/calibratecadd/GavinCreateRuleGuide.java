package org.molgenis.calibratecadd;

import org.molgenis.data.annotation.entity.impl.gavin.GavinEntry;

import java.io.File;
import java.util.HashMap;

/**
 * Created by joeri on 9/22/16.
 */
public class GavinCreateRuleGuide {

    HashMap<String, GavinEntry> gavinData;
    public GavinCreateRuleGuide()
    {
        this.gavinData = loadGAVIN(predictionToolPath + File.separatorChar + "GAVIN_calibrations_"+version+".tsv").getGeneToEntry();
    }


}
