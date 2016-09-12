package org.molgenis.calibratecadd.structs;

import org.molgenis.calibratecadd.support.EntityPlus;
import org.molgenis.data.annotation.entity.impl.gavin.GavinEntry;

import java.util.HashMap;
import java.util.List;

/**
 * Created by joeri on 9/12/16.
 */
public class GeneCalibResult {

    public GavinEntry ge;
    public List<EntityPlus> matchedVariants;

    public GeneCalibResult(GavinEntry ge, List<EntityPlus> matchedVariants) {
        this.ge = ge;
        this.matchedVariants = matchedVariants;
    }
}
