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

    @Override
    public String toString() {
        // Gene	Category	Chr	Start	End	NrOfPopulationVariants	NrOfPathogenicVariants	NrOfOverlappingVariants	NrOfFilteredPopVariants	PathoMAFThreshold	PopImpactHighPerc   PopImpactModeratePerc	PopImpactLowPerc	PopImpactModifierPerc	PathoImpactHighPerc	PathoImpactModeratePerc	PathoImpactLowPerc	PathoImpactModifierPerc PopImpactHighEq	PopImpactModerateEq	PopImpactLowEq	PopImpactModifierEq


        StringBuffer sb = new StringBuffer();
        sb.append(this.ge.gene + "\t");
        sb.append(this.ge.category + "\t");
        sb.append(this.ge.chromosome + "\t");
        sb.append(this.ge.start + "\t");
        sb.append(this.ge.end + "\t");
        sb.append((this.ge.NrOfPopulationVariants == null ? "" : this.ge.NrOfPopulationVariants) + "\t");
        sb.append((this.ge.NrOfPathogenicVariants == null ? "" : this.ge.NrOfPathogenicVariants) + "\t");
        sb.append((this.ge.NrOfOverlappingVariants == null ? "" : this.ge.NrOfOverlappingVariants) + "\t");
        sb.append((this.ge.NrOfFilteredPopVariants == null ? "" : this.ge.NrOfFilteredPopVariants) + "\t");
        sb.append((this.ge.PathoMAFThreshold == null ? "" : this.ge.PathoMAFThreshold) + "\t");
        sb.append((this.ge.PopImpactHighPerc == null ? "" : this.ge.PopImpactHighPerc) + "\t");
        sb.append((this.ge.PopImpactModeratePerc == null ? "" : this.ge.PopImpactModeratePerc) + "\t");
        sb.append((this.ge.PopImpactLowPerc == null ? "" : this.ge.PopImpactLowPerc) + "\t");
        sb.append((this.ge.PopImpactModifierPerc == null ? "" : this.ge.PopImpactModifierPerc)+ "\t");
        sb.append((this.ge.PathoImpactHighPerc == null ? "" : this.ge.PathoImpactHighPerc) + "\t");
        sb.append((this.ge.PathoImpactModeratePerc == null ? "" : this.ge.PathoImpactModeratePerc)+ "\t");
        sb.append((this.ge.PathoImpactLowPerc == null ? "" : this.ge.PathoImpactLowPerc) + "\t");
        sb.append((this.ge.PathoImpactModifierPerc == null ? "" : this.ge.PathoImpactModifierPerc) + "\t");
        sb.append((this.ge.PopImpactHighEq == null ? "" : this.ge.PopImpactHighEq) + "\t");
        sb.append((this.ge.PopImpactModerateEq == null ? "" : this.ge.PopImpactModerateEq)+ "\t");
        sb.append((this.ge.PopImpactLowEq == null ? "" : this.ge.PopImpactLowEq)+ "\t");
        sb.append((this.ge.PopImpactModifierEq == null ? "" : this.ge.PopImpactModifierEq));

        return sb.toString();
    }
}
