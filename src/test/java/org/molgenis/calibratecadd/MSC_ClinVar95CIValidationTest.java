package org.molgenis.calibratecadd;

import org.testng.annotations.Test;

/**
 * Created by joeri on 9/1/16.
 */
public class MSC_ClinVar95CIValidationTest {

    @Test
    public void testPredictionTool() throws Exception {
        new AbstractValidationTest("MSC_ClinVar95CI").testPredictionTool();
    }

}
