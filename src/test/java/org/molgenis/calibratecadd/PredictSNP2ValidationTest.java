package org.molgenis.calibratecadd;

import org.testng.annotations.Test;

/**
 * Created by joeri on 9/1/16.
 */
public class PredictSNP2ValidationTest {

    @Test
    public void testPredictionTool() throws Exception {
        new AbstractValidationTest("SIFT").testPredictionTool();
    }

}
