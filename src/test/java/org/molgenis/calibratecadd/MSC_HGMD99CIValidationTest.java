package org.molgenis.calibratecadd;

import org.testng.annotations.Test;

/**
 * Created by joeri on 9/1/16.
 */
public class MSC_HGMD99CIValidationTest {

    @Test
    public void testPredictionTool() throws Exception {
        new AbstractValidationTest("MSC_HGMD99CI").testPredictionTool();
    }

}
