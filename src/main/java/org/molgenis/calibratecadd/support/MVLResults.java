package org.molgenis.calibratecadd.support;

/*
 * Results per MVL, e.g. CAR, DYS, etc.
 */
public class MVLResults {

	//in the original MVL
	public int nrOf_B_LB;
	public int nrOf_VOUS;
	public int nrOf_P_LP;
	
	//classified via CCGG
	public int nrOf_B_LB_judgedAs_B_LB;
	public int nrOf_B_LB_judgedAs_P_LP;
	public int nrOf_VOUS_judgedAs_B_LB;
	public int nrOf_VOUS_judgedAs_P_LP;
	public int nrOf_P_LP_judgedAs_B_LB;
	public int nrOf_P_LP_judgedAs_P_LP;
	
	public MVLResults()
	{
		nrOf_B_LB = 0;
		nrOf_VOUS = 0;
		nrOf_P_LP = 0;
		
		nrOf_B_LB_judgedAs_B_LB = 0;
		nrOf_B_LB_judgedAs_P_LP = 0;
		
		nrOf_VOUS_judgedAs_B_LB = 0;
		nrOf_VOUS_judgedAs_P_LP = 0;
		
		nrOf_P_LP_judgedAs_B_LB = 0;
		nrOf_P_LP_judgedAs_P_LP = 0;
		
	}

	@Override
	public String toString() {
		return "MVLResults [nrOf_B_LB=" + nrOf_B_LB + ", nrOf_VOUS=" + nrOf_VOUS + ", nrOf_P_LP=" + nrOf_P_LP
				+ ", nrOf_B_LB_judgedAs_B_LB=" + nrOf_B_LB_judgedAs_B_LB + ", nrOf_B_LB_judgedAs_P_LP="
				+ nrOf_B_LB_judgedAs_P_LP + ", nrOf_VOUS_judgedAs_B_LB=" + nrOf_VOUS_judgedAs_B_LB
				+ ", nrOf_VOUS_judgedAs_P_LP=" + nrOf_VOUS_judgedAs_P_LP + ", nrOf_P_LP_judgedAs_B_LB="
				+ nrOf_P_LP_judgedAs_B_LB + ", nrOf_P_LP_judgedAs_P_LP=" + nrOf_P_LP_judgedAs_P_LP + "]";
	}
	
	/**
	 * Proportion of benign/LB variants classified by CCGG as pathogenic/LP, so we make a FP type of mistake according to the human expert.
	 * @return
	 */
	public double getPercOfFalsePositives()
	{
		return ((double)nrOf_B_LB_judgedAs_P_LP/(double)nrOf_B_LB_judgedAs_B_LB)*100;
	}
	
	/**
	 * Proportion of pathogenic/LP variants classified by CCGG as benign/LB, so we make a FN type of mistake according to the human expert.
	 * @return
	 */
	public double getPercOfFalseNegatives()
	{
		return ((double)nrOf_P_LP_judgedAs_B_LB/(double)nrOf_P_LP)*100;
	}
}
