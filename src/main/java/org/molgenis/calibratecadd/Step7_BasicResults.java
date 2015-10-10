package org.molgenis.calibratecadd;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Scanner;

/**
 * 
 * Read:
 * 
 * gene	chr	pos	ref	alt	group	cadd
 * IFT172	2	27680627	A	T	PATHOGENIC	28.0
 * IFT172	2	27700177	A	T	PATHOGENIC	15.99
 * IFT172	2	27693963	C	T	PATHOGENIC	36
 * 
 * Write:
 * 
 * gene	nPath	nPopul	medianPatho	medianPopul	medianDiff
 * IFT172	10	14	22.35	25.69	3.34
 *
 *(more? highest, lowest, averages...)
 *
 */
public class Step7_BasicResults {

	public static void main(String[] args) throws FileNotFoundException {
		Scanner res = new Scanner(new File(args[0]));
		PrintWriter out = new PrintWriter(new File(args[1]));
		
		res.nextLine();
		String line = null;
		while(res.hasNextLine())
		{
			line = res.nextLine();
			
		}

	}

}
