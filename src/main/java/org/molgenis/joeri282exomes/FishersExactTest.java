package org.molgenis.joeri282exomes;

import java.math.BigDecimal;
import java.util.Arrays;

/**
 * 
 * TAKEN FROM https://github.com/bcl2group/GroupData/blob/c051f5408c758c4c18d260921db23f6fe1d2edf2/skanda/code/
 * pipeline_src_compilation_related/pipeline-with-src-april-2012/src/utilities/FishersExactTest.java
 * 
 * @author Kent
 * 
 */
public class FishersExactTest
{

	/**
	 * The scale of the BigDecimal objects used in calculating the hypergeometric distribution. Note that the default
	 * value of 100 should be sufficient for most purposes. If computation is too slow, then reduce the scale;
	 * conversely, if additional accuracy is desired, then increase it.
	 */
	public static int scale = 100;

	/**
	 * Computes the probability of a contingency table with the given values, as given by the hypergeometric
	 * distribution. It is assumed that the values are ordered as follows: <br />
	 * <table border="1">
	 * <tr>
	 * <td>a</td>
	 * <td>b</td>
	 * </tr>
	 * <tr>
	 * <td>c</td>
	 * <td>d</td>
	 * </tr>
	 * </table>
	 * <br />
	 * 
	 * However, any permutation of the rows or columns or reflection across the diagonal is also admissible.
	 * 
	 * @param a
	 *            an entry in the table as described above.
	 * @param b
	 *            an entry in the table as described above.
	 * @param c
	 *            an entry in the table as described above.
	 * @param d
	 *            an entry in the table as described above.
	 * @return the probability that one row occurs given another row.
	 */
	public static double prob(int a, int b, int c, int d)
	{
		BigDecimal result = (new BigDecimal(1)).setScale(scale);
		int numValues = 1;
		for (int i = a + 1; i <= a + b; i++)
		{
			result = result.multiply(new BigDecimal(i));
			result = result.divide(new BigDecimal(numValues++), BigDecimal.ROUND_HALF_EVEN);
		}
		for (int i = b + 1; i <= b + d; i++)
		{
			result = result.multiply(new BigDecimal(i));
			result = result.divide(new BigDecimal(numValues++), BigDecimal.ROUND_HALF_EVEN);
		}
		for (int i = c + 1; i <= c + a; i++)
		{
			result = result.multiply(new BigDecimal(i));
			result = result.divide(new BigDecimal(numValues++), BigDecimal.ROUND_HALF_EVEN);
		}
		for (int i = d + 1; i <= d + c; i++)
		{
			result = result.multiply(new BigDecimal(i));
			result = result.divide(new BigDecimal(numValues++), BigDecimal.ROUND_HALF_EVEN);
		}
		return result.doubleValue();
	}

	/**
	 * Calculates the p-value of the given contingency table using Fisher's exact test. See {@link #FishersExactTest()}
	 * for a description of the input variables.
	 * 
	 * @param a
	 *            an entry in the table as described in {@link #FishersExactTest()}.
	 * @param b
	 *            an entry in the table as described in {@link #FishersExactTest()}.
	 * @param c
	 *            an entry in the table as described in {@link #FishersExactTest()}.
	 * @param d
	 *            an entry in the table as described in {@link #FishersExactTest()}.
	 * @return the p-value.
	 */
	public static double run(int a, int b, int c, int d/* , String alternative */)
	{
		return method2(a, b, c, d);
	}

	private static double method2(int a, int b, int c, int d)
	{
		double oddsRatio = oddsRatio(a, b, c, d);
		double result = 0.0;

		for (int aa = 0; aa <= a + c; ++aa)
		{
			int cc = a + c - aa;
			int bb = a + b - aa;
			int dd = b + d - bb;
			if (bb < 0 || dd < 0 || oddsRatio(aa, bb, cc, dd) < oddsRatio)
			{
				continue;
			}
			result += prob(aa, bb, cc, dd);
		}

		return result;
	}

	private static double oddsRatio(int a, int b, int c, int d)
	{
		return (double) a / b / c * d;
	}

	private static double method1(int a, int b, int c, int d)
	{
		int[] matrix =
		{ a, b, c, d };

		int r1 = a + b;
		int r2 = c + d;
		int c1 = a + c;
		int c2 = b + d;
		int[] sums =
		{ r1, r2, c1, c2 };
		Arrays.sort(sums);
		int min = sums[0];

		if (min == r2)
		{
			swap(matrix, 0, 2);
			swap(matrix, 1, 3);
		}
		else if (min == c1)
		{
			swap(matrix, 1, 2);
		}
		else if (min == c2)
		{
			swap(matrix, 0, 3);
		}

		int minSum = matrix[0] + matrix[1];
		c1 = matrix[0] + matrix[2];
		c2 = matrix[1] + matrix[3];

		double ref = 0.0;
		double result = 0.0;
		double[] pValues = new double[minSum + 1];
		for (int i = 0; i <= minSum; i++)
		{
			int aa = i;
			int bb = minSum - i;
			int cc = c1 - aa;
			int dd = c2 - bb;
			pValues[i] = prob(aa, bb, cc, dd);
			if (aa == matrix[0] && bb == matrix[1])
			{
				ref = pValues[i];
			}
		}

		for (int i = 0; i <= minSum; i++)
		{
			if (pValues[i] <= ref)
			{
				result += pValues[i];
			}
		}
		return result;
	}

	/**
	 * Swaps the positions of the specified elements in the given integer array in place.
	 * 
	 * @param array
	 *            the array of integers.
	 * @param i1
	 *            the index of the first element.
	 * @param i2
	 *            the index of the second element.
	 */
	private static void swap(int[] array, int i1, int i2)
	{
		int temp = array[i1];
		array[i1] = array[i2];
		array[i2] = temp;
	}
}
