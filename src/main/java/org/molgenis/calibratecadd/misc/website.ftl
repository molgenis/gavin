<!doctype html>
<html>
	<head>
		<meta charset="utf-8">
		<script src="jquery.js"></script>
		<script>
		
		$(document).ready(function() {
		
			$('table td:nth-child(2)').each(function() {
				var val = $(this).text();
				
				if (val == 'N1' || val == 'N2' || val == 'N3') { 											$(this).css('backgroundColor', '#ff0000'); }
				else if (val == 'T1' || val == 'T2' || val == 'C3' || val == 'C4' || val == 'C5') { 		$(this).css('backgroundColor', '#ffff00'); }
				else if (val == 'I1' || val == 'I2'  || val == 'I3' ) { 									$(this).css('backgroundColor', '#00ffff'); }
				else if (val == 'C1' || val == 'C2' ) { 													$(this).css('backgroundColor', '#00ff00'); }
			});
			$('table td:nth-child(4)').each(function() {
				var pval = $(this).text();
				
				if ((pval >= 0) && (pval < (0.01))) {		$(this).css('backgroundColor', '#99ff9c'); }
				else if((pval >= 0.01) && (pval < 0.05)) {	$(this).css('backgroundColor', '#fff599'); }
				else if((pval >= 0.05) && (pval <= 1)) {	$(this).css('backgroundColor', '#dddddd'); }
				else {										$(this).css('backgroundColor', '#ffffff'); }
			});
			$('table td').filter(':nth-child(7), :nth-child(8), :nth-child(9), :nth-child(10)').each(function() {
				//use: http://www.bretttolbert.com/projects/colorscale/  max/min hue: 0 - 135, max/min sat: 0.4, max/min val: 1
				var val = $(this).text();
				
				if (val == 'n/a') { 					$(this).css('backgroundColor', '#ffffff'); }
				else if ((val >= 0) && (val < 5)) { 	$(this).css('backgroundColor', '#ffffff'); }
				else if((val >= 5) && (val < 10)) {		$(this).css('backgroundColor', '#99ff9c'); }
				else if((val >= 10) && (val < 15)) {	$(this).css('backgroundColor', '#adff99'); }
				else if((val >= 15) && (val < 20)) {	$(this).css('backgroundColor', '#c4ff99'); }
				else if((val >= 20) && (val < 25)) {	$(this).css('backgroundColor', '#dbff99'); }
				else if((val >= 25) && (val < 30)) {	$(this).css('backgroundColor', '#f2ff99'); }
				else if((val >= 30) && (val < 35)) {	$(this).css('backgroundColor', '#fff599'); }
				else if((val >= 35) && (val < 40)) {	$(this).css('backgroundColor', '#ffde99'); }
				else if((val >= 40) && (val < 45)) {	$(this).css('backgroundColor', '#ffc799'); }
				else if((val >= 45) && (val <= 50)) {	$(this).css('backgroundColor', '#ffb099'); }
				else {									$(this).css('backgroundColor', '#ff9999'); }
			});
		});

		</script>
				
		<style type="text/css">
			table {border: 1px solid black; padding: 0;}
			th {background-color: #3f73c1; color: white; border: 1px solid black; padding: 2px;}
			td {padding: 4px; margin: 0; border: 1px solid black;}
		</style>
		
		<title>Calibrated CADD Gene Guide</title>
</head>

<body>
<div style="margin: 0px auto; display: table; text-align: center; font-family: 'Helvetica', 'Arial', sans-serif;">

<h1>Calibrated CADD Gene Guide - ALPHA 3</h1>
<img src="umcg.jpg" height="50" alt="UMCG" />
<img src="rug.jpg" height="50" alt="RUG" />
<img src="5gpm.png" height="50" alt="5GPM" />
<h5>© 2015, Genomics Coordination Center, Department of Genetics, UMCG Groningen</h5>

<br>

<h2>How to use</h2>
<table style="width: 80%; margin: 0 auto;"><tr><td style="border: 0px solid black; text-align: left;">
	<ol>
		<li>Hit CTRL+F or CMD+F to enable your browser search and find your gene of interest.</li>
		<li>Check the p-value find out if an association between CADD score and pathogenic vs. population variants has been detected:
		<ul>
			<li>If the p-value is <div style="background-color:#99ff9c; display: inline;">green</div>, the results are significant at p-val &lt; 0.01 and usually quite meaningful.</li>
			<li>If the p-value is <div style="background-color:#fff599; display: inline;">yellow</div>, the results are at p-val 0.01 to 0.05, still useful but not very convincing.</li>
			<li>If the p-value is <div style="background-color:#dddddd; display: inline;">gray</div>, there is no proven association, usually due to low sample size.</li>
		</ul>
		</li>
		<li>Use the group means or sensitivity/specificity thresholds as a means to interpret your variants.</li>
		<li>Inspect the plots to get a better feel for the genomic context. Perhaps the CADD scores are more informative for one part of a gene than for another.</li>
		<li>Batch downloads are available for <a href="downloads/clinvar.patho.fix.snpeff.exac.genesumm.tsv">gene summary</a> and <a href="downloads/clinvar.patho.fix.snpeff.exac.withcadd.tsv">variant level</a> data. Analysis source code is available on <a href="https://github.com/molgenis/molgenis-data-cadd" target="_blank">GitHub</a>.</li>
	</ol>
</td></tr></table>

<br>

<h2>Legend</h3>
<table style="font-size: 12px; width: 80%; margin: 0 auto;"><tr><td style="border: 0px solid black; text-align: left;">
	<b>Gene</b> = The HGNC symbol of the gene being tested in this row. Click to view plot. <a href="http://www.genenames.org/" target="_blank">more info</a><br>
	<b>P-value</b> = Mann-Whitney <i>U</i> test p-value. Probability of population variants being a different from pathogenic variants. <a href="https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test" target="_blank">more info</a><br>
	<b>nPatho</b> = Number of pathogenic variants used (ClinVar). <a href="http://www.ncbi.nlm.nih.gov/clinvar/" target="_blank">more info</a><br>
	<b>nPop</b> = Number of matched population variants used (ExAC). <a href="http://exac.broadinstitute.org/" target="_blank">more info</a><br>
	<b>MeanPatho</b> = Mean CADD score of all variants considered pathogenic. These variants are obtained from ClinVar sources.<br>
	<b>MeanPop</b> = Mean CADD score of all population variants. These variants are closely matched to the properties of known pathogenic variants.<br>
	<b>95% sens</b> = Threshold for 95% sensitivity, set at the 5th percentile of pathogenic variants. Use this to filter out false positives ('benign'), because you will only remove a few true positives (pathogenic) below this threshold.<br>
	<b>95% spec</b> = Threshold for 95% specificity, set at the 95th percentile of population variants. Use this to identify true positives (pathogenic), because you will only find a few false positives ('benign') above this threshold.<br>
	<b>Variants</b> = Click link to download the pathogenic and population variants with their CADD scores used in the assessment. <a href="http://cadd.gs.washington.edu/" target="_blank">more info</a><br>
	<b>Category</b> =
	
	<p><i>Not enough data:</i>
	<ul>
		<li><b>N1</b> = Not enough ClinVar variants available for calibration (less than 2).</li>
		<li><b>N2</b> = No ExAC variants available in ClinVar variant interval.</li>
		<li><b>N3</b> = CADD calibration failed, e.g. due to being unable to score variants.</li>
	</ul></p>
	
	
	<p><i>Pathogenic MAF threshold available:</i>
	<ul>
		<li><b>T1</b> = Pathogenic MAF could be established and >0 ExAC variants in ClinVar interval, but no there are no ExAC variants found below this threshold.</li>
		<li><b>T2</b> = Pathogenic MAF established and >0 ExAC variants in ClinVar interval that are below this threshold, but protein impact distributions could not be equalized. There is also no impact seperation point that can be used as a threshold.</li>
	</ul></p>

	<p><i>Pathogenic MAF and 'protein impact' threshold available:</i>
	<ul>
		<li><b>I1</b> = Like T2, but there is a separation point: pathogenic variants have >0 HIGH, while ExAC has MODERATE or less.</li>
		<li><b>I2</b> = Like T2, but there is a separation point: pathogenic variants have >0 MODERATE, while ExAC has LOW or less.</li>
		<li><b>I3</b> = Like T2, but there is a separation point: pathogenic variants have >0 LOW, while ExAC only has MODIFIER.</li>
	</ul></p>

	<p><i>Results of CADD calibration on pathogenic MAF-filtered, 'protein impact equalized' variants sets. Pathogenic MAF is available, CADD scores may not be informative.</i>
	<ul>
		<li><b>C1</b> = CADD scores are highly informative for this gene (patho > popul , P < 0.01).</li>
		<li><b>C2</b> = CADD scores are somewhat informative for this gene (patho > popul , 0.01 ≤ P ≤ 0.05).</li>
		<li><b>C3</b> = We don't know if CADD scores are informative. (P > 0.05 and < 5 samples in one or both groups)</li>
		<li><b>C4</b> = Enough samples (> 5 in each group), but CADD scores are not informative for this gene (P > 0.05)</li>
		<li><b>C5</b> = Artifact, e.g. P ≤ 0.05 but patho mean greather than population mean.</li>
	</ul></p>
</td></tr></table>

<br><br>


<table style="text-align: center; width: 80%; margin: 0 auto;">
<tr style="text-align: center;"><th>Gene</th><th>Category</th><th>Info</th><th>P-value</th><th>nPatho</th><th>nPop</th><th>MeanPatho<th>MeanPop</th><th>95% sens</th><th>95% spec</th><th>Variants</th></tr>

<#list genes as gene>
<#--tr><td><a href="plots/BRCA2.png" target="_blank">BRCA2</a></td><td>6.45630831985457e-68</td><td>2964</td><td>0.7735326157</td><td>16.12</td><td>0.856502242152466</td><td>0.911347517730497</td><td>0.952518262206844</td><td>289</td><td>216</td><td><a href="data/BRCA2.tsv">download</a></td></tr-->
<tr><td><#if gene[8] == 'n/a'>${gene[0]}<#else><a href="plots/${gene[0]}.png" target="_blank">${gene[0]}</a></#if></td><td>${gene[1]}</td><td>todo</td><td><#if gene[8] == 'n/a'>n/a<#else>${gene[8]?number?string["0.####"]}</#if></td><td>${gene[3]}</td><td>${gene[4]}</td><td>${gene[5]}</td><td>${gene[6]}</td><td>${gene[9]}</td><td>${gene[10]}</td><td><#if gene[8] == 'n/a'>n/a<#else><a href="data/${gene[0]}.tsv">download</a></#if></td></tr>
</#list>
</table>

</div>

</body>
</html>