<!doctype html>
<html>
	<head>
		<meta charset="utf-8">
		<script src="jquery.js"></script>
		<script>
		
		$(document).ready(function() {
			$('table td:nth-child(2)').each(function() {
				var pval = $(this).text();
				if ((pval >= 0) && (pval < (0.01))) {		$(this).css('backgroundColor', '#99ff9c'); }
				else if((pval >= 0.01) && (pval < 0.05)) {	$(this).css('backgroundColor', '#fff599'); }
				else if((pval >= 0.05) && (pval <= 1)) {	$(this).css('backgroundColor', '#dddddd'); }
				else {										$(this).css('backgroundColor', '#ffffff'); }
			});
			$('table td').filter(':nth-child(5), :nth-child(6), :nth-child(7), :nth-child(8)').each(function() {
				//use: http://www.bretttolbert.com/projects/colorscale/  max/min hue: 0 - 135, max/min sat: 0.4, max/min val: 1
				var val = $(this).text();
				
				if ((val >= 0) && (val < 5)) { 			$(this).css('backgroundColor', '#ffffff'); }
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
			event.preventDefault();
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
<div style="margin: 0px auto; display: table; text-align: center;">

<h1>Calibrated CADD Gene Guide - BETA 1</h1>
<img src="umcg.jpg" height="50" alt="UMCG" />
<img src="rug.jpg" height="50" alt="RUG" />
<img src="5gpm.png" height="50" alt="5GPM" />
<h5>© 2015, Genomics Coordination Center, Department of Genetics, UMCG Groningen</h5>

<br>

<h2>How to use:</h2>
<table style="width: 80%; margin: 0 auto;"><tr><td style="border: 0px solid black; text-align: left;">
	<ol>
		<li>Hit CTRL+F or CMD+F to enable your browser search and find your gene of interest.</li>
		<li>Check the U-test p-value find out if an association between CADD score and pathogenic vs. population variants has been detected:
		<ul>
			<li>If the p-value is <div style="background-color:#99ff9c; display: inline;">green</div>, the results are significant at p-val <0.01 and usually quite meaningful.</li>
			<li>If the p-value is <div style="background-color:#fff599; display: inline;">yellow</div>, the results are at p-val <0.05, useful but not very convincing.</li>
			<li>If the p-value is <div style="background-color:#dddddd; display: inline;">gray</div>, there is no proven association, usually due to low sample size.</li>
		</ul>
		</li>
		<li>Open the plot by clicking on the gene name. Inspect the optimal group separation made by Youden's cutpoint.</li>
		<li>A high Youden's index make the cutpoint more reliable. Other properties (e.g. PPV and NPV at cutpoint) may also help you decide to trust it.</li>
		<li>You can always eyeball the plot to find a different (e.g. more stringent) cutpoint that you trust more.</li>
	</ol>
</td></tr></table>

<br>

<h3>Legend</h3>
<table style="font-size: 14px; width: 80%; margin: 0 auto;"><tr><td style="border: 0px solid black; text-align: left;">
	<b>Gene</b> = The HGNC symbol of the gene being tested in this row. Click to view plot. <a href="http://www.genenames.org/">more info</a><br>
	<b>P-value</b> = Mann-Whitney <i>U</i> test p-value. Probability of population variants being a different from pathogenic variants. <a href="https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test">more info</a><br>
	<b>nPatho</b> = Number of pathogenic variants used (ClinVar). <a href="http://www.ncbi.nlm.nih.gov/clinvar/">more info</a><br>
	<b>nPop</b> = Number of matched population variants used (ExAC). <a href="http://exac.broadinstitute.org/">more info</a><br>
	<b>MeanPatho</b> = xx<br>
	<b>MeanPop</b> = xx<br>
	<b>95% sens</b> = xx<br>
	<b>95% spec</b> = xx<br>
	<b>Variants</b> = Click link to download the pathogenic and population variants with their CADD scores used in the assessment. <a href="http://cadd.gs.washington.edu/">more info</a><br>
</td></tr></table>

<br><br>

<table style="text-align: left;">
<tr style="text-align: center;"><th>Gene</th><th>P-value</th><th>nPatho</th><th>nPop</th><th>MeanPatho><th>MeanPop</th><th>95% sens</th><th>95% spec</th><th>Variants</th></tr>

<#list genes as gene>
<#--tr><td><a href="plots/BRCA2.png" target="_blank">BRCA2</a></td><td>6.45630831985457e-68</td><td>2964</td><td>0.7735326157</td><td>16.12</td><td>0.856502242152466</td><td>0.911347517730497</td><td>0.952518262206844</td><td>289</td><td>216</td><td><a href="data/BRCA2.tsv">download</a></td></tr-->
<tr><td><a href="plots/${gene[0]}.png" target="_blank">${gene[0]}</td><td>${gene[6]}</td><td>${gene[1]}</td><td>${gene[2]}</td><td>${gene[3]}</td><td>${gene[4]}</td><td>${gene[7]}</td><td>${gene[8]}</td><td><a href="data/${gene[0]}.tsv">download</a></td></tr>
</#list>
</table>

</div>

</body>
</html>