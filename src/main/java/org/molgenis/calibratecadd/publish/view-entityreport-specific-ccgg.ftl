<#-- modal header -->			
<div class="modal-header">
	<button type="button" class="close" data-dismiss="modal" aria-hidden="true">&times;</button>
	<h4 class="modal-title">DataSet: ${entityMetadata.getLabel()?html}</h4>
</div>

<#-- modal body -->

<div class="modal-body">

	<div class="row">
		<div class="col-md-4">
			<h1>${entity.get('Gene')}</h1><br>
			<#if entity.get('UTestPvalue')??>
				<a href="https://molgenis26.target.rug.nl/ccgg_resources/plots/${entity.get('Gene')}.png" target="_blank"><img src="https://molgenis26.target.rug.nl/ccgg_resources/plots/${entity.get('Gene')}.png" height="238" width="324"/></a>
				<br>
				<a href="https://molgenis26.target.rug.nl/ccgg_resources/data/${entity.get('Gene')}.tsv" target="_blank">download data</a>
			<#else>
				<i>No plot available</i>
			</#if>
		</div>
		<div class="col-md-4">
			<h1>Recommendation</h1><br>
			<h4>${entity.get('Recommendation')?html}</h4>
		</div>
		<div class="col-md-4">
			<div id="modalGenomeBrowser"></div>
		</div>
	</div>

    <div class="row">
	    <div class="col-md-12">
		<#-- Generic entity information split into three columns -->
			<#assign counter = 0 />
			<table class="table">
				<tbody>
					<tr>
						<#list entity.getEntityMetaData().getAtomicAttributes().iterator() as atomicAttribute>
	
	                        <#assign key = atomicAttribute.getName()>
	
							<#if counter == 3>
								</tr>
								<tr>
								<#assign counter = 0>
							</#if>
								
							<th>${key?html}</th>
							<#if entity.get(key)??>
								<#if entity.get(key)?is_sequence>
									<td>
									<#list entity.get(key) as value>
										${value!?html}<#if value_has_next>, </#if>
									</#list>
									</td>
								<#else>
									<td>${entity.getString(key)!?html}</td>
								</#if>
							<#else>
								<td>&nbsp;</td>
							</#if>
							
							<#assign counter = counter + 1>
						</#list>
						
						<#-- fill last row with empty data -->
						<#assign counter = 3 - counter>
						<#list 1..counter as i>
							<th>&nbsp;</th>
							<td>&nbsp;</td>
						</#list>
					</tr>
				</tbody>
			</table>
		</div>
	</div>
</div>

<#-- modal footer -->
<div class="modal-footer">
	<button type="button" class="btn btn-default" data-dismiss="modal">close</button>
</div>

<script>
    molgenis.dataexplorer.data.createGenomeBrowser({
        pageName: 'modalGenomeBrowser', 
        noPersist: true, 
        chr: '${entity.getString("Chr")}', 
        viewStart: ${entity.getString("Start")} - 10000, 
        viewEnd: ${entity.getString("End")} + 10000
    });

    setTimeout(function(){
        $('.modal-body').animate({scrollTop:0},0);
    }, 10);
</script>
