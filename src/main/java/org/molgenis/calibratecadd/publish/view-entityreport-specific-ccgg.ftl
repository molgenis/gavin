<#-- modal header -->           
<div class="modal-header">
    <button type="button" class="close" data-dismiss="modal" aria-hidden="true">&times;</button>
    <h4 class="modal-title"><#-- Custom entity report title goes here --></h4>
</div>

<#-- modal body -->
<div class="modal-body">
    <#-- Custom entity report contents go here -->

        chr: ${entity.getString("Chrom")}, 
        viewStart: ${entity.getString("Start")}
        viewEnd: ${entity.getString("Stop")}
</div>


<div id="modalGenomeBrowser"></div>



<#-- modal footer -->
<div class="modal-footer">
    <button type="button" class="btn btn-default" data-dismiss="modal">close</button>
</div>

