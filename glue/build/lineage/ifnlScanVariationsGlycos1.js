
// Load EVE reference data from tab file 
var scanResult;
glue.inMode("alignment/AL_IFNL_MAMMAL_B", function() {
	scanResult = glue.tableToObjects(glue.command(["variation", "member", "scan", "-r", "REF_IFNL_Mammal_b_MASTER", "-f", "glycosylation-region-1", "-v", "n-linked-glycosylation-o1", "-t"]));
	//glue.log("INFO", "load result was:", scanResult);
});

// Iterate through results
 _.each(scanResult, function(resultObj) {

	 var sequenceID = resultObj["sequenceID"];
	 var queryAAs   = resultObj["queryAAs"];
	 
	 //glue.log("INFO", "Reference name result was:", referenceName);
	 //glue.log("INFO", "Query AAs result was:", queryAAs);
	 
    // Update the variation field
    glue.inMode("sequence/ncbi-curated-mammalia/"+sequenceID, function() {
    
		glue.command(["set", "field", "glycosylation_o1_present", queryAAs]);       
    });

});
