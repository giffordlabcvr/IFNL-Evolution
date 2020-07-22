// List of genes to scan for
var alignments = ["IFNLb1", "IFNLb2" ];

// Get list of variations to scan for 
// list variation -w "name like 'n-linked-glycosylation%'" 
var variationTable;



glue.logInfo("Processing sequence ID: "+id);

glue.inMode("/alignment/"+alignment+"/"+alignmentName, function() {

	var rowObjs = glue.tableToObjects(glue.command(["variation","member","scan","-r","REF_IFNLb2_Mammal_MASTER","-f","glycosylation-region-1", "-v", "n-linked-glycosylation-o1"]));
	for(var i = 0; i < rowObjs.length; i++) {
  
		var alignmentName       = rowObjs[i].alignmentName;
		var sequenceID          = rowObjs[i].sequenceID;
		var sufficientCoverage  = rowObjs[i].sufficientCoverage;
		var present             = rowObjs[i].present;

		glue.logInfo("Processing sequence ID: "+sequenceID+"COVERAGE: "+sufficientCoverage+", PRESENT: "+present);
	
	}		
});
