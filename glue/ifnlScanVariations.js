
// Get list of variations to scan for 
var variationTable = {};
get_variations(variationTable);

// Scan alignments for variations
var variationResultTable = {};
scan_for_variations(variationTable)

// Enter results in table
scan_for_variations(variationTable)



// Save table result to file



// Get feature names for a give reference sequence
function get_variations(variationTable) {

    // list variation -w "name like 'n-linked-glycosylation%'" 
	var rowObjs = glue.tableToObjects(glue.command(["list","variation","-w","name like 'n-linked-glycosylation%'"]));
	
	for(var i = 0; i < rowObjs.length; i++) {
  
        var rowObj = rowObjs[i]
		var variationName = rowObj.name;
		var featureName   = rowObj["featureLoc.feature.name"]
		var refSeqName    = rowObj["featureLoc.referenceSequence.name"];
		glue.logInfo("Processing variation: "+variationName+" in reference: "+refSeqName+", feature: "+featureName);
	
	    variationTable.featureName = rowObj;

	}		
}

// Scan for variations
function scan_for_variations(variationTable) {

	for(var i = 0; i < refseqs_to_scan.length; i++) {
  
        var rowObj = rowObjs[i]
		var variationName = rowObj.name;
		var featureName   = rowObj["featureLoc.feature.name"]
		var refSeqName    = rowObj["featureLoc.referenceSequence.name"];
		glue.logInfo("Processing variation: "+variationName+" in reference: "+refSeqName+", feature: "+featureName);
	
		// Scan in alignment for variation
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

	}		
}

