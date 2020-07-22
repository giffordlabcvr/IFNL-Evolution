// List of genes to scan for
var alignments = ["IFNLb1", "IFNLb2" ];

// Get list of variations to scan for 
// list variation -w "name like 'n-linked-glycosylation%'" 
var variationTable;

glue.logInfo("Processing sequence ID: "+id);

glue.inMode("/alignment/"+alignment+"/"+alignmentName, function() {

	var rowObjs = glue.tableToObjects(glue.command(["amino-acid","-f",gene,"-r"+refsequence]));
	for(var i = 0; i < rowObjs.length; i++) {

		var codon = rowObjs[i].codonNts;
		var codonLabel = rowObjs[i].codonLabel;
		var aminoAcid  = rowObjs[i].aminoAcid;
		var relRefNt   = rowObjs[i].relRefNt;
		glue.logInfo("Processing sequence ID: "+id+"Processing gene: "+gene+", got codon: "+codon+" Amino Acid "+aminoAcid);
		var resultObj = rowObjs[i];
		var codonIndex = gene + codonLabel;
		geneResultMap[codonIndex] = resultObj;
	
	}		
});
