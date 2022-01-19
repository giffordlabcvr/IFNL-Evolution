function calculateEntropy(refsequence, alignment, feature) {

	// set values for ambiguity codes 
	var codes = {N:['G','A','C','T'], B:['C', 'G', 'T'], D:['A', 'G', 'T'], H:['A','C','T'], V:['A', 'C', 'G'], K:['G', 'T'], Y:['C','T'], S:['C', 'G'], W:['A', 'T'], R:['A', 'G'], M:['A', 'C']};

	// create object to store counts
	var freqs = {};


	// export alignment from GLUE
	glue.inMode("module/fastaAlignmentExporter", function(){
	
		glue.logInfo(alignment);
	
		var sequences = glue.command(["export",alignment,"-r",refsequence,"-f",feature,"-a","-p"]);
		var list = sequences.nucleotideFasta.sequences;
	
	
		//loop through each sequence in the alignment
		_.each(list, function(seq){
	
		//loop through each position in the current sequence
			for (var i=0; i < seq.sequence.length; i++) {
				var key = "pos_" + (i + 1);
				var base = seq.sequence.charAt(i);
			
				//check the current base is a recognised character or exit with error
				if("ATCG-NBDHVKYSWRM".indexOf(base) < 0) {

					throw new Error("Wasn't expecting base = '"+base+"'");
				}
				if( base == "-" ) {
					base = "gap";
				}
				// set values for this position
				var positionFreqs = freqs[key];
				if(positionFreqs == null) {
					positionFreqs = {A:0, T:0, C:0, G:0, gap:0};
					freqs[key] = positionFreqs;
				}
			
				// if the current base is an ambiguity code
				if (base in codes) {
					options = codes[base];
					value = 1 / (codes[base].length);
				
					// increment the value for the current base
					_.each(options, function(opt) {
				
						positionFreqs[opt] += value;					
					});
				
				} else {
					positionFreqs[base] ++;
				}
			}
		});
	});

	entropyValues = {}
	i = 0

	// loop through stored position counts and calculate entropy
	_.each (freqs, function(pos){

		var total = pos['A'] + pos['G'] + pos['C'] + pos['T'];
	
		var v1 = -(pos['A'] / total) * Math.log(pos['A'] / total);
		var v2 = -(pos['G'] / total) * Math.log(pos['G'] / total);
		var v3 = -(pos['C'] / total) * Math.log(pos['C'] / total);
		var v4 = -(pos['T'] / total) * Math.log(pos['T'] / total);
	
		// set value to 0 if variable does not contain a valid number (i.e. ln(0) = -infinity)
		v1 = v1 || 0 ;
		v2 = v2 || 0 ;
		v3 = v3 || 0 ;
		v4 = v4 || 0 ;
	
		var entropy = (v1 + v2 + v3 + v4).toFixed(3);
	
		var key2 = i + 1;
		entropyValues[key2] = entropy;
		i++;

	});

	//show entropy output in GLUE
	var entropyArray = [];
	_.each(_.pairs(entropyValues), function(pair) {
		var key = pair[0];
		var value = pair[1];
		entropyArray.push({
			position: key,
			entropy: value
		});
	});

	return entropyArray;
}



