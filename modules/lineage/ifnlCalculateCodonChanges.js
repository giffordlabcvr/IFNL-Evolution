var refsequence = "REF_MASTER_B_AF033819";
var refname = 'AF033819|HIV-1.M.B'
var alignment   = "AL_HIV-1M_B";
var feature = 'vpr';

var aa_codons = {ATT:'I', ATC:'I', ATA:'I', CTT:'L', CTC:'L', CTA:'L', CTG:'L', TTA:'L', TTG:'L', GTT:'V', GTC:'V', GTA:'V', GTG:'V', TTT:'F', TTC:'F', ATG:'M', TGT:'C', TGC:'C', GCT:'A', GCC:'A', GCA:'A', GCG:'A', GGT:'G', GGC:'G', GGA:'G', GGG:'G', CCT:'P', CCC:'P', CCA:'P', CCG:'P', ACT:'T', ACC:'T', ACA:'T', ACG:'T', TCT:'S', TCC:'S', TCA:'S', TCG:'S', AGT:'S', AGC:'S', TAT:'Y', TAC:'Y', TGG:'W', CAA:'Q', CAG:'Q', AAT:'N', AAC:'N', CAT:'H', CAC:'H', GAA:'E', GAG:'E', GAT:'D', GAC:'D', AAA:'K', AAG:'K', CGT:'R', CGC:'R', CGA:'R', CGG:'R', AGA:'R', AGG:'R', TAA:'*', TAG:'*', TGA:'*'};
var ref_seq = '';
var codon_counts = {};

// export alignment from GLUE
glue.inMode("module/hivFastaAlignmentExporter", function(){	
	var sequences = glue.command(["export",alignment,"-r",refsequence,"-f",feature,"-a","-p"]);
	var list = sequences.nucleotideFasta.sequences;
	
	// loop through each sequence in the alignment to find reference sequence
	_.each(list, function(seq){			
		if (seq.id == refname ) {
			ref_seq = seq.sequence;
		}	
	});			
// loop through each sequence in the alignment and compare with ref
	_.each(list, function(seq){
		var position = 1;
		for (var i=0; i < seq.sequence.length; i+=3) {
			//var key = "codon_" + position;
			var key = position;
			
			var seq_codon = seq.sequence.charAt(i) +seq.sequence.charAt(i+1) +seq.sequence.charAt(i+2);
			var ref_codon = ref_seq.charAt(i) + ref_seq.charAt(i+1) + ref_seq.charAt(i+3);
		
			if (seq_codon != ref_codon) {				
				
				//glue.logInfo(key + "\t" + seq_codon + "\t" + ref_codon);
				
				var seq_aa = aa_codons[seq_codon];
				var ref_aa = aa_codons[ref_codon];

				//glue.logInfo(seq.id + "\t" + ref_aa + "\t" + seq_aa);
							
				var positionFreqs = codon_counts[key];
				if (positionFreqs == null) {
					positionFreqs = {S:0, NS:0, T:0};
					codon_counts[key] = positionFreqs;
				}
				opt = ""
				if (seq_aa != ref_aa){
					//increment counter for NS
					opt = "NS";
				} else if (seq_aa == ref_aa) {
					//increment counter for S
					opt = "S";
				}
				positionFreqs[opt]++;
				positionFreqs["T"]++;
			}

			position++;
		}
	});

//output = codon | synonymous changes  | non-synonymous changes  | total changes
glue.logInfo("codon_counts", codon_counts);	
});