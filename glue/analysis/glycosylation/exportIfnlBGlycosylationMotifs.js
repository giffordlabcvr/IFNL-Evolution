//
// exportIfnlBGlycosylationMotifs.js
//
// Extract the exact N-linked glycosylation motifs detected within the
// six predefined IFNL-B glycosylation windows.
//
// The aligned protein sequence is exported relative to
// REF_IFNL_Mammal_b_MASTER, so each character position corresponds
// directly to a labeled ORF codon on the MASTER reference.
//
// Outputs:
//
//   tabular/analysis/ifnl_b_glycosylation_motifs.tsv
//       One row per detected motif.
//
//   tabular/analysis/ifnl_b_glycosylation_motif_summary.tsv
//       Counts of each region / codon / motif combination.
//
// Run from:
//
//   /project/ifnl
//
// with:
//
//   run script glue/analysis/exportIfnlBGlycosylationMotifs.js
//

var ALIGNMENT_NAME = "AL_IFNL_MAMMAL_B";
var REFERENCE_NAME = "REF_IFNL_Mammal_b_MASTER";
var FEATURE_NAME = "orf";

var EXPORTER_MODULE = "fastaProteinAlignmentExporter";

var OUTPUT_PATH =
    "tabular/analysis/ifnl_b_glycosylation_motifs.tsv";

var SUMMARY_PATH =
    "tabular/analysis/ifnl_b_glycosylation_motif_summary.tsv";


var WINDOWS = [
    {
        region: "o1",
        variation: "n-linked-glycosylation-o1",
        codonStart: 30,
        codonEnd: 39
    },
    {
        region: "o2",
        variation: "n-linked-glycosylation-o2",
        codonStart: 40,
        codonEnd: 50
    },
    {
        region: "o3",
        variation: "n-linked-glycosylation-o3",
        codonStart: 55,
        codonEnd: 65
    },
    {
        region: "o4",
        variation: "n-linked-glycosylation-o4",
        codonStart: 70,
        codonEnd: 80
    },
    {
        region: "o5",
        variation: "n-linked-glycosylation-o5",
        codonStart: 100,
        codonEnd: 109
    },
    {
        region: "o6",
        variation: "n-linked-glycosylation-o6",
        codonStart: 110,
        codonEnd: 120
    }
];


var File = Java.type("java.io.File");
var FileWriter = Java.type("java.io.FileWriter");
var BufferedWriter = Java.type("java.io.BufferedWriter");


// -----------------------------------------------------------------------------
// Utilities
// -----------------------------------------------------------------------------

function ensureParentDir(path) {

    var file = new File(path);
    var parent = file.getParentFile();

    if(parent !== null && !parent.exists()) {
        parent.mkdirs();
    }
}


function tsvEscape(value) {

    if(value === null || value === undefined) {
        return "";
    }

    return String(value)
        .replace(/\t/g, " ")
        .replace(/\r?\n/g, " ");
}


function writeLine(writer, fields) {

    writer.write(
        fields.map(tsvEscape).join("\t")
    );

    writer.newLine();
}


//
// Exporter IDs have the form:
//
//   AL_IFNL_MAMMAL_B.ncbi-curated-mammalia-b.HQ201955
//
// Strip the alignment prefix, then split once to recover sourceName
// and sequenceID.
//
function parseExporterId(id) {

    var prefix = ALIGNMENT_NAME + ".";

    var remainder = id;

    if(remainder.indexOf(prefix) === 0) {
        remainder = remainder.substring(prefix.length);
    }

    var firstDot = remainder.indexOf(".");

    if(firstDot < 0) {

        return {
            sourceName: "",
            sequenceID: remainder
        };
    }

    return {
        sourceName: remainder.substring(0, firstDot),
        sequenceID: remainder.substring(firstDot + 1)
    };
}


//
// Test whether three aligned amino-acid characters form an
// N-linked glycosylation consensus:
//
//     N-X-S/T
//
// where X is not proline.
//
// Gaps are not allowed.
// The middle position may contain any represented amino-acid
// character except P or '-'.
//
function isGlycosylationMotif(a1, a2, a3) {

    if(a1 !== "N") {
        return false;
    }

    if(a2 === "P" || a2 === "-") {
        return false;
    }

    if(a3 !== "S" && a3 !== "T") {
        return false;
    }

    return true;
}


// -----------------------------------------------------------------------------
// Extract sequence records from exporter result
// -----------------------------------------------------------------------------

function getExportedSequences(result) {

    if(result === null || result === undefined) {
        throw "Protein alignment exporter returned no result";
    }

    if(result.aminoAcidFasta === undefined) {
        throw "Unexpected exporter result: aminoAcidFasta not found";
    }

    if(result.aminoAcidFasta.sequences === undefined) {
        throw "Unexpected exporter result: aminoAcidFasta.sequences not found";
    }

    return result.aminoAcidFasta.sequences;
}


// -----------------------------------------------------------------------------
// Scan windows
// -----------------------------------------------------------------------------

var motifRows = [];


//
// Key:
//
//   region|motifStart|motif
//
// Value:
//   count
//
var summaryCounts = {};


_.each(WINDOWS, function(window) {

    glue.logInfo(
        "Exporting " +
        window.region +
        " (codons " +
        window.codonStart +
        "-" +
        window.codonEnd +
        ")"
    );

	var result;

	glue.inMode(
		"module/" + EXPORTER_MODULE,
		function() {

			result = glue.command([
				"export",
				ALIGNMENT_NAME,
				"-r", REFERENCE_NAME,
				"-f", FEATURE_NAME,
				"-l",
				String(window.codonStart),
				String(window.codonEnd),
				"-a",
				"-p"
			]);
		}
	);

    var sequences = getExportedSequences(result);

    glue.logInfo(
        "Received " +
        sequences.length +
        " aligned protein sequences for " +
        window.region
    );


    _.each(sequences, function(record) {

        var parsedId =
            parseExporterId(record.id);

        var aaSequence =
            String(record.sequence).toUpperCase();

        var expectedLength =
            window.codonEnd -
            window.codonStart + 1;

        if(aaSequence.length !== expectedLength) {

            glue.logInfo(
                "WARNING: " +
                record.id +
                " returned " +
                aaSequence.length +
                " aa for " +
                window.region +
                "; expected " +
                expectedLength
            );
        }


        //
        // Keep the aligned sequence intact.
        //
        // Character i corresponds to:
        //
        //     MASTER codon = codonStart + i
        //
        // Do NOT remove gaps.
        //
        for(var i = 0; i <= aaSequence.length - 3; i++) {

            var a1 = aaSequence.charAt(i);
            var a2 = aaSequence.charAt(i + 1);
            var a3 = aaSequence.charAt(i + 2);

            if(!isGlycosylationMotif(a1, a2, a3)) {
                continue;
            }

            var motif =
                a1 + a2 + a3;

            var motifStart =
                window.codonStart + i;

            var motifEnd =
                motifStart + 2;


            motifRows.push({
                sourceName: parsedId.sourceName,
                sequenceID: parsedId.sequenceID,
                region: window.region,
                variation: window.variation,
                motifStart: motifStart,
                motifEnd: motifEnd,
                motif: motif,
                alignedWindow: aaSequence
            });


            var summaryKey =
                window.region + "|" +
                motifStart + "|" +
                motif;

            if(summaryCounts[summaryKey] === undefined) {

                summaryCounts[summaryKey] = {
                    region: window.region,
                    motifStart: motifStart,
                    motifEnd: motifEnd,
                    motif: motif,
                    count: 0
                };
            }

            summaryCounts[summaryKey].count++;
        }
    });
});


// -----------------------------------------------------------------------------
// Sort motif rows
// -----------------------------------------------------------------------------

motifRows.sort(function(a, b) {

    if(a.region !== b.region) {
        return a.region.localeCompare(b.region);
    }

    if(a.motifStart !== b.motifStart) {
        return a.motifStart - b.motifStart;
    }

    if(a.motif !== b.motif) {
        return a.motif.localeCompare(b.motif);
    }

    if(a.sourceName !== b.sourceName) {
        return a.sourceName.localeCompare(b.sourceName);
    }

    return a.sequenceID.localeCompare(b.sequenceID);
});


// -----------------------------------------------------------------------------
// Write long-form motif table
// -----------------------------------------------------------------------------

ensureParentDir(OUTPUT_PATH);

var writer =
    new BufferedWriter(
        new FileWriter(OUTPUT_PATH)
    );

try {

    writeLine(writer, [
        "source_name",
        "sequence_id",
        "region",
        "variation",
        "motif_start",
        "motif_end",
        "motif",
        "aligned_window"
    ]);


    _.each(motifRows, function(row) {

        writeLine(writer, [
            row.sourceName,
            row.sequenceID,
            row.region,
            row.variation,
            row.motifStart,
            row.motifEnd,
            row.motif,
            row.alignedWindow
        ]);
    });

} finally {

    writer.close();
}


// -----------------------------------------------------------------------------
// Write motif summary
// -----------------------------------------------------------------------------

var summaryRows =
    _.values(summaryCounts);


summaryRows.sort(function(a, b) {

    if(a.region !== b.region) {
        return a.region.localeCompare(b.region);
    }

    if(a.motifStart !== b.motifStart) {
        return a.motifStart - b.motifStart;
    }

    if(a.motif !== b.motif) {
        return a.motif.localeCompare(b.motif);
    }

    return b.count - a.count;
});


ensureParentDir(SUMMARY_PATH);

var summaryWriter =
    new BufferedWriter(
        new FileWriter(SUMMARY_PATH)
    );

try {

    writeLine(summaryWriter, [
        "region",
        "motif_start",
        "motif_end",
        "motif",
        "count"
    ]);


    _.each(summaryRows, function(row) {

        writeLine(summaryWriter, [
            row.region,
            row.motifStart,
            row.motifEnd,
            row.motif,
            row.count
        ]);
    });

} finally {

    summaryWriter.close();
}


// -----------------------------------------------------------------------------
// Report
// -----------------------------------------------------------------------------

glue.logInfo(
    "Detected " +
    motifRows.length +
    " glycosylation motif occurrence(s)"
);

glue.logInfo(
    "Observed " +
    summaryRows.length +
    " distinct region / coordinate / motif combinations"
);

glue.logInfo(
    "Motif table: " +
    OUTPUT_PATH
);

glue.logInfo(
    "Summary table: " +
    SUMMARY_PATH
);

