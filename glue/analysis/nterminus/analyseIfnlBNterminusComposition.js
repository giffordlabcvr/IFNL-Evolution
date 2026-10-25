//
// analyseIfnlBNterminusComposition.js
//
// Exploratory analysis of cysteine / tryptophan composition in the
// N-terminal region of mammalian IFNL-B proteins.
//
// Each IFNL-B sequence has been promoted to a reference with a complete
// ORF feature-location. The script therefore uses:
//
//   reference <REF> feature-location orf amino-acid
//
// to obtain the translated protein directly from each reference.
//
// Required sequence-table fields (create once in schema-project mode):
//   table sequence
//   create field mammalian_b_cw_profile VARCHAR 10
//   create field mammalian_b_cw_counts VARCHAR 10
//   exit
//
// Database output (first 20 translated ORF residues only):
//   mammalian_b_cw_profile: C_only, W_only, C_and_W, neither, unknown
//   mammalian_b_cw_counts: exact counts, e.g. C3_W0; unknown if uncallable.
//
// These are composition summaries, not N-terminal phylogenetic types or
// gene-conversion calls. An annotated ORF start is assumed, as in the original.
// Short windows and windows containing any nonstandard amino acid are unknown.
// Sequences outside the reference query are not changed.
//
// Set WRITE_DATABASE=false to produce the TSV without database writes.
// All analysis and duplicate checks complete before any sequence is updated.
// Writeback is per command, not an all-or-nothing transaction; reruns are safe.
//
// Output:
//   tabular/analysis/ifnl_b_n_terminus_composition.tsv
//
// Run from:
//   /project/ifnl
//
// with:
//   run script glue/analysis/analyseIfnlBNterminusComposition.js
//

var N_TERMINAL_LENGTH = 20;
var WRITE_DATABASE = true;

var PROFILE_FIELD = "mammalian_b_cw_profile";
var COUNTS_FIELD = "mammalian_b_cw_counts";

var OUTPUT_PATH =
    "tabular/analysis/ifnl_b_n_terminus_composition.tsv";

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


function countResidue(sequence, residue) {

    var count = 0;

    for(var i = 0; i < sequence.length; i++) {

        if(sequence.charAt(i) === residue) {
            count++;
        }
    }

    return count;
}


// Build profiles only for complete, unambiguous 20-aa windows.
// No thresholds for "richness" or predominance are imposed.
function cwProfile(nTerminus) {
    if(nTerminus.length !== N_TERMINAL_LENGTH) {
        return {profile: "unknown", counts: "unknown", status: "short_window"};
    }
    if(!/^[ACDEFGHIKLMNPQRSTVWY]+$/.test(nTerminus)) {
        return {profile: "unknown", counts: "unknown", status: "ambiguous_window"};
    }
    var c = countResidue(nTerminus, "C");
    var w = countResidue(nTerminus, "W");
    var profile = c > 0 ? (w > 0 ? "C_and_W" : "C_only") :
                         (w > 0 ? "W_only" : "neither");
    return {profile: profile, counts: "C" + c + "_W" + w, status: "called"};
}


// -----------------------------------------------------------------------------
// Get mammalian IFNL-B references
// -----------------------------------------------------------------------------

var references = glue.tableToObjects(
    glue.command([
        "list", "reference",
        "name",
        "sequence.source.name",
        "sequence.sequenceID",
        "sequence.species",
        "sequence.species_order",
        "sequence." + PROFILE_FIELD,
        "sequence." + COUNTS_FIELD,
        "-w",
        "sequence.gene_name = 'IFNL-Mammalia-B'"
    ])
);

glue.logInfo(
    "Found " +
    references.length +
    " mammalian IFNL-B references"
);


// -----------------------------------------------------------------------------
// Analyse each reference
// -----------------------------------------------------------------------------

var outputRows = [];

_.each(references, function(refRow) {

    var referenceName =
        refRow["name"];

    var sourceName =
        refRow["sequence.source.name"];

    var sequenceID =
        refRow["sequence.sequenceID"];

    var species =
        refRow["sequence.species"];

    var speciesOrder =
        refRow["sequence.species_order"];


    glue.logInfo(
        "Analysing " +
        referenceName +
        " (" +
        sequenceID +
        ")"
    );


    var aaRows;

    glue.inMode(
        "reference/" +
        referenceName +
        "/feature-location/orf",
        function() {

            aaRows = glue.tableToObjects(
                glue.command([
                    "amino-acid"
                ])
            );
        }
    );


    //
    // Reconstruct translated ORF.
    //
    var protein = "";
    var nTerminalTranslationIssue = false;

    _.each(aaRows, function(row) {

        var aa = row["aminoAcid"];
        if(protein.length < N_TERMINAL_LENGTH &&
           (aa === null || aa === undefined || String(aa).length !== 1)) {
            nTerminalTranslationIssue = true;
        }

        if(aa !== null &&
           aa !== undefined &&
           String(aa).length > 0) {

            protein += String(aa);
        }
    });


    //
    // Extract first N residues.
    //
    var actualNTermLength =
        Math.min(
            N_TERMINAL_LENGTH,
            protein.length
        );

    var nTerminus =
        protein.substring(
            0,
            actualNTermLength
        );


    //
    // Count C and W.
    //
    var cCount =
        countResidue(
            nTerminus,
            "C"
        );

    var wCount =
        countResidue(
            nTerminus,
            "W"
        );

    var cwTotal =
        cCount + wCount;

    var cMinusW =
        cCount - wCount;


    //
    // Also calculate counts in the rest of the protein.
    //
    var remainder =
        protein.substring(
            actualNTermLength
        );

    var remainderCCount =
        countResidue(
            remainder,
            "C"
        );

    var remainderWCount =
        countResidue(
            remainder,
            "W"
        );


    var profileResult = cwProfile(nTerminus);
    if(nTerminalTranslationIssue) {
        profileResult = {profile: "unknown", counts: "unknown",
                         status: "missing_or_invalid_translation_row"};
    }

    outputRows.push({
        cwProfile: profileResult.profile,
        cwCounts: profileResult.counts,
        cwStatus: profileResult.status,

        referenceName: referenceName,
        sourceName: sourceName,
        sequenceID: sequenceID,

        species: species,
        speciesOrder: speciesOrder,

        proteinLength: protein.length,
        nTerminalLength: actualNTermLength,
        nTerminalSequence: nTerminus,

        cCount: cCount,
        wCount: wCount,
        cwTotal: cwTotal,
        cMinusW: cMinusW,

        remainderCCount: remainderCCount,
        remainderWCount: remainderWCount
    });
});


// -----------------------------------------------------------------------------
// Sort output
// -----------------------------------------------------------------------------

outputRows.sort(function(a, b) {

    var orderA =
        (a.speciesOrder || "") +
        "|" +
        (a.species || "") +
        "|" +
        a.sequenceID;

    var orderB =
        (b.speciesOrder || "") +
        "|" +
        (b.species || "") +
        "|" +
        b.sequenceID;

    return orderA.localeCompare(orderB);
});


// Deduplicate references to the same sequence. Conflicting results abort before
// database writes. JSON keys avoid collisions between source/sequence ID pairs.
var uniqueRows = [];
var seenSequences = {};
_.each(outputRows, function(row) {
    var key = JSON.stringify([row.sourceName, row.sequenceID]);
    var prior = seenSequences[key];
    if(prior !== undefined) {
        if(prior.cwProfile !== row.cwProfile ||
           prior.cwCounts !== row.cwCounts ||
           prior.nTerminalSequence !== row.nTerminalSequence) {
            throw new Error("Conflicting reference translations for " + key);
        }
    } else {
        seenSequences[key] = row;
        uniqueRows.push(row);
    }
});


// -----------------------------------------------------------------------------
// Write TSV
// -----------------------------------------------------------------------------

ensureParentDir(
    OUTPUT_PATH
);

var writer =
    new BufferedWriter(
        new FileWriter(
            OUTPUT_PATH
        )
    );

try {

    writeLine(writer, [

        "reference_name",
        "source_name",
        "sequence_id",

        "species",
        "species_order",

        "protein_length",
        "n_terminal_length",
        "n_terminal_sequence",

        "c_count",
        "w_count",
        "cw_total",
        "c_minus_w",

        "remainder_c_count",
        "remainder_w_count",
        "cw_profile",
        "cw_counts",
        "cw_status"
    ]);


    _.each(
        outputRows,
        function(row) {

            writeLine(writer, [

                row.referenceName,
                row.sourceName,
                row.sequenceID,

                row.species,
                row.speciesOrder,

                row.proteinLength,
                row.nTerminalLength,
                row.nTerminalSequence,

                row.cCount,
                row.wCount,
                row.cwTotal,
                row.cMinusW,

                row.remainderCCount,
                row.remainderWCount,
                row.cwProfile,
                row.cwCounts,
                row.cwStatus
            ]);
        }
    );

} finally {

    writer.close();
}


// Write to the sequence record, not the reference record.
// Unknown results overwrite old values so stale profiles cannot survive reruns.
var written = 0;
if(WRITE_DATABASE) {
    _.each(uniqueRows, function(row) {
        glue.inMode("sequence/" + row.sourceName + "/" + row.sequenceID, function() {
            glue.command(["set", "field", PROFILE_FIELD, row.cwProfile]);
            glue.command(["set", "field", COUNTS_FIELD, row.cwCounts]);
        });
        written++;
    });
}
glue.logInfo("Database writeback: " +
    (WRITE_DATABASE ? written + " unique sequences updated" : "disabled"));


// -----------------------------------------------------------------------------
// Simple summary
// -----------------------------------------------------------------------------

var cOnly = 0;
var wOnly = 0;
var both = 0;
var neither = 0;
var unknown = 0;

_.each(outputRows, function(row) {

    if(row.cwProfile === "unknown") {
        unknown++;
        return;
    }

    if(row.cCount > 0 &&
       row.wCount === 0) {

        cOnly++;

    } else if(row.wCount > 0 &&
              row.cCount === 0) {

        wOnly++;

    } else if(row.cCount > 0 &&
              row.wCount > 0) {

        both++;

    } else {

        neither++;
    }
});


glue.logInfo(
    "N-terminal analysis complete"
);

glue.logInfo(
    "C present / W absent: " +
    cOnly
);

glue.logInfo(
    "W present / C absent: " +
    wOnly
);

glue.logInfo(
    "Both C and W present: " +
    both
);

glue.logInfo(
    "Neither C nor W present: " +
    neither
);

glue.logInfo(
    "Output: " +
    OUTPUT_PATH
);


glue.logInfo("Unknown (short or ambiguous N terminus): " + unknown);
