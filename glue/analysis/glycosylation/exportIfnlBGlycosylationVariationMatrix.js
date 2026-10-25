//
// glue/analysis/exportIfnlBGlycosylationVariationMatrix.js
//
// Scan predefined N-linked glycosylation Variations across the
// mammalian IFNL-B alignment and export a presence/absence matrix.
//
// Output values:
//   1  = variation present
//   0  = variation absent, with sufficient coverage
//   NA = insufficient coverage
//
// Run from:
//
//   /project/ifnl
//
// with:
//
//   run script glue/analysis/exportIfnlBGlycosylationVariationMatrix.js
//

var ALIGNMENT_NAME = "AL_IFNL_MAMMAL_B";
var REFERENCE_NAME = "REF_IFNL_Mammal_b_MASTER";
var FEATURE_NAME = "orf";

var OUTPUT_PATH = "tabular/analysis/ifnl_b_glycosylation_variations.tsv";

var VARIATIONS = [
    "n-linked-glycosylation-o1",
    "n-linked-glycosylation-o2",
    "n-linked-glycosylation-o3",
    "n-linked-glycosylation-o4",
    "n-linked-glycosylation-o5",
    "n-linked-glycosylation-o6"
];

var File = Java.type("java.io.File");
var FileWriter = Java.type("java.io.FileWriter");
var BufferedWriter = Java.type("java.io.BufferedWriter");


// -----------------------------------------------------------------------------
// Utility functions
// -----------------------------------------------------------------------------

function ensureParentDir(path) {
    var file = new File(path);
    var parent = file.getParentFile();

    if(parent !== null && !parent.exists()) {
        parent.mkdirs();
    }
}


function openWriter(path) {
    return new BufferedWriter(new FileWriter(path));
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


// -----------------------------------------------------------------------------
// Load alignment members
// -----------------------------------------------------------------------------

var members = [];

glue.inMode("alignment/" + ALIGNMENT_NAME, function() {

    members = glue.tableToObjects(
        glue.command([
            "list", "member",
            "sequence.source.name",
            "sequence.sequenceID"
        ])
    );
});

glue.logInfo(
    "Found " + members.length +
    " members in " + ALIGNMENT_NAME
);


// -----------------------------------------------------------------------------
// Build index
// -----------------------------------------------------------------------------
//
// Key:
//     sourceName/sequenceID
//
// Value:
//     {
//         sourceName: "...",
//         sequenceID: "...",
//         variation calls...
//     }
//

var memberIndex = {};

_.each(members, function(member) {

    var sourceName = member["sequence.source.name"];
    var sequenceID = member["sequence.sequenceID"];

    var key = sourceName + "/" + sequenceID;

    memberIndex[key] = {
        sourceName: sourceName,
        sequenceID: sequenceID
    };
});


// -----------------------------------------------------------------------------
// Scan each glycosylation Variation
// -----------------------------------------------------------------------------

_.each(VARIATIONS, function(variationName) {

    glue.logInfo(
        "Scanning variation " + variationName
    );

    var scanRows = [];

    glue.inMode("alignment/" + ALIGNMENT_NAME, function() {

        scanRows = glue.tableToObjects(
            glue.command([
                "variation", "member", "scan",
                "-r", REFERENCE_NAME,
                "-f", FEATURE_NAME,
                "-v", variationName
            ])
        );
    });

    glue.logInfo(
        "Variation " + variationName +
        ": received " + scanRows.length + " scan rows"
    );

    _.each(scanRows, function(scanRow) {

        var sourceName = scanRow["sourceName"];
        var sequenceID = scanRow["sequenceID"];

        var key = sourceName + "/" + sequenceID;

        var memberObj = memberIndex[key];

        if(memberObj === undefined) {

            glue.logInfo(
                "WARNING: scan returned unknown member " + key
            );

            return;
        }

        var sufficientCoverage =
            scanRow["sufficientCoverage"];

        var present =
            scanRow["present"];

        //
        // GLUE may return booleans either as actual booleans
        // or as string values depending on result handling.
        //

        var coverageOK =
            sufficientCoverage === true ||
            String(sufficientCoverage).toLowerCase() === "true";

        var isPresent =
            present === true ||
            String(present).toLowerCase() === "true";

        if(!coverageOK) {
            memberObj[variationName] = "NA";
        } else if(isPresent) {
            memberObj[variationName] = "1";
        } else {
            memberObj[variationName] = "0";
        }
    });

});


// -----------------------------------------------------------------------------
// Fill in any missing scan values
// -----------------------------------------------------------------------------
//
// Normally every member should appear in each scan.
// If it does not, record NA rather than silently treating it as absent.
//

_.each(_.values(memberIndex), function(memberObj) {

    _.each(VARIATIONS, function(variationName) {

        if(memberObj[variationName] === undefined) {
            memberObj[variationName] = "NA";
        }
    });
});


// -----------------------------------------------------------------------------
// Sort output
// -----------------------------------------------------------------------------

var outputRows = _.values(memberIndex);

outputRows.sort(function(a, b) {

    var keyA =
        a.sourceName + "/" + a.sequenceID;

    var keyB =
        b.sourceName + "/" + b.sequenceID;

    return keyA.localeCompare(keyB);
});


// -----------------------------------------------------------------------------
// Write TSV
// -----------------------------------------------------------------------------

ensureParentDir(OUTPUT_PATH);

var writer = openWriter(OUTPUT_PATH);

try {

    var header = [
        "source_name",
        "sequence_id"
    ].concat(VARIATIONS);

    writeLine(writer, header);

    _.each(outputRows, function(row) {

        var fields = [
            row.sourceName,
            row.sequenceID
        ];

        _.each(VARIATIONS, function(variationName) {
            fields.push(row[variationName]);
        });

        writeLine(writer, fields);
    });

} finally {

    writer.close();
}


// -----------------------------------------------------------------------------
// Summary
// -----------------------------------------------------------------------------

glue.logInfo(
    "Exported glycosylation variation matrix for " +
    outputRows.length + " IFNL-B alignment members"
);

glue.logInfo(
    "Output: " + OUTPUT_PATH
);
