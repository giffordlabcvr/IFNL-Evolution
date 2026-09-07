//
// populateIfnlBGlycosylationSites.js
//
// Populate sequence-level predicted N-linked glycosylation sequon fields
// for members of AL_IFNL_MAMMAL_B using position-specific amino-acid
// variation scans.
//
// Run from:
//
//   /project/ifnl
//
// with:
//
//   run script glue/analysis/populateIfnlBGlycosylationSites.js
//

var ALIGNMENT_NAME = "AL_IFNL_MAMMAL_B";
var REFERENCE_NAME = "REF_IFNL_Mammal_b_MASTER";
var FEATURE_NAME = "orf";

var VARIATIONS = [

    {
        variationName: "n-linked-glycosylation-o1-35",
        fieldName: "glycosylation_o1_35"
    },

    {
        variationName: "n-linked-glycosylation-o2-42",
        fieldName: "glycosylation_o2_42"
    },

    {
        variationName: "n-linked-glycosylation-o3-59",
        fieldName: "glycosylation_o3_59"
    },

    {
        variationName: "n-linked-glycosylation-o4-72",
        fieldName: "glycosylation_o4_72"
    },

    {
        variationName: "n-linked-glycosylation-o4-74",
        fieldName: "glycosylation_o4_74"
    },

    {
        variationName: "n-linked-glycosylation-o4-76",
        fieldName: "glycosylation_o4_76"
    },

    {
        variationName: "n-linked-glycosylation-o5-103",
        fieldName: "glycosylation_o5_103"
    },

    {
        variationName: "n-linked-glycosylation-o6-111",
        fieldName: "glycosylation_o6_111"
    },

    {
        variationName: "n-linked-glycosylation-o6-115",
        fieldName: "glycosylation_o6_115"
    },

    {
        variationName: "n-linked-glycosylation-o6-116",
        fieldName: "glycosylation_o6_116"
    }
];


//
// Convert a GLUE result value to boolean.
//
function isTrue(value) {

    if(value === true) {
        return true;
    }

    if(value === false ||
       value === null ||
       value === undefined) {
        return false;
    }

    return String(value).toLowerCase() === "true";
}


//
// Populate one field on a sequence.
//
function setSequenceField(sourceName, sequenceID, fieldName, value) {

    glue.inMode(
        "sequence/" + sourceName + "/" + sequenceID,
        function() {

            glue.command([
                "set", "field",
                fieldName,
                value
            ]);
        }
    );
}


//
// Scan each position-specific glycosylation variation.
//
_.each(VARIATIONS, function(variationDef) {

    var variationName = variationDef.variationName;
    var fieldName = variationDef.fieldName;

    glue.logInfo(
        "Scanning " +
        variationName +
        " -> sequence." +
        fieldName
    );

    var scanRows;

    glue.inMode(
        "alignment/" + ALIGNMENT_NAME,
        function() {

            scanRows = glue.tableToObjects(
                glue.command([
                    "variation", "member", "scan",
                    "-r", REFERENCE_NAME,
                    "-f", FEATURE_NAME,
                    "-v", variationName
                ])
            );
        }
    );

    var presentCount = 0;
    var absentCount = 0;
    var insufficientCount = 0;

    _.each(scanRows, function(row) {

        var sourceName = row["sourceName"];
        var sequenceID = row["sequenceID"];

        var sufficientCoverage =
            isTrue(row["sufficientCoverage"]);

        var present =
            isTrue(row["present"]);

        var value;

        if(!sufficientCoverage) {

            value = "insufficient_coverage";
            insufficientCount++;

        } else if(present) {

            value = "present";
            presentCount++;

        } else {

            value = "absent";
            absentCount++;
        }

        setSequenceField(
            sourceName,
            sequenceID,
            fieldName,
            value
        );
    });

    glue.logInfo(
        variationName +
        ": present=" + presentCount +
        ", absent=" + absentCount +
        ", insufficient_coverage=" + insufficientCount
    );
});


glue.logInfo(
    "Finished populating IFNL-B glycosylation sequon fields"
);