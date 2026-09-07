//
// summariseIfnlBGlycosylationContext.js
//
// Summarise the genomic, structural and functional context of the
// predefined IFNL-B N-linked glycosylation variation windows.
//
// Uses:
//   REF_IFNL_Mammal_b_MASTER
//
// Glycosylation variations are defined in ORF codon-label coordinates.
// Existing feature-location annotations are defined in reference nucleotide
// coordinates.
//
// The script therefore:
//
//   1. obtains the ORF codon -> reference nucleotide mapping from GLUE;
//   2. obtains reference-coordinate segments for annotated features;
//   3. determines which features overlap each glycosylation window;
//   4. exports a TSV summary.
//
// Run from:
//
//   /project/ifnl
//
// with:
//
//   run script glue/analysis/summariseIfnlBGlycosylationContext.js
//

var REFERENCE_NAME = "REF_IFNL_Mammal_b_MASTER";
var ORF_FEATURE = "orf";

var OUTPUT_PATH =
    "tabular/analysis/ifnl_b_glycosylation_context.tsv";


//
// Glycosylation variation windows.
// These should correspond exactly to the definitions in the build.
//
var GLYCO_WINDOWS = [
    {
        name: "n-linked-glycosylation-o1",
        shortName: "o1",
        codonStart: 30,
        codonEnd: 39
    },
    {
        name: "n-linked-glycosylation-o2",
        shortName: "o2",
        codonStart: 40,
        codonEnd: 50
    },
    {
        name: "n-linked-glycosylation-o3",
        shortName: "o3",
        codonStart: 55,
        codonEnd: 65
    },
    {
        name: "n-linked-glycosylation-o4",
        shortName: "o4",
        codonStart: 70,
        codonEnd: 80
    },
    {
        name: "n-linked-glycosylation-o5",
        shortName: "o5",
        codonStart: 100,
        codonEnd: 109
    },
    {
        name: "n-linked-glycosylation-o6",
        shortName: "o6",
        codonStart: 110,
        codonEnd: 120
    }
];


//
// Features whose context we want to report.
//
var FEATURE_GROUPS = {

    exon: [
        "mRNA_exon1",
        "mRNA_exon2",
        "mRNA_exon3",
        "mRNA_exon4",
        "mRNA_exon5",
        "mRNA_exon6"
    ],

    protein_region: [
        "SP",
        "N-helix",
        "helix-b",
        "helix-c",
        "helix-d",
        "helix-e",
        "helix-f",
        "loop1",
        "loop2",
        "loop3",
        "loop4",
        "loop5",
        "loop6",
        "c-term"
    ],

    receptor_context: [
        "receptor_binding",
        "non_receptor_binding",
        "ifnlr1_binding_surface",
        "ifnlr2_binding_surface",
        "il10r1r2_binding_surface"
    ],

    structural_context: [
        "exposed_residues",
        "exposed_nonbinding_residues",
        "buried_residues"
    ]
};


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


function writeLine(writer, values) {

    writer.write(
        values.map(tsvEscape).join("\t")
    );

    writer.newLine();
}


function rangesOverlap(start1, end1, start2, end2) {

    return start1 <= end2 && start2 <= end1;
}


// -----------------------------------------------------------------------------
// Obtain ORF amino-acid coordinate mapping
// -----------------------------------------------------------------------------

var aaRows;

glue.inMode(
    "reference/" + REFERENCE_NAME +
    "/feature-location/" + ORF_FEATURE,
    function() {

        aaRows = glue.tableToObjects(
            glue.command([
                "amino-acid"
            ])
        );
    }
);

glue.logInfo(
    "Retrieved " + aaRows.length +
    " ORF amino-acid rows from " + REFERENCE_NAME
);


//
// Index rows by codon label.
//
var codonIndex = {};

_.each(aaRows, function(row) {

    var codonLabel = parseInt(row["codonLabel"]);
    var refNt = parseInt(row["refNt"]);

    if(!isNaN(codonLabel) && !isNaN(refNt)) {

        codonIndex[codonLabel] = {
            codonLabel: codonLabel,
            refNt: refNt,
            aminoAcid: row["aminoAcid"]
        };
    }
});


// -----------------------------------------------------------------------------
// Convert each glycosylation codon window to reference nucleotide coordinates
// -----------------------------------------------------------------------------

_.each(GLYCO_WINDOWS, function(window) {

    var firstCodon = codonIndex[window.codonStart];
    var lastCodon = codonIndex[window.codonEnd];

    if(firstCodon === undefined) {
        throw "No ORF mapping found for codon " +
            window.codonStart;
    }

    if(lastCodon === undefined) {
        throw "No ORF mapping found for codon " +
            window.codonEnd;
    }

    //
    // refNt is the first reference nucleotide of the codon.
    // The final codon therefore extends through refNt + 2.
    //
    window.refStart = firstCodon.refNt;
    window.refEnd = lastCodon.refNt + 2;

    glue.logInfo(
        window.shortName +
        ": codons " +
        window.codonStart + "-" + window.codonEnd +
        " -> reference nts " +
        window.refStart + "-" + window.refEnd
    );
});


// -----------------------------------------------------------------------------
// Retrieve feature-location segments
// -----------------------------------------------------------------------------

var featureSegments = {};

function loadFeatureSegments(featureName) {

    var rows;

    glue.inMode(
        "reference/" + REFERENCE_NAME +
        "/feature-location/" + featureName,
        function() {

            rows = glue.tableToObjects(
                glue.command([
                    "list", "segment"
                ])
            );
        }
    );

    featureSegments[featureName] = [];

    _.each(rows, function(row) {

        var refStart = parseInt(row["refStart"]);
        var refEnd = parseInt(row["refEnd"]);

        if(!isNaN(refStart) && !isNaN(refEnd)) {

            featureSegments[featureName].push({
                refStart: refStart,
                refEnd: refEnd
            });
        }
    });
}


//
// Load all features once.
//
_.each(_.keys(FEATURE_GROUPS), function(groupName) {

    _.each(FEATURE_GROUPS[groupName], function(featureName) {

        loadFeatureSegments(featureName);

        glue.logInfo(
            "Loaded " +
            featureSegments[featureName].length +
            " segment(s) for feature " +
            featureName
        );
    });
});


// -----------------------------------------------------------------------------
// Determine overlaps
// -----------------------------------------------------------------------------

function overlappingFeatures(window, featureNames) {

    var hits = [];

    _.each(featureNames, function(featureName) {

        var segments = featureSegments[featureName];

        var overlaps = _.some(segments, function(segment) {

            return rangesOverlap(
                window.refStart,
                window.refEnd,
                segment.refStart,
                segment.refEnd
            );
        });

        if(overlaps) {
            hits.push(featureName);
        }
    });

    return hits;
}


// -----------------------------------------------------------------------------
// Export summary
// -----------------------------------------------------------------------------

ensureParentDir(OUTPUT_PATH);

var writer =
    new BufferedWriter(
        new FileWriter(OUTPUT_PATH)
    );

try {

    writeLine(writer, [
        "variation",
        "short_name",
        "codon_start",
        "codon_end",
        "ref_start",
        "ref_end",
        "exon",
        "protein_region",
        "receptor_context",
        "structural_context"
    ]);

    _.each(GLYCO_WINDOWS, function(window) {

        var exonHits =
            overlappingFeatures(
                window,
                FEATURE_GROUPS.exon
            );

        var proteinHits =
            overlappingFeatures(
                window,
                FEATURE_GROUPS.protein_region
            );

        var receptorHits =
            overlappingFeatures(
                window,
                FEATURE_GROUPS.receptor_context
            );

        var structuralHits =
            overlappingFeatures(
                window,
                FEATURE_GROUPS.structural_context
            );

        writeLine(writer, [
            window.name,
            window.shortName,
            window.codonStart,
            window.codonEnd,
            window.refStart,
            window.refEnd,
            exonHits.join(","),
            proteinHits.join(","),
            receptorHits.join(","),
            structuralHits.join(",")
        ]);

        glue.logInfo(
            window.shortName +
            ": exon=[" + exonHits.join(",") +
            "]; protein=[" + proteinHits.join(",") +
            "]; receptor=[" + receptorHits.join(",") +
            "]; structure=[" + structuralHits.join(",") + "]"
        );
    });

} finally {

    writer.close();
}


glue.logInfo(
    "Wrote glycosylation context table to " +
    OUTPUT_PATH
);