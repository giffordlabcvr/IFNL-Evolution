function calculateEntropy(refsequence, alignment, source, feature) {

    var alignmentMode = "alignment/" + alignment;
    var whereClause = "sequence.source.name = '" + source + "'";
    var countsByCodon = {};

    glue.inMode(alignmentMode, function() {

        var frequencyResult = glue.command([
            "amino-acid",
            "frequency",
            "-w", whereClause,
            "-r", refsequence,
            "-f", feature
        ]);

        var resultTable =
            frequencyResult["alignmentAminoAcidFrequencyResult"];

        if (resultTable == null || resultTable["row"] == null) {
            throw new Error(
                "The amino-acid frequency command returned no result rows"
            );
        }

        _.each(resultTable["row"], function(row) {

            /*
             * Current GLUE result layout:
             * value[1] = labelled codon
             * value[2] = amino-acid character
             * value[3] = member count
             */
            var codon = Number(row["value"][1]);
            var aminoAcid = String(row["value"][2]);
            var count = Number(row["value"][3]);

            if (isNaN(codon)) {
                throw new Error(
                    "Invalid codon value in amino-acid frequency result: " +
                    row["value"][1]
                );
            }

            if (isNaN(count)) {
                throw new Error(
                    "Invalid count at codon " + codon +
                    " for amino acid " + aminoAcid
                );
            }

            if (countsByCodon[codon] == null) {
                countsByCodon[codon] = [];
            }

            countsByCodon[codon].push({
                aminoAcid: aminoAcid,
                count: count
            });
        });
    });

    var entropyResults = [];

    _.each(countsByCodon, function(aminoAcidCounts, codon) {

        var total = 0;

        _.each(aminoAcidCounts, function(aminoAcidCount) {
            total += aminoAcidCount.count;
        });

        if (total <= 0) {
            return;
        }

        var entropy = 0;

        _.each(aminoAcidCounts, function(aminoAcidCount) {

            if (aminoAcidCount.count > 0) {
                var frequency = aminoAcidCount.count / total;
                entropy -= frequency * Math.log(frequency);
            }
        });

        entropyResults.push({
            codon: Number(codon),
            entropy: Number(entropy.toFixed(3))
        });
    });

    entropyResults.sort(function(result1, result2) {
        return result1.codon - result2.codon;
    });

    return entropyResults;
}

