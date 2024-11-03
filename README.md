IFNL-GLUE: Comparative Genomics of Interferon Lambda Genes
----------------------------------------------------------

Background
----------

IFNL-GLUE is a sequence-oriented resource for the comparative genomic analysis of interferon lambda (IFNL) genes, developed using the GLUE software framework.

**Interferon lambda (IFNL)** is a family of cytokines, including IFN-λ1, IFN-λ2, IFN-λ3, and IFN-λ4, that play a crucial role in mucosal immunity by activating antiviral responses primarily in epithelial cells. Unlike type I interferons, which act broadly across cell types, IFNL's activity is limited to cells expressing the specific IFN-λ receptor (IFN-λR1), found primarily in mucosal tissues such as the respiratory and gastrointestinal tracts. This focused activity allows IFNL to control infections at localized sites without widespread immune activation, making it a critical factor in immune responses to pathogens at mucosal surfaces. Variations in IFN-λ4, for example, have been associated with differences in disease susceptibility and treatment outcomes, emphasizing the importance of IFNL diversity in host-pathogen interactions.

**GLUE** is an open, integrated software toolkit designed for storing and interpreting sequence data. It supports the creation of custom projects that incorporate essential data items for comparative genomic analysis, including sequences, multiple sequence alignments, gene feature annotations, and other associated metadata. Projects are loaded into the GLUE "engine," creating a relational database that represents the semantic relationships between data items. This structure enables systematic comparative analyses and the development of robust, sequence-based resources for genomic research.

Features
--------

-   **Comprehensive Database**: Contains a curated collection of IFNL sequences with detailed metadata, supporting research into the evolutionary and functional diversity of interferon lambda genes.

-   **Genotyping and Visualization Tool**: Facilitates genotyping and visualization of submitted IFNL sequences, allowing users to analyze and compare variants across the gene family.

-   **Pre-built Multiple Sequence Alignments**: Includes pre-built alignments of IFNL genes, available for download to support customized comparative analyses.

-   **Phylogenetic Structure**: Organizes IFNL sequence data in a phylogenetically structured format, allowing users to explore evolutionary relationships and genetic diversity within the IFNL gene family.

-   **Automated Genotyping**: Uses GLUE's maximum likelihood clade assignment (MLCA) algorithm to classify IFNL sequences by subtype and lineage, supporting systematic classification across the gene family.

-   **Rich Annotations**: Features annotated reference sequences for comparative analysis, focusing on conservation, structural context, and functional implications of genetic variation.

-   **Web User Interface**: The IFNL-GLUE-WEB extension provides a web interface for browsing the IFNL-GLUE database, as well as tools for sequence analysis and visualization.

-   **Research and Clinical Applicability**: Supports both exploratory research and clinical studies on IFNL genes, enabling insights into immune response variability and potential therapeutic targets.

Usage
-----

GLUE provides an interactive command line environment for bioinformaticians to develop and use GLUE projects efficiently. This environment includes productivity-oriented features such as command completion, command history, and interactive paging through tabular data.

For detailed instructions on using IFNL-GLUE for comparative genomic analysis, refer to the [GLUE reference documentation](http://glue-tools.cvr.gla.ac.uk/).

Data Sources
------------

IFNL-GLUE relies on the following data sources:

-   [NCBI Nucleotide](https://www.ncbi.nlm.nih.gov/nuccore)
-   [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy)

Contributing
------------

We welcome contributions from the community! If you're interested in contributing to IFNL-GLUE, please review our [Contribution Guidelines](./md/CONTRIBUTING.md).

[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](./md/code_of_conduct.md)

License
-------

The project is licensed under the [GNU Affero General Public License v. 3.0](https://www.gnu.org/licenses/agpl-3.0.en.html)

Contact
-------

For questions, issues, or feedback, please contact us at gluetools@gmail.com or open an issue on the [GitHub repository](https://github.com/giffordlabcvr/IFNL-GLUE/issues).

* * * * *
