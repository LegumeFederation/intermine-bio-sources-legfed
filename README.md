# intermine_legfed
Legume Federation data source; install as intermine/bio/sources/legfed.

This data source directory contains several data sources:

**/legfed-chado-db:**

**SequenceProcessor.java** - a spin of Kim Rutherford's chado SequenceProcessor to support the fields present in the LIS chado database. It is heavily modified from the original under chado-db. It handles genomic, not genetic relationships, and works primarily with the feature, featureloc and feature_relationship tables.

**GeneticProcessor.java** - stores the various genetic data and relationships from the LIS chado database. These come from featuremap, featurepos, featureloc, feature_relationship and, of course, feature.

**GOProcessor.java** - parses the GO:terms out of the chado.gene.description field and associates them with the gene. Run the go data source to fill out the GO terms.

**FeaturePropProcessor.java** - loads various properties of features from chado.featureprop into class attributes.

**/legfed-cmap-file:**

**CMapFileConverter.java** - stores the genetic markers, linkage groups and QTLs from a CMap tab-delimited file export.

**/legfed-geneticmarker-gff:**

**GeneticMarkerGFFConverter.java** - stores genetic marker genomic data from a GFF3 file.

**/legfed-qtlmarker-file:**

**QTLMarkerFileConverter.java** - stores QTL-genetic marker relations from a tab-delimited file. Columns are QTLid, QTL name, marker name.

**/legfed-synteny-gff:**

**SyntenyGFFConverter.java** - stores the synteny regions from a GFF dump with records formatted by DAGchainer.

