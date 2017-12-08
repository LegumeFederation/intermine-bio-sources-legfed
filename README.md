# intermine_legfed
Legume Federation data sources; install as intermine/bio/sources/legfed.

This data source directory contains several data sources:

#### legfed-chado-db
A custom version of FlyMine's chado-db source written by Kim Rutherford.

**SequenceProcessor.java** - a modification to support the fields present in the LIS chado database. It still handles strictly genomic, not genetic items, and works primarily with the feature, featureloc and feature_relationship tables.

**GeneticProcessor.java** - stores genetic data and relationships from the LIS chado database. These come from featuremap, featurepos, featureloc, feature_relationship and, of course, feature. Populates GeneFamily, GeneticMap, LinkageGroup, GeneticMarker and QTL objects.

**FeaturePropProcessor.java** - loads various properties of features from chado.featureprop. Doesn't create any new objects, but fills in many attributes.

#### legfed-cmap-file

**CMapFileConverter.java** - reads and stores genetic markers, linkage groups and QTLs from a CMap tab-delimited file export. Extends BioFileConverter. Creates GeneticMarker, LinkageGroup and QTL objects if they don't already exist.

#### legfed-geneticmarker-gff

**GeneticMarkerGFFConverter.java** - reads and stores genetic marker genomic data from a GFF3 file. Extends BioFileConverter. Creates GeneticMarker objects if they don't already exist.

#### legfed-qtlmarker-file

**QTLMarkerFileConverter.java** - reads and stores QTL-genetic marker relations from a tab-delimited file. Columns are QTLid, QTL name, marker name. Extends BioFileConverter. Creates GeneticMarker and QTL objects if they don't already exist.

#### legfed-synteny-gff

**SyntenyGFFConverter.java** - reads and stores synteny blocks and synteny regions from a GFF dump with records formatted by DAGchainer. Creates SyntenyBlock, SyntenicRegion objects.

#### legfed-expression

**ExpressionFileConverter.java** - reads a tab-delimited expression file of custom format into ExpressionSource, ExpressionSample and ExpressionValue objects.
