# intermine_legfed
Legume Federation data source; install as intermine/bio/sources/legfed.

This data source contains several processors of chado data:

**SequenceProcessor.java** - a spin of Kim Rutherford's chado SequenceProcessor to support the fields present in the LIS chado database. It is heavily modified from the original under chado-db. It handles genomic, not genetic relationships, and works primarily with the feature, featureloc and feature_relationship tables.

**GeneticProcessor.java** - stores the various genetic data and relationships from the LIS chado database. These come from featuremap, featurepos, featureloc, feature_relationship and, of course, feature.

**GOProcessor.java** - parses the GO:terms out of the chado.gene.description field and associates them with the gene. Run the go data source to fill out the GO terms.

**FeaturePropProcessor.java** - loads various properties of features from chado.featureprop into class attributes.

And processors of data from flat files:

**CMapFileProcessor.java** - stores the genetic markers, linkage groups and QTLs from a CMap tab-delimited file export.

**GeneticMarkerGFFProcessor.java** - stores genetic marker genomic data from a GFF3 file.

**QTLMarkerFileProcessor.java** - stores QTL-genetic marker relations from a tab-delimited file. Columns are QTLid, QTL name, marker name.

**SyntenyProcessor.java** - stores the synteny regions from a GFF dump with records formatted by DAGchainer.

Plus some handy utility classes:

**ChadoFeature.java** -  incorporates the fields from the chado feature table into a single object for convenience, along with handy methods for populating Items.

**CMapRecord.java** - encapsulates a single tab-delimited CMap file record.

**GFFRecord.java** - encapsulates a single GFF3 record, parsing out the standard, Ensembl and DAGchainer GFF3 attributes.

**PubMedSearch.java** - performs a search of NCBI PubMed on journal, year and authors. This is used to get the PubMed ID from the information provided in the chado database for genetic maps.

A few chado-db classes are changed to allow more input from project.xml:

**ChadoDBConverter.java** - added setDagChainerFile() and getDagChainerFile(), used by SyntenyProcessor (it's a GFF3 file with DAGchainer annotation).

And finally we have the post-processor for these data:

**LegfedPostProcess.java** - finds QTLs whose associated genetic markers span genes and stores the results in the QTL.spannedGenes and Gene.spanningQTLs collections.
