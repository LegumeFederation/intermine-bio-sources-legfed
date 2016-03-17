# intermine_legfed
Legume Federation data source; install as intermine/bio/sources/legfed.

This data source contains several processors:

**SequenceProcessor.java** - a spin of Kim Rutherford's chado SequenceProcessor to support the fields present in the LIS chado database. It is heavily modified from the original under chado-db. It handles genomic, not genetic relationships, and works primarily with the feature, featureloc and feature_relationship tables.

**GeneticProcessor.java** - stores the various genetic data and relationships from the LIS chado database. These come from featuremap, featurepos, featureloc, feature_relationship and, of course, feature.

**CMapProcessor.java** - stores the genetic marker and QTL data from a CMap tab-delimited file along with a GFF3 export and a QTL-marker flat file from David Grant at Soybase. Currently hardcoded for soybean.

**SyntenyProcessor.java** - stores the synteny regions from a GFF dump with records formatted by DAGchainer. Currently hardcoded for Pv and Gm.

**GOProcessor.java** - parses the GO:terms out of chado gene.description field and associates them with the gene. Run the go data source to fill out the GO terms.

And some handy utility classes:

**ChadoFeature.java** -  incorporates the fields from the chado feature table into a single object for convenience, along with handy methods for populating Items.

**CMapRecord.java** - encapsulates a single tab-delimited CMap file record.

**GFFRecord.java** - encapsulates a single GFF3 record, parsing out the standard, Ensembl and DAGchainer GFF3 attributes.

**PubMedSearch.java** - performs a search of NCBI PubMed on journal, year and authors. This is used to get the PubMed ID from the information provided in the chado database for genetic maps.

A few chado-db classes are changed:

**ChadoDBConverter.java** - added setDagChainerFile() and getDagChainerFile(), used by SyntenyProcessor (it's a GFF3 file with DAGchainer annotation).

