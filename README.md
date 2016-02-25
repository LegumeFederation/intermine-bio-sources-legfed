# intermine_legfed
Legume Federation data source; install as intermine/bio/sources/legfed.

This data source contains several processors as well as utility classes:

**SequenceProcessor.java** - a spin of Kim Rutherford's chado SequenceProcessor to support the fields present in the LIS chado database. It is heavily modified from the original under chado-db. It handles genomic, not genetic relationships, and works primarily with the feature, featureloc and feature_relationship tables.

**GeneticProcessor.java** - stores the various genetic data and relationships from the LIS chado database. These come from featuremap, featurepos, featureloc, feature_relationship and, of course, feature.

**CMapProcessor.java** - stores the genetic marker and QTL data from a CMap tab-delimited file along with a GFF3 export. Currently hardcoded to support a couple files exported from Soybase.

**ChadoFeature.java** -  incorporates the fields from the chado feature table into a single object for convenience, along with handy methods.

**CMapRecord.java** - encapsulates a single tab-delimited CMap file record.

**GFFRecord.java** - encapsulates a single GFF3 record, parsing out the standard Ensembl GFF3 attributes.

**PubMedSearch.java** - performs a search of NCBI PubMed on journal, year and authors. This is used to get the PubMed ID from the information provided in the chado database for genetic maps.

The other classes are unchanged from chado-db (I think).

