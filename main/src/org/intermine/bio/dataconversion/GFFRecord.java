package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2015-2016 NCGR
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */

import java.util.Comparator;

import org.intermine.xml.full.Item;

/**
 * Encapsulates a single GFF3 record, parsing out extra attributes from the attributes field:
 *   Ensembl
 *   DAGchainer
 *
 * Native comparator is based on seqid, strand, start and end.
 * Alternative comparators are provided for alternative sorting.
 *
 * Missing string attributes are null.
 *
 * @author Sam Hokin, NCGR
 */
public class GFFRecord implements Comparable {

    // GFF3 fields
    public String seqid;
    public String source;
    public String type;
    public int start;
    public int end;
    public String score;
    public char strand;
    public char phase;
    public String attributes;

    // standard attributes
    String attributeID;
    String attributeName;
    String attributeParent;

    // Ensembl attributes
    String attributeBiotype;
    String attributeDescription;
    String attributeDbxref;
    String attributeOntologyTerm;
    int attributeVersion;
    int attributeConstitutive;
    int attributeRank;

    // DAGchainer attributes
    String attributeTarget;
    double attributeMedianKs;
    
    /**
     * Instantiate from a line from a GFF file. Do nothing if it's a comment.
     */
    public GFFRecord(String line) {
        if (!line.startsWith("#")) {
            // GBrowse compatibility
            line = line.replace("NAME","Name");
            line = line.replace("PARENT","Parent");
            // parse out fields
            String[] chunks = line.split("\t");
            if (chunks.length==9) {
                // main fields
                seqid = chunks[0];
                source = chunks[1];
                type = chunks[2];
                start = Integer.parseInt(chunks[3]);
                end = Integer.parseInt(chunks[4]);
                score = chunks[5];
                strand = chunks[6].charAt(0);
                phase = chunks[7].charAt(0);
                attributes = chunks[8];
                // parse out known attributes
                parseAttributes();
            }
        }
    }

    /**
     * Instantiate from a full set of values
     */
    public GFFRecord(String seqid, String source, String type, int start, int end, String score, char strand, char phase, String attributes) {
        this.seqid = seqid;
        this.source = source;
        this.type = type;
        this.start = start;
        this.end = end;
        this.score = score;
        this.strand = strand;
        this.phase = phase;
        this.attributes = attributes;
        parseAttributes();
    }

    /**
     * Native comparator: compare on seqid, strand, start and end.
     */
    public int compareTo(Object o) {
        GFFRecord that = (GFFRecord)o;
        // alpha seqid
        if (!this.seqid.equals(that.seqid)) return this.seqid.compareTo(that.seqid);
        // + strand ahead of - strand
        if (this.strand!='.' && that.strand!='.' && this.strand!=that.strand) {
            if (this.strand=='+') return 1; else return -1;
        }
        // left start before right start
        if (this.start!=that.start) return this.start - that.start;
        // left end before right end
        if (this.end!=that.end) return this.end - that.end;
        // it's a tie!
        return 0;
    }
    
    /**
     * Equal if both have same seqid, strand, start and end (has to be consistent with compareTo above).
     */
    public boolean equals(Object o) {
        GFFRecord that = (GFFRecord)o;
        return this.seqid.equals(that.seqid) && this.strand==that.strand && this.start==that.start && this.end==that.end;
    }

    /**
     * Output a tab-delimited line containing this instance's values
     */
    public String toString() {
        return seqid+'\t'+source+'\t'+type+'\t'+start+'\t'+end+'\t'+score+'\t'+strand+'\t'+phase+'\t'+attributes;
    }

    /**
     * Parse out known attributes; if Name is missing, set it to ID and vice-versa.
     */
    void parseAttributes() {
        String[] chunks = attributes.split(";");
        for (int i=0; i<chunks.length; i++) {
            // standard
            if (chunks[i].toLowerCase().startsWith("id=")) attributeID = getAttributeValue(chunks[i]);
            if (chunks[i].toLowerCase().startsWith("name=")) attributeName = getAttributeValue(chunks[i]);
            if (chunks[i].toLowerCase().startsWith("parent=")) attributeParent = getAttributeValue(chunks[i]);
            // Ensembl
            if (chunks[i].toLowerCase().startsWith("biotype=")) attributeBiotype = getAttributeValue(chunks[i]);
            if (chunks[i].toLowerCase().startsWith("description=")) attributeDescription = getAttributeValue(chunks[i]);
            if (chunks[i].toLowerCase().startsWith("dbxref=")) attributeDbxref = getAttributeValue(chunks[i]);
            if (chunks[i].toLowerCase().startsWith("ontology_term=")) attributeOntologyTerm = getAttributeValue(chunks[i]);
            if (chunks[i].toLowerCase().startsWith("version=")) attributeVersion = Integer.parseInt(getAttributeValue(chunks[i]));
            if (chunks[i].toLowerCase().startsWith("constitutive=")) attributeConstitutive = Integer.parseInt(getAttributeValue(chunks[i]));
            if (chunks[i].toLowerCase().startsWith("rank=")) attributeRank = Integer.parseInt(getAttributeValue(chunks[i]));
            // DAGchainer
            if (chunks[i].toLowerCase().startsWith("target=")) attributeTarget = getAttributeValue(chunks[i]);
            if (chunks[i].toLowerCase().startsWith("median_ks=")) attributeMedianKs = Double.parseDouble(getAttributeValue(chunks[i]));
        }
        // set Name to ID if missing
        if (attributeName==null && attributeID!=null) updateAttributeName(attributeID);
        // set ID to Name if missing
        if (attributeID==null && attributeName!=null) updateAttributeID(attributeName);
    }

    /**
     * Parse the value from a single attribute string
     */
    String getAttributeValue(String attributeString) {
        String[] parts = attributeString.split("=");
        return parts[1];
    }

    /**
     * Update the attributes string with a new attribute Name, or sets it if it's not already there
     */
    public void updateAttributeName(String newName) {
        String[] chunks = attributes.split(";");
        attributes = "";
        boolean nameSet = false;
        for (int i=0; i<chunks.length; i++) {
            if (chunks[i].toLowerCase().startsWith("name=")) {
                attributeName = newName;
                attributes += "Name="+newName+";";
                nameSet = true;
            } else {
                attributes += chunks[i]+";";
            }
        }
        if (!nameSet) {
            attributeName = newName;
            attributes = "Name="+newName+";"+attributes;
        }
    }

    /**
     * Update the attributes string with a new attribute ID, or sets it if it's not already there.
     */
    public void updateAttributeID(String newID) {
        String[] chunks = attributes.split(";");
        attributes = "";
        boolean idSet = false;
        for (int i=0; i<chunks.length; i++) {
            if (chunks[i].toLowerCase().startsWith("id=")) {
                attributeID = newID;
                attributes += "ID="+newID+";";
                idSet = true;
            } else {
                attributes += chunks[i]+";";
            }
        }
        if (!idSet) {
            attributeID = newID;
            attributes = "ID="+newID+";"+attributes;
        }
    }

    /**
     * Update the attributes string with a new attribute ID, or sets it if it's not already there.
     */
    public void updateAttributeID(int newID) {
        String[] chunks = attributes.split(";");
        attributes = "";
        boolean idSet = false;
        for (int i=0; i<chunks.length; i++) {
            if (chunks[i].toLowerCase().startsWith("id=")) {
                attributeID = ""+newID;
                attributes += "ID="+newID+";";
                idSet = true;
            } else {
                attributes += chunks[i]+";";
            }
        }
        if (!idSet) {
            attributeID = ""+newID;
            attributes = "ID="+newID+";"+attributes;
        }
    }

    /**
     * Calculate the gap between two sequences; GFF always has start<end
     */
    public static int gap(GFFRecord gff1, GFFRecord gff2) {
        if (gff1.end<gff2.start) {
            return gff2.start - gff1.end;  // 1 left of 2
        } else {
            return gff1.start - gff2.end;  // 2 left of 1
        }
    }

    /**
     * Compare based on ID attribute, else default
     */
    public static final Comparator<GFFRecord> IDComparator = new Comparator<GFFRecord>() {
            @Override
            public int compare(GFFRecord gff1, GFFRecord gff2) {
                if (gff1.attributeID!=null && gff2.attributeID!=null) {
                    return gff1.attributeID.compareTo(gff2.attributeID);
                } else {
                    return gff1.compareTo(gff2);
                }
            }
        };
	
    /**
     * Compare based on attributeName and score, assuming score can be parsed to an int; DESCENDING on score so first in remains
     */
    public static final Comparator<GFFRecord> ScoreComparator = new Comparator<GFFRecord>() {
            @Override
            public int compare(GFFRecord gff1, GFFRecord gff2) {
                if (gff1.attributeName!=null && gff1.score!=null && gff2.attributeName!=null && gff2.score!=null) {
                    if (!gff1.attributeName.equals(gff2.attributeName)) {
                        return gff1.attributeName.compareTo(gff2.attributeName);
                    } else {
                        double score1 = Double.parseDouble(gff1.score);
                        double score2 = Double.parseDouble(gff2.score);
                        return (int)(score2-score1); // DESCENDING
                    }
                } else {
                    return gff1.compareTo(gff2);
                }
            }
        };

    /**
     * Compare based on attributeName; else default
     */
    public static final Comparator<GFFRecord> NameComparator = new Comparator<GFFRecord>() {
            @Override
            public int compare(GFFRecord gff1, GFFRecord gff2) {
                if (gff1.attributeName!=null && gff2.attributeName!=null) {
                    return gff1.attributeName.compareTo(gff2.attributeName);
                } else {
                    return gff1.compareTo(gff2);
                }
            }
        };

    /**
     * Return the DAGchainer target strand from this instance
     */
    public char getTargetStrand() {
        return attributeName.charAt(attributeName.length()-1);
    }

    /**
     * Return the DAGchainer target chromosome from this instance. Return null if not a DAGchainer record.
     */
    public String getTargetChromosome() {
        if (attributeTarget!=null) {
            String[] chunks = attributeTarget.split("%20");
            return chunks[0];
        } else {
            return null;
        }
    }

    /**
     * Return the DAGchainer target sequence start from this instance. Return 0 if not a DAGchainer record.
     */
    public int getTargetStart() {
        if (attributeTarget!=null) {
            String[] chunks = attributeTarget.split("%20");
            return Integer.parseInt(chunks[1]);
        } else {
            return 0;
        }
    }

    /**
     * Return the DAGchainer target sequence end from this instance. Return 0 if not a DAGchainer record.
     */
    public int getTargetEnd() {
        if (attributeTarget!=null) {
            String[] chunks = attributeTarget.split("%20");
            return Integer.parseInt(chunks[2]);
        } else {
            return 0;
        }
    }

    /**
     * Populate the attributes of a SequenceFeature Item with this GFFRecord's data; Organism, Chromosome and ChromosomeLocation Items must be passed in
     */
    public void populateSequenceFeature(Item sequenceFeature, Item organism, Item chromosome, Item chromosomeLocation) {
        sequenceFeature.setAttribute("primaryIdentifier", attributeName);
        sequenceFeature.setAttribute("length", String.valueOf(end-start+1));
        sequenceFeature.setAttribute("score", String.valueOf(score));
        sequenceFeature.setReference("organism", organism);
        sequenceFeature.setReference("chromosome", chromosome);
        sequenceFeature.setReference("chromosomeLocation", chromosomeLocation);
        chromosomeLocation.setAttribute("start", String.valueOf(start));
        chromosomeLocation.setAttribute("end", String.valueOf(end));
        chromosomeLocation.setAttribute("strand", String.valueOf(strand));
        chromosomeLocation.setReference("feature", sequenceFeature);
        chromosomeLocation.setReference("locatedOn", chromosome);
    }

    /**
     * Populate the attributes of a SequenceFeature Item with this GFFRecord's DAGchainer data; Organism, Chromosome and ChromosomeLocation Items must be passed in.
     * Does nothing if not a DAGchainer record.
     */
    public void populateDAGchainerRegion(Item sequenceFeature, Item organism, Item chromosome, Item chromosomeLocation) {
        if (attributeTarget!=null) {
            sequenceFeature.setAttribute("primaryIdentifier", attributeName);
            sequenceFeature.setAttribute("length", String.valueOf(getTargetEnd()-getTargetStart()+1));
            sequenceFeature.setReference("organism", organism);
            sequenceFeature.setReference("chromosome", chromosome);
            sequenceFeature.setReference("chromosomeLocation", chromosomeLocation);
            chromosomeLocation.setAttribute("start", String.valueOf(getTargetStart()));
            chromosomeLocation.setAttribute("end", String.valueOf(getTargetEnd()));
            chromosomeLocation.setAttribute("strand", String.valueOf(getTargetStrand()));
            chromosomeLocation.setReference("feature", sequenceFeature);
            chromosomeLocation.setReference("locatedOn", chromosome);
        }
    }

}
