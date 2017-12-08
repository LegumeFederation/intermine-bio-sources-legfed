package org.intermine.bio.dataconversion;

/**
 * Copyright (C) 2015-2017 NCGR
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */

/**
 * Encapsulates a single tab-delimited SNP VCF file record describing a SNP marker's position on the genome.
 *
 * <pre>
 * CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
 * Vu01	74363	2_37329	A	G	999	.	DP=378
 * </pre>
 *
 * Comparator is based on ID.
 *
 * @author Sam Hokin
 */
public class SNPVCFRecord implements Comparable {
    
    // VCF fields
    String chromosome;
    int pos;
    String id;
    String ref; // usually a single character like A
    String alt; // usually a single character like G
    int qual;
    String filter;
    String info;

    /**
     * Instantiate from a line from a VCF file. Do nothing if it's a header or comment line.
     */
    public SNPVCFRecord(String line) {
        
        try {
            
            // parse line
            String[] parts = line.split("\t");
            chromosome = parts[0];
            pos = Integer.parseInt(parts[1]);
            id = parts[2];
            ref = parts[3];
            alt = parts[4];
            if (!parts[5].equals(".")) qual = Integer.parseInt(parts[5]);
            if (!parts[6].equals(".")) filter = parts[6];
            info = parts[7];
            
        } catch (Exception ex) {
            throw new RuntimeException("Error parsing VCF file line:\n"+
                                       ex.toString()+"\n"+
                                       "chromosome="+chromosome+"\n"+
                                       "pos="+pos+"\n"+
                                       "id="+id+"\n"+
                                       "ref="+ref+"\n"+
                                       "alt="+alt+"\n"+
                                       "qual="+qual+"\n"+
                                       "filter="+filter+"\n"+
                                       "info="+info);
        }
        
    }

    /**
     * For sorting by chromosome and position
     */
    public int compareTo(Object o) {
        SNPVCFRecord that = (SNPVCFRecord)o;
	if (!this.chromosome.equals(that.chromosome)) {
	    return this.chromosome.compareTo(that.chromosome);
	} else {
	    return this.pos-that.pos;
	}
    }

}
