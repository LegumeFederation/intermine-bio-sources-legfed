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

/**
 * Encapsulates a single tab-delimited SNP Marker file record.
 * <pre>
 * ID	DesignSequence	Alleles	Source	BeadType	StepDescription	3702	3885	...
 * 1_0002	TAC..[A/G]..AGA	A/G	COPA	0	List 1 inside	AT1G53750.1	Phvul.010G122200	...
 * </pre>
 * Comparator is based on marker ID.
 *
 * @author Sam Hokin, NCGR
 */
public class SNPMarkerRecord implements Comparable {
    
    String marker;
    String designSequence;
    String alleles;
    String source;
    Integer beadType;
    String stepDescription;
    String[] associatedGenes;

    /**
     * Instantiate from a line from a LinkageGroup file. Do nothing if it's a header or comment line.
     */
    public SNPMarkerRecord(String line) {

        try {
            // parse line
            String[] parts = line.split("\t");
            marker = parts[0];
            designSequence = parts[1];
            if (parts.length>2 && parts[2].length()>0) alleles = parts[2];
            if (parts.length>3 && parts[3].length()>0) source = parts[3];
            if (parts.length>4 && parts[4].length()>0) beadType = new Integer(parts[4]);
            if (parts.length>5 && parts[5].length()>0) stepDescription = parts[5];
            if (parts.length>6) {
                associatedGenes = new String[parts.length-6];
                for (int i=6; i<parts.length; i++) {
                    if (parts[i].length()>0) associatedGenes[i-6] = parts[i];
                }
            }
        } catch (Exception ex) {
            throw new RuntimeException("Error parsing genetic map file line:\n"+
                                       ex.toString()+"\n"+
                                       "marker="+marker+"\n"+
                                       "designSequence="+designSequence+"\n"+
                                       "alleles="+alleles+"\n"+
                                       "source="+source+"\n"+
                                       "beadType="+beadType+"\n"+
                                       "stepDescription="+stepDescription+"\n"+
                                       "associatedGenes="+associatedGenes);
        }
        
    }

    /**
     * For sorting by linkage group, position, marker.
     */
    public int compareTo(Object o) {
        SNPMarkerRecord that = (SNPMarkerRecord)o;
        return this.marker.compareTo(that.marker);
    }

}
