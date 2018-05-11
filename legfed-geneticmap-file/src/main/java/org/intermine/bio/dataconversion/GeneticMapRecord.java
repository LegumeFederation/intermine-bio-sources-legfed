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
 * Encapsulates a single tab-delimited LinkageGroup file record of the form:
 * <pre>
 * Marker	LG	Type	Pos (cM)	QTL	Traits
 * </pre>
 * Comparator is based on Linkage group, Position and Marker.
 * QTL and Traits associated are individually optional.
 *
 * @author Sam Hokin, NCGR
 */
public class GeneticMapRecord implements Comparable {
    
    String marker;
    int lg;
    String type;
    double position;
    String qtl;    // optional
    String traits; // optional

    /**
     * Instantiate from a line from a LinkageGroup file. Do nothing if it's a header or comment line.
     */
    public GeneticMapRecord(String line) {

        if (!line.startsWith("#") && !line.startsWith("Marker")) {
            try {
                // parse line
                String[] parts = line.split("\t");
                marker = parts[0];
                lg = Integer.parseInt(parts[1]);
                type = parts[2];
                position = Double.parseDouble(parts[3]);
                if (parts.length>4) qtl = parts[4];
                if (parts.length>5) traits = parts[5];
            } catch (Exception ex) {
                throw new RuntimeException("Error parsing genetic map file line:\n"+
                                           ex.toString()+"\n"+
                                           line+"\n"+
                                           "marker="+marker+"\n"+
                                           "lg="+lg+"\n"+
                                           "type="+type+"\n"+
                                           "position="+position+"\n"+
                                           "qtl="+qtl+"\n"+
                                           "traits="+traits);
            }
        }
        
    }

    /**
     * For sorting by linkage group, position, marker.
     */
    public int compareTo(Object o) {
        GeneticMapRecord that = (GeneticMapRecord)o;
        if (this.lg!=that.lg) {
            return this.lg - that.lg;
        } else if (this.position!=that.position) {
            return (int) ((this.position-that.position)*100);
        } else {
            return this.marker.compareTo(that.marker);
        }
    }

}
