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
 * Encapsulates a single tab-delimited GT file record.
 * Comparator is based on genetic marker name.
 *
 * @author Sam Hokin, NCGR
 */
public class GTRecord implements Comparable {
    
    public String marker;
    public String[] values;

    /**
     * Instantiate from a record from a genotype file. Assume it's not a header line. (They may be a bit inconsistent.)
     */
    public GTRecord(String record) {
        
        // parse record
        String[] parts = record.split("\t");
        marker = parts[0];
        int num = parts.length - 1; // number of markers
        values = new String[num];
        for (int i=0; i<num; i++) {
            values[i] = parts[i+1];
        }

    }

    /**
     * For sorting by marker
     */
    public int compareTo(Object o) {
        GTRecord that = (GTRecord)o;
        return this.marker.compareTo(that.marker);
    }

}
