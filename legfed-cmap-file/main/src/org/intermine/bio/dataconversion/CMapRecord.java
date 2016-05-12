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
 * Encapsulates a single tab-delimited CMap file record.
 * Comparator is based on map_name, feature_start and feature_stop.
 *
 * @author Sam Hokin, NCGR
 */
public class CMapRecord implements Comparable {
    
    String map_acc;
    String map_name;
    double map_start;
    double map_stop;
    String feature_acc;
    String feature_name;
    String feature_aliases;
    double feature_start;
    double feature_stop;
    String feature_type_acc;
    boolean is_landmark;

    /**
     * Instantiate from a line from a CMap file. Do nothing if it's a header line.
     */
    public CMapRecord(String line) {

        if (!line.startsWith("map_acc")) {
            // parse line
            String[] parts = line.split("\t");
            map_acc = parts[0];
            map_name = parts[1];
            map_start = Double.parseDouble(parts[2]);
            map_stop = Double.parseDouble(parts[3]);
            feature_acc = parts[4];
            feature_name = parts[5];
            feature_aliases = parts[6];
            feature_start = Double.parseDouble(parts[7]);
            feature_stop = Double.parseDouble(parts[8]);
            feature_type_acc = parts[9];
            is_landmark = Integer.parseInt(parts[10])==1;
        }

    }

    /**
     * For sorting; map name, or if equal, feature start (int)
     */
    public int compareTo(Object o) {
        CMapRecord that = (CMapRecord)o;
        if (this.map_name.equals(that.map_name)) {
            return (int)(this.feature_start - that.feature_start);
        } else {
            return this.map_name.compareTo(that.map_name);
        }
    }

    /**
     * Return true if record is a QTL
     */
    boolean isQTL() {
        return feature_type_acc!=null && feature_type_acc.startsWith("QTL");
    }

    /**
     * Return true if record is a genetic marker
     */
    boolean isMarker() {
        return !isQTL();
    }

}
