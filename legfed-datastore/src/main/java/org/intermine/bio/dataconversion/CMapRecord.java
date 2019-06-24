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
 *
 * map_acc                 map_name map_start map_stop feature_acc                         feature_name feature_aliases feature_start feature_stop feature_type_acc is_landmark
 * PvCookUCDavis2009_Pv09  Pv09     0         66.1     PvCookUCDavis2009_Pv09_Pv_TOG913042 Pv_TOG913042 TOG913042       66.1          66.1         SNP              0
 *
 * @author Sam Hokin, NCGR
 */
public class CMapRecord {
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
        if (!line.startsWith("#") && !line.startsWith("map_acc")) {
            try {
                // parse line
                String[] parts = line.split("\t");
                map_acc = parts[0];
                map_name = parts[1];
                map_start = Double.parseDouble(parts[2]);
                map_stop = Double.parseDouble(parts[3]);
                feature_acc = parts[4].replace("\"", "");
                feature_name = parts[5].replace("\"", "");
                feature_aliases = parts[6];
                feature_start = Double.parseDouble(parts[7]);
                feature_stop = Double.parseDouble(parts[8]);
                feature_type_acc = parts[9];
                if (parts.length>10) is_landmark = Integer.parseInt(parts[10])==1;
            } catch (Exception ex) {
                throw new RuntimeException("Error parsing CMap line:\n"+
                                           line+"\n"+
                                           "map_acc="+map_acc+"\n"+
                                           "map_name="+map_name+"\n"+
                                           "map_start="+map_start+"\n"+
                                           "map_stop="+map_stop+"\n"+
                                           "feature_acc="+feature_acc+"\n"+
                                           "feature_name="+feature_name+"\n"+
                                           "feature_aliases="+feature_aliases+"\n"+
                                           "feature_start="+feature_start+"\n"+
                                           "feature_stop="+feature_stop+"\n"+
                                           "feature_type_acc="+feature_type_acc+"\n"+
                                           "is_landmark="+is_landmark);
            }
        }
    }

    /**
     * Return true if this record actually has data.
     */
    boolean hasData() {
        return map_acc!=null;
    }

    /**
     * Return true if record is a QTL.
     */
    boolean isQTL() {
        return feature_type_acc!=null && feature_type_acc.startsWith("QTL");
    }

    /**
     * Return true if record is a genetic marker, i.e. is not a QTL.
     */
    boolean isMarker() {
        return !isQTL();
    }

    /**
     * Return true if record is a SNP marker.
     */
    boolean isSNP() {
        return feature_type_acc!=null && feature_type_acc.equals("SNP");
    }
}
