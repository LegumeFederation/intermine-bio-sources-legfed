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
 * Encapsulates a single tab-delimited QTLMarker file record.
 * Comparator is based on QTL name and marker name.
 *
 * File format:
 * QTLID[tab]QTLName[tab]LocusName
 *
 * @author Sam Hokin, NCGR
 */
public class QTLMarkerRecord implements Comparable {
    
    String qtlId;
    String qtlName;
    String markerName;

    /**
     * Instantiate from a line from a QTLMarker file. Do nothing if it's a header line.
     */
    public QTLMarkerRecord(String line) {

        if (!line.startsWith("QTLID")) {
            // parse line
            String[] parts = line.split("\t");
            qtlId = parts[0];
            qtlName = parts[1];
            markerName = parts[2];
        }

    }

    /**
     * For sorting; map name, or if equal, feature start (int)
     */
    public int compareTo(Object o) {
        QTLMarkerRecord that = (QTLMarkerRecord)o;
        if (this.qtlName.equals(that.qtlName)) {
            return this.markerName.compareTo(that.markerName);
        } else {
            return this.qtlName.compareTo(that.qtlName);
        }
    }

}
