package org.intermine.bio.dataconversion;

import org.apache.log4j.Logger;

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
 * Only Marker and QTL are required; Trait and comma-separated TOTerms are optional.
 *
 * #Marker       QTL        Trait          TOTerms
 * 1_0531        QRk-vu11.1 RKN resistance TO:0002704,TO:0002706
 *
 * @author Sam Hokin, NCGR
 */
public class QTLMarkerRecord implements Comparable {

    private static final Logger LOG = Logger.getLogger(QTLMarkerRecord.class);

    String markerName;
    String qtlName;
    String trait;
    String[] toTerms;

    /**
     * Instantiate from a line from a QTLMarker file. Do nothing if it's a header line.
     */
    public QTLMarkerRecord(String line) {

        if (!line.startsWith("#")) {
	    String[] parts = line.split("\t");
	    if (parts.length>=2) {
		markerName = parts[0];
		qtlName = parts[1];
		if (parts.length>=3) {
		    trait = parts[2];
		    if (parts.length>=4)
			toTerms = parts[3].split(",");
		}
	    } else {
		LOG.info("Ignored line:"+line);
	    }
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
