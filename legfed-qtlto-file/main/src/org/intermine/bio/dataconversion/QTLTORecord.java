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
 * #QTL       TOTerms
 * QRk-vu11.1 TO:0002704,TO:0002706
 *
 * @author Sam Hokin, NCGR
 */
public class QTLTORecord implements Comparable {

    private static final Logger LOG = Logger.getLogger(QTLTORecord.class);

    String qtlName;
    String[] toTerms;

    /**
     * Instantiate from a line from a QTLTO file. Do nothing if it's a header line.
     * This allows extra ignored columns to exist.
     */
    public QTLTORecord(String line) {
        if (!line.startsWith("#")) {
	    String[] parts = line.split("\t");
	    if (parts.length>=2) {
		qtlName = parts[0];
                toTerms = parts[1].split(",");
                // clear toTerms if empty
                boolean ok = true;
                for (int i=0; i<toTerms.length; i++) {
                    if (toTerms[i]==null || toTerms[i].trim().length()==0) ok = false;
                }
                if (!ok) toTerms = new String[0];
            }
	}
    }

    /**
     * For sorting; map name, or if equal, feature start (int)
     */
    public int compareTo(Object o) {
        QTLTORecord that = (QTLTORecord)o;
        return this.qtlName.compareTo(that.qtlName);
    }

}
