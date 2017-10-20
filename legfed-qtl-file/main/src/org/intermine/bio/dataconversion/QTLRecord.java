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
 * Encapsulates a single tab-delimited QTL file record.
 *
 * QTLName Parent_1 Parent_2 TraitName Journal Year Volume Page Title
 *
 * @author Sam Hokin, NCGR
 */
public class QTLRecord implements Comparable {

    private static final Logger LOG = Logger.getLogger(QTLRecord.class);

    String qtlName;
    String parent1;
    String parent2;
    String traitName;
    String pubJournal;
    String pubYear;
    String pubVolume;
    String pubPages;
    String pubTitle;
    int pubPMID;

    /**
     * Instantiate from a line from a QTL file. Do nothing if it's a comment.
     */
    public QTLRecord(String line) {
        if (!line.startsWith("#")) {
	    String[] parts = line.split("\t");
            qtlName = parts[0];
            parent1 = parts[1];
            parent2 = parts[2];
            traitName = parts[3];
            pubJournal = parts[4];
            pubYear = parts[5];
            pubVolume = parts[6];
            pubPages = parts[7];
            pubTitle = parts[8];
            if (parts.length>9) {
                try {
                    pubPMID = Integer.parseInt(parts[9]);
                } catch (Exception e) {
                    // do nothing, probably some stray extra column data
                }
            }
	}
    }

    /**
     * For alpha sorting on QTL name
     */
    public int compareTo(Object o) {
        QTLRecord that = (QTLRecord)o;
        return this.qtlName.compareTo(that.qtlName);
    }

}
