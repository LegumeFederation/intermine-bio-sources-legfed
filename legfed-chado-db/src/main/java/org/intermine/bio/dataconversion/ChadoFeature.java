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

import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Timestamp;

import java.util.Map;
import java.util.HashMap;

import org.intermine.xml.full.Item;

/**
 * Incorporates the fields from the chado feature table into a single object for convenience, along with handy methods.
 * equals method is implemented which compares feature_id values.
 *
 * @author Sam Hokin
 */

public class ChadoFeature {

    // chado.feature fields
    public int feature_id;
    public int dbxref_id;
    public int organism_id;
    public String name;
    public String uniquename; 
    public String residues;
    public int seqlen;
    public String md5checksum;
    public int type_id;
    public boolean is_analysis;
    public boolean is_obsolete;
    public Timestamp timeaccessioned;
    public Timestamp timelastmodified;

    /**
     * Construct a ChadoFeature from a populated ResultSet
     */
    public ChadoFeature(ResultSet rs) throws SQLException {
        feature_id = rs.getInt("feature_id");
        dbxref_id = rs.getInt("dbxref_id");
        organism_id = rs.getInt("organism_id");
        name = rs.getString("name");
        uniquename = rs.getString("uniquename"); // never null and length>0 always
        residues = rs.getString("residues");
        seqlen = rs.getInt("seqlen");
        md5checksum = rs.getString("md5checksum");
        type_id = rs.getInt("type_id");
        is_analysis = rs.getBoolean("is_analysis");
        is_obsolete = rs.getBoolean("is_obsolete");
        timeaccessioned = rs.getTimestamp("timeaccessioned");
        timelastmodified = rs.getTimestamp("timelastmodified");
    }

    /**
     * Two features are equal if they have the same feature_id
     */
    public boolean equals(ChadoFeature that) {
        return this.feature_id==that.feature_id;
    }
   

    /**
     * Populate the attributes of a BioEntity Item with this feature's data.
     */
    public void populateBioEntity(Item bioEntity, Item organism) {
        bioEntity.setAttribute("chadoFeatureId", String.valueOf(feature_id));
        bioEntity.setAttribute("primaryIdentifier", uniquename);
        if (name!=null && name.length()>0) {
            bioEntity.setAttribute("secondaryIdentifier", name);
        }
        bioEntity.setReference("organism", organism);
    }

    /**
     * Populate the attributes of a SequenceFeature Item with this feature's data; related Sequence Item must be passed in
     */
    public void populateSequenceFeature(Item sequenceFeature, Item sequence, Item organism) {
        populateBioEntity(sequenceFeature, organism);
        if (seqlen>0) sequenceFeature.setAttribute("length", String.valueOf(seqlen));
        if (residues!=null) {
            // populate the Sequence
            sequence.setAttribute("residues", residues);
            if (seqlen>0) sequence.setAttribute("length", String.valueOf(seqlen));
            if (md5checksum!=null) sequence.setAttribute("md5checksum", md5checksum);
            sequenceFeature.setReference("sequence", sequence);
        }
    }
        
}
        
