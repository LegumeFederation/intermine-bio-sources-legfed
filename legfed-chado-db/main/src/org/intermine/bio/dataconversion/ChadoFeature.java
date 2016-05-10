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
     * Populate the attributes of a SequenceFeature Item with this feature's data; related sequence Item must be passed in
     */
    public void populateSequenceFeature(Item sequenceFeature, Item sequence, Item organism) {
        populateBioEntity(sequenceFeature, organism);
        sequenceFeature.setAttribute("length", String.valueOf(seqlen));
        if (residues!=null) {
            // populate the Sequence
            sequence.setAttribute("residues", residues);
            sequence.setAttribute("length", String.valueOf(seqlen));
            if (md5checksum!=null) sequence.setAttribute("md5checksum", md5checksum);
            sequenceFeature.setReference("sequence", sequence);
        }
    }

    /**
     * Embarrassing hack: return the NCBI Taxon ID for a Legume chado organism_id. Yes, I'm deeply ashamed.
     */
    public int getTaxonId(int orgId) {
        Map<Integer,Integer> orgToTaxon = new HashMap<Integer,Integer>();
        orgToTaxon.put(new Integer(1556), new Integer(13333));
        orgToTaxon.put(new Integer(31), new Integer(31));
        orgToTaxon.put(new Integer(6), new Integer(3702));
        orgToTaxon.put(new Integer(1571), new Integer(130453));
        orgToTaxon.put(new Integer(13), new Integer(3818));
        orgToTaxon.put(new Integer(1572), new Integer(130454));
        orgToTaxon.put(new Integer(14), new Integer(3821));
        orgToTaxon.put(new Integer(28), new Integer(53854));
        orgToTaxon.put(new Integer(15), new Integer(3827));
        orgToTaxon.put(new Integer(1570), new Integer(3398));
        orgToTaxon.put(new Integer(16), new Integer(3847));
        orgToTaxon.put(new Integer(17), new Integer(3864));
        orgToTaxon.put(new Integer(18), new Integer(34305));
        orgToTaxon.put(new Integer(29), new Integer(3870));
        orgToTaxon.put(new Integer(19), new Integer(3871));
        orgToTaxon.put(new Integer(20), new Integer(3879));
        orgToTaxon.put(new Integer(21), new Integer(3880));
        orgToTaxon.put(new Integer(9), new Integer(4530));
        orgToTaxon.put(new Integer(23), new Integer(3886));
        orgToTaxon.put(new Integer(24), new Integer(3885));
        orgToTaxon.put(new Integer(22), new Integer(3888));
        orgToTaxon.put(new Integer(1557), new Integer(3760));
        orgToTaxon.put(new Integer(1558), new Integer(4081));
        orgToTaxon.put(new Integer(30), new Integer(57577));
        orgToTaxon.put(new Integer(1584), new Integer(3899));
        orgToTaxon.put(new Integer(25), new Integer(3906));
        orgToTaxon.put(new Integer(26), new Integer(157791));
        orgToTaxon.put(new Integer(27), new Integer(3917));
        orgToTaxon.put(new Integer(1559), new Integer(29760));
        orgToTaxon.put(new Integer(1560), new Integer(4577));
        Integer taxId = orgToTaxon.get(new Integer(orgId));
        if (taxId!=null) {
            return taxId.intValue();
        } else {
            return 0;
        }
    }
        
}
        
