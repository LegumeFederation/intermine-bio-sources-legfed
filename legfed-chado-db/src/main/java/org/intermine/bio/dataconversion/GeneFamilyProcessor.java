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

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.TreeSet;
import java.util.Set;

import org.apache.log4j.Logger;

import org.intermine.bio.util.OrganismData;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Attribute;
import org.intermine.xml.full.Item;
import org.intermine.xml.full.Reference;

/**
 * Create and store GeneFamily, ConsensusRegion, Gene.geneFamily by querying the chado feature, featureprop and phylotree tables.
 *
 * @author Sam Hokin, NCGR
 */
public class GeneFamilyProcessor extends ChadoProcessor {
	
    private static final Logger LOG = Logger.getLogger(GeneFamilyProcessor.class);

    /**
     * Create a new GeneFamilyProcessor
     *
     * @param chadoDBConverter the ChadoDBConverter that is controlling this processor
     */
    public GeneFamilyProcessor(ChadoDBConverter chadoDBConverter) {
        super(chadoDBConverter);
    }

    /**
     * {@inheritDoc}
     * We process the chado database by reading the phylotree, phylonode, feature and feature_relationship tables
     */
    @Override
    public void process(Connection connection) throws SQLException, ObjectStoreException {
        
        // initialize our DB stuff
        Statement stmt = connection.createStatement();
        ResultSet rs;
        
        Map<Integer,Item> geneMap = new HashMap<>();
        Map<String,Item> geneFamilyMap = new HashMap<>();

        // CV term IDs
        int geneFamilyTypeId = getCVTermId(stmt, "gene family");
        int consensusRegionTypeId = getCVTermId(stmt, "consensus_region");
        int geneTypeId = getCVTermId(stmt, "gene");

        // now grab the gene families from featureprop and put them in a map
        // NOTE: gene families have no chadoId, they are simply names
        rs = stmt.executeQuery("SELECT DISTINCT value FROM featureprop WHERE type_id="+geneFamilyTypeId);
        while (rs.next()) {
            String name = rs.getString("value");
            Item geneFamily = getChadoDBConverter().createItem("GeneFamily");
            geneFamily.setAttribute("primaryIdentifier", name);
            geneFamilyMap.put(name, geneFamily);
        }
        rs.close();
	LOG.info("Created "+geneFamilyMap.size()+" gene families from featureprop query.");

        // now drill through phylotree populating the existing gene family descriptions
        rs = stmt.executeQuery("SELECT DISTINCT name,comment FROM phylotree");
        while (rs.next()) {
            String name = rs.getString("name");
            String description = rs.getString("comment");
            if (geneFamilyMap.containsKey(name)) {
                Item geneFamily = geneFamilyMap.get(name);
                geneFamily.setAttribute("description", description);
            }
        }
        rs.close();
	LOG.info("Gene family descriptions added from phylotree.");

        // Associate consensus regions with gene families
        // HACK: we assume that consensus regions are named by appending "-consensus" to their gene family name!!!!
        rs = stmt.executeQuery("SELECT feature_id,uniquename,name FROM feature WHERE type_id="+consensusRegionTypeId);
        while (rs.next()) {
            int feature_id = rs.getInt("feature_id");
            String uniquename = rs.getString("uniquename");
            String name = rs.getString("name");
            String[] parts = uniquename.split("-"); // [consensus region] = [gene family]-consensus
            String geneFamilyName = parts[0];
            // only store consensus regions associated with gene families we've retrieved
            if (geneFamilyMap.containsKey(geneFamilyName)) {
                Item geneFamily = geneFamilyMap.get(geneFamilyName);
                Item consensusRegion = getChadoDBConverter().createItem("ConsensusRegion");
                consensusRegion.setAttribute("chadoId", String.valueOf(feature_id));
                consensusRegion.setAttribute("primaryIdentifier", uniquename);
                consensusRegion.setAttribute("secondaryIdentifier", name);
                consensusRegion.setAttribute("chadoUniqueName", uniquename);
                consensusRegion.setAttribute("chadoName", name);
                consensusRegion.setReference("geneFamily", geneFamily);
                store(consensusRegion);
            }
        }
        rs.close();
	LOG.info("Consensus regions associated with gene families.");

        // Spin through the gene families and associate their genes with them
        for (String geneFamilyName : geneFamilyMap.keySet()) {
            Item geneFamily = geneFamilyMap.get(geneFamilyName);
            String query = "SELECT feature.feature_id,feature.uniquename,feature.name" +
                " FROM feature,featureprop" +
                " WHERE feature.feature_id=featureprop.feature_id" +
                " AND feature.type_id="+geneTypeId +
                " AND featureprop.type_id="+geneFamilyTypeId +
                " AND featureprop.value='"+geneFamilyName+"'";
            rs = stmt.executeQuery(query);
            while (rs.next()) {
                int chadoId = rs.getInt("feature_id");
                String uniquename = rs.getString("uniquename");
                String name = rs.getString("name");
                Item gene;
                if (geneMap.containsKey(chadoId)) {
                    gene = geneMap.get(chadoId);
                } else {
                    gene = getChadoDBConverter().createItem("Gene");
                    gene.setAttribute("chadoId", String.valueOf(chadoId));
                    gene.setAttribute("primaryIdentifier", uniquename);
                    gene.setAttribute("chadoUniqueName", uniquename);
                    gene.setAttribute("chadoName", name);
                    gene.setReference("geneFamily", geneFamily);
                    geneMap.put(chadoId, gene);
                }
            }
            rs.close();
        }

        // store stuff stored in maps
        store(geneMap.values());
        store(geneFamilyMap.values());
    }

    /**
     * Get the CVTerm ID for a given CVTerm name.
     * @param stmt the database connection statement, initialized to the chado database
     * @param name the desired CV term name
     * @return the CV term id
     * @throws SQLException
     */
     protected int getCVTermId(Statement stmt, String name) throws SQLException {
        ResultSet rs = stmt.executeQuery("SELECT cvterm_id FROM cvterm WHERE name='"+name+"'");
        if (rs.next()) {
            int cvtermId = rs.getInt("cvterm_id");
            rs.close();
            return cvtermId;
        } else {
            throw new RuntimeException("Could not determine CV term id for '"+name+"'.");
        }
     }

    
}
