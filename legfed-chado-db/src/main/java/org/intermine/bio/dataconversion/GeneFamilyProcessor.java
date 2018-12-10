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
 * Create and store GeneFamily, ConsensusRegion, Gene.geneFamily, Homologue and Gene.homologues by querying the chado feature, featureprop and phylotree tables.
 *
 * Homologous genes are defined as genes that share a gene family.
 *
 * project.xml parameters:
 *   organisms "3920"
 *   strains "IT97K-499-35"
 *   homologue.organisms "3920 3702 3847 3880"
 *   homologue.strains "IT97K-499-35 Col0 Williams82  Mt4.0"
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
        Statement stmt1 = connection.createStatement();
        Statement stmt2 = connection.createStatement();
        ResultSet rs1;
        ResultSet rs2;
        
        // same gene multiple times
        Map<String,Item> geneMap = new HashMap<String,Item>();

        // same gene family multiple times
        Map<String,Item> geneFamilyMap = new HashMap<String,Item>();

        // get the desired chado source organism_ids
        Set<Integer> sourceOrganismIds = getChadoDBConverter().getDesiredChadoOrganismIds();
	
	// get the desired chado target (homologue) organism_ids
        Set<Integer> targetOrganismIds = getChadoDBConverter().getDesiredChadoHomologueOrganismIds();

        // CV term IDs
        int geneFamilyTypeId = getCVTermId(stmt1, "gene family");
        int consensusRegionTypeId = getCVTermId(stmt1, "consensus_region");
        int geneTypeId = getCVTermId(stmt1, "gene");

        // now grab the gene families from featureprop and put them in a map
        rs1 = stmt1.executeQuery("SELECT DISTINCT value FROM featureprop WHERE type_id="+geneFamilyTypeId);
        while (rs1.next()) {
            String name = rs1.getString("value");
            Item geneFamily = getChadoDBConverter().createItem("GeneFamily");
            geneFamily.setAttribute("primaryIdentifier", name);
            geneFamilyMap.put(name, geneFamily);
        }
        rs1.close();
	LOG.info("Created "+geneFamilyMap.size()+" gene families from featureprop query.");

        // now drill through phylotree populating the existing gene family descriptions
        rs1 = stmt1.executeQuery("SELECT DISTINCT name,comment FROM phylotree");
        while (rs1.next()) {
            String name = rs1.getString("name");
            String description = rs1.getString("comment");
            if (geneFamilyMap.containsKey(name)) {
                Item geneFamily = geneFamilyMap.get(name);
                geneFamily.setAttribute("description", description);
            }
        }
        rs1.close();
	LOG.info("Gene family descriptions added from phylotree.");

        // Associate consensus regions with gene families
        // HACK: we assume that consensus regions are named by appending "-consensus" to their gene family name!!!!
        rs1 = stmt1.executeQuery("SELECT uniquename,name FROM feature WHERE type_id="+consensusRegionTypeId);
        while (rs1.next()) {
            String uniquename = rs1.getString("uniquename");
            String name = rs1.getString("name");
            String[] parts = uniquename.split("-"); // [consensus region] = [gene family]-consensus
            String geneFamilyName = parts[0];
            // only store consensus regions associated with gene families we've retrieved
            if (geneFamilyMap.containsKey(geneFamilyName)) {
                Item geneFamily = geneFamilyMap.get(geneFamilyName);
                Item consensusRegion = getChadoDBConverter().createItem("ConsensusRegion");
                consensusRegion.setAttribute("primaryIdentifier", uniquename);
                consensusRegion.setAttribute("secondaryIdentifier", name);
                consensusRegion.setAttribute("chadoUniqueName", uniquename);
                consensusRegion.setAttribute("chadoName", name);
                consensusRegion.setReference("geneFamily", geneFamily);
                store(consensusRegion);
            }
        }
	LOG.info("Consensus regions associated with gene families.");

        // Spin through the gene families and find source-target homologues within each one
        for (String geneFamilyName : geneFamilyMap.keySet()) {
            Item geneFamily = geneFamilyMap.get(geneFamilyName);
            // seem to have to store it here rather than later
            store(geneFamily);
            // store the gene names and their organism ID in maps
            Map<String,Integer> sourceGeneMap = new HashMap<String,Integer>();
            Map<String,Integer> targetGeneMap = new HashMap<String,Integer>();
            String query = "SELECT feature.organism_id,feature.uniquename,feature.name" +
                " FROM feature,featureprop " +
                " WHERE feature.feature_id=featureprop.feature_id" +
                " AND feature.type_id="+geneTypeId +
                " AND featureprop.type_id="+geneFamilyTypeId +
                " AND featureprop.value='"+geneFamilyName+"'";
            rs1 = stmt1.executeQuery(query);
            while (rs1.next()) {
                int organism_id = rs1.getInt("organism_id");
                String uniquename = rs1.getString("uniquename");
		String name = rs1.getString("name");
		String key = uniquename+"xxx"+name;
                if (sourceOrganismIds.contains(organism_id)) sourceGeneMap.put(key,organism_id);
                if (targetOrganismIds.contains(organism_id)) targetGeneMap.put(key,organism_id);
            }
            rs1.close();
            // create and store the desired genes and homologs from this gene family
            for (String sourceKey : sourceGeneMap.keySet()) {
                Integer sourceOrganismId = sourceGeneMap.get(sourceKey);
                Item sourceOrganism = getChadoDBConverter().getOrganismItem(sourceOrganismId);
		Item sourceStrain = getChadoDBConverter().getStrainItem(sourceOrganismId);
                Item sourceGene;
                if (geneMap.containsKey(sourceKey)) {
                    sourceGene = geneMap.get(sourceKey);
                } else {
		    String[] parts = sourceKey.split("xxx");
		    String uniquename = parts[0];
		    String name = parts[1];
                    sourceGene = getChadoDBConverter().createItem("Gene");
                    sourceGene.setAttribute("primaryIdentifier", uniquename);
		    sourceGene.setAttribute("secondaryIdentifier", name);
                    sourceGene.setAttribute("chadoUniqueName", uniquename);
		    sourceGene.setAttribute("chadoName", name);
                    sourceGene.setReference("organism", sourceOrganism);
		    if (sourceStrain!=null) sourceGene.setReference("strain", sourceStrain);
                    sourceGene.setReference("geneFamily", geneFamily);
                    geneMap.put(sourceKey, sourceGene);
                }
                for (String targetKey : targetGeneMap.keySet()) {
                    if (!targetKey.equals(sourceKey)) {
                        Integer targetOrganismId = targetGeneMap.get(targetKey);
                        Item targetOrganism = getChadoDBConverter().getHomologueOrganismItem(targetOrganismId);
			Item targetStrain = getChadoDBConverter().getHomologueStrainItem(targetOrganismId);
                        Item targetGene;
                        if (geneMap.containsKey(targetKey)) {
                            targetGene = geneMap.get(targetKey);
                        } else {
			    String[] parts = targetKey.split("xxx");
			    String uniquename = parts[0];
			    String name = parts[1];
                            targetGene = getChadoDBConverter().createItem("Gene");
                            targetGene.setAttribute("primaryIdentifier", uniquename);
			    targetGene.setAttribute("secondaryIdentifier", name);
                            targetGene.setAttribute("chadoUniqueName", uniquename);
			    targetGene.setAttribute("chadoName", name);
                            targetGene.setReference("organism", targetOrganism);
			    if (targetStrain!=null) targetGene.setReference("strain", targetStrain);
                            targetGene.setReference("geneFamily", geneFamily);
                            geneMap.put(targetKey, targetGene);
                        }
                        // homologue
                        Item homologue = getChadoDBConverter().createItem("Homologue");
                        String type = "sameGeneFamily"; // not really known which sort of homologue they are
                        homologue.setAttribute("type", type);
                        homologue.setReference("geneFamily", geneFamily);
                        homologue.setReference("gene", sourceGene);
                        homologue.setReference("homologue", targetGene);
                        store(homologue); // store now, not going to touch it later
                    }
                }
            }
        }

        // store the genes since it doesn't seem to work if they're stored right after creation
        LOG.info("Storing "+geneMap.size()+" Gene items.");
        for (Item gene : geneMap.values()) {
            store(gene);
        }

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
