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
 * Create and store GeneFamily, Gene.geneFamily, Homologue and Gene.homologues by querying the chado phylotree and phylonode tables, as well as feature.
 * This means that homologous genes are defined as genes that share a gene family.
 *
 * phylotree.name = gene family
 * phylotree.phylotree_id -- phylonode.phylotree_id
 * phylonode.feature_id = polypeptide corresponding to a homolog in this family
 * phylonode.feature_id = feature_relationship.subject_id gives object_id = mrna_id of the corresponding mRNA
 * mrna_id = feature_relationship.subject_id gives object_id = gene_id of the corresponding gene
 *
 * Since this processor deals only with chado data, Items are stored in maps with Integer keys equal to
 * the chado feature.feature_id.
 *
 * project.xml parameters:
 *
 * organisms = taxon IDs of desired organisms to process as source genes for Homologue.gene
 * homologue.organisms = taxon IDs of desired organisms to process as homologue genes for Homologue.homologue
 *
 * @author Sam Hokin, NCGR
 */
public class HomologyProcessor extends ChadoProcessor {
	
    private static final Logger LOG = Logger.getLogger(HomologyProcessor.class);

    /**
     * Create a new HomologyProcessor
     * @param chadoDBConverter the ChadoDBConverter that is controlling this processor
     */
    public HomologyProcessor(ChadoDBConverter chadoDBConverter) {
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
        Statement stmt3 = connection.createStatement();
        ResultSet rs1;
        ResultSet rs2;
        ResultSet rs3;
        
        // -----------------------------------------------------
        // ---------------- PARAMETER CHECKS -------------------
        // -----------------------------------------------------

        // Phytozome version must be specified in phytozome.version
        String phytozomeVersion = getChadoDBConverter().getPhytozomeVersion();
        if (phytozomeVersion==null || phytozomeVersion.length()==0) {
            throw new RuntimeException("Property phytozome.version must be specified in project.xml.");
        }

        // ---------------------------------------------------------
        // ---------------- INITIAL DATA LOADING -------------------
        // ---------------------------------------------------------

	// store the full set of organisms in one place for storage since source and target organisms can overlap
	Map<Integer,Item> organismMap = new HashMap<Integer,Item>();
        
        // build the Organism map from the supplied taxon IDs
        Map<Integer,Item> sourceOrganismMap = new HashMap<Integer,Item>();
        Map<Integer,OrganismData> chadoToOrgData = getChadoDBConverter().getChadoIdToOrgDataMap();
        for (Map.Entry<Integer,OrganismData> entry : chadoToOrgData.entrySet()) {
            Integer organismId = entry.getKey();
            OrganismData organismData = entry.getValue();
            int taxonId = organismData.getTaxonId();
            String variety = organismData.getVariety();
            Item organism = getChadoDBConverter().createItem("Organism");
            BioStoreHook.setSOTerm(getChadoDBConverter(), organism, "organism", getChadoDBConverter().getSequenceOntologyRefId());
            organism.setAttribute("taxonId", String.valueOf(taxonId));
            if (variety!=null) organism.setAttribute("variety", variety);
            sourceOrganismMap.put(organismId, organism);
	    organismMap.put(organismId, organism);
        }
        LOG.info("Created "+sourceOrganismMap.size()+" source organism Items.");
        if (sourceOrganismMap.size()==0) {
            throw new RuntimeException("Property organisms must contain at least one taxon ID in project.xml.");
        }

        // build the homologue Organism map from the requested taxon IDs
        Map<Integer,Item> targetOrganismMap = new HashMap<Integer,Item>();
        Map<Integer,OrganismData> chadoToHomologueOrgData = getChadoDBConverter().getChadoIdToHomologueOrgDataMap();
        for (Map.Entry<Integer,OrganismData> entry : chadoToHomologueOrgData.entrySet()) {
            Integer organismId = entry.getKey();
            OrganismData organismData = entry.getValue();
            int taxonId = organismData.getTaxonId();
            String variety = organismData.getVariety();
            Item organism = getChadoDBConverter().createItem("Organism");
            BioStoreHook.setSOTerm(getChadoDBConverter(), organism, "organism", getChadoDBConverter().getSequenceOntologyRefId());
            organism.setAttribute("taxonId", String.valueOf(taxonId));
            if (variety!=null) organism.setAttribute("variety", variety);
            targetOrganismMap.put(organismId, organism);
	    if (!organismMap.containsKey(organismId)) organismMap.put(organismId, organism);
        }
        LOG.info("Created "+targetOrganismMap.size()+" target organism Items.");
        if (targetOrganismMap.size()==0) {
            throw new RuntimeException("Property homologue.organisms must contain at least one taxon ID in project.xml.");
        }

	// store the organisms
	for (Item organism : organismMap.values()) {
            store(organism);
	}

        // create a comma-separated list of organism IDs for WHERE organism_id IN query clauses
        String organismSQL = "(";
        boolean first = true;
        for (Integer organismId : sourceOrganismMap.keySet()) {
            if (first) {
                organismSQL += organismId;
                first = false;
            } else {
                organismSQL += ","+organismId;
            }
        }
        for (Integer organismId : targetOrganismMap.keySet()) {
            organismSQL += ","+organismId;
        }
        organismSQL += ")";
        LOG.info("Querying phylonode records with organism_id in "+organismSQL);

        // run through the phylotree, creating and storing Homologue items for each gene family
        // geneMap contains ALL genes for final storage
        // also store gene families and consensus regions for relation to gene family
        HashMap<String,Item> geneMap = new HashMap<String,Item>();
        rs1 = stmt1.executeQuery("SELECT * FROM phylotree WHERE name LIKE '"+phytozomeVersion+".%'");
        while (rs1.next()) {
            int phylotree_id = rs1.getInt("phylotree_id");
            String name = rs1.getString("name");
            String description = rs1.getString("comment");
            Item geneFamily = getChadoDBConverter().createItem("GeneFamily");
            geneFamily.setAttribute("primaryIdentifier", name);
            geneFamily.setAttribute("description", description);
            // ASSUME that consensus region primaryIdentifier = geneFamily.primaryIdentifier+"-consensus" to relate to consensus region
            Item consensusRegion = getChadoDBConverter().createItem("ConsensusRegion");
            consensusRegion.setAttribute("primaryIdentifier", name+"-consensus");
            consensusRegion.setReference("geneFamily", geneFamily);
            geneFamily.setReference("consensusRegion", consensusRegion);
            // store 'em
            store(geneFamily);
            store(consensusRegion);
            // query the members of this gene family, polypeptides, that are in the desired organisms by referencing the feature table
            Set<Item> sourceSet = new HashSet<Item>();
            Set<Item> targetSet = new HashSet<Item>();
            rs2 = stmt2.executeQuery("SELECT feature.* FROM phylonode,feature WHERE phylonode.feature_id=feature.feature_id AND organism_id IN "+organismSQL+" AND phylotree_id="+phylotree_id);
            while (rs2.next()) {
                // organism
                Integer organismId = new Integer(rs2.getInt("organism_id"));
                Item organism;
                if (sourceOrganismMap.containsKey(organismId)) {
                    organism = sourceOrganismMap.get(organismId);
                } else {
                    organism = targetOrganismMap.get(organismId);
                }
                // extract the gene via two feature_relationship relations: polypeptide -> mRNA -> gene
                int polypeptide_id = rs2.getInt("feature_id");
                rs3 = stmt3.executeQuery("SELECT * FROM feature WHERE feature_id=(" +
                                         "SELECT object_id FROM feature_relationship WHERE subject_id=(" +
                                         "SELECT object_id FROM feature_relationship WHERE subject_id="+polypeptide_id+"))");
                if (rs3.next()) {
                    String geneName = rs3.getString("uniquename");
                    if (geneMap.containsKey(geneName)) {
                        // get the already stored gene
                        Item gene = geneMap.get(geneName);
                        if (sourceOrganismMap.containsKey(organismId)) sourceSet.add(gene);
                        if (targetOrganismMap.containsKey(organismId)) targetSet.add(gene);
                    } else {
                        // store this new gene
                        Item gene = getChadoDBConverter().createItem("Gene");
                        Item sequence = getChadoDBConverter().createItem("Sequence");
                        ChadoFeature cf = new ChadoFeature(rs3);
                        cf.populateSequenceFeature(gene, sequence, organism);
                        // associate this gene family with this gene (reverse-referenced to GeneFamily)
                        gene.setReference("geneFamily", geneFamily);
                        // store 'em
                        store(gene);
                        store(sequence);
                        // put in maps and sets
                        geneMap.put(geneName, gene);
                        if (sourceOrganismMap.containsKey(organismId)) sourceSet.add(gene);
                        if (targetOrganismMap.containsKey(organismId)) targetSet.add(gene);
                    }
                }
                // // no such gene in the feature table, we'll make one up from our derived name with no other info
                //         Item gene = getChadoDBConverter().createItem("Gene");
                //         gene.setAttribute("primaryIdentifier", geneName);
                //         gene.setReference("organism", organism);
                //         // associate this gene family with this gene (reverse-referenced to GeneFamily)
                //         gene.setReference("geneFamily", geneFamily);
                //         store(gene);
                //         geneMap.put(geneName, gene);
                //         if (sourceOrganismMap.containsKey(organismId)) sourceSet.add(gene);
                //         if (targetOrganismMap.containsKey(organismId)) targetSet.add(gene);
                rs3.close();
            }
            rs2.close();
            // create a Homologue Item for every combination of sourceSet and targetSet items
            for (Item sourceGene : sourceSet) {
                String sourceGeneId = sourceGene.getAttribute("primaryIdentifier").getValue();
                String sourceOrganismRefId = sourceGene.getReference("organism").getRefId();
                for (Item targetGene : targetSet) {
                    String targetGeneId = targetGene.getAttribute("primaryIdentifier").getValue();
                    String targetOrganismRefId = targetGene.getReference("organism").getRefId();
                    if (!sourceGeneId.equals(targetGeneId)) {
                        String type = "orthologue";
                        if (sourceOrganismRefId.equals(targetOrganismRefId)) type = "paralogue";
                        // forward reference source to target
                        Item forward = getChadoDBConverter().createItem("Homologue");
                        forward.setAttribute("type", type);
                        forward.setReference("geneFamily", geneFamily);
                        forward.setReference("gene", sourceGene);
                        forward.setReference("homologue", targetGene);
                        store(forward);
                        // reverse reference target to source
                        Item reverse = getChadoDBConverter().createItem("Homologue");
                        reverse.setAttribute("type", type);
                        reverse.setReference("geneFamily", geneFamily);
                        reverse.setReference("gene", targetGene);
                        reverse.setReference("homologue", sourceGene);
                        store(reverse);
                    }
                }
            }
        }
        rs1.close();

    }

    /**
     * Store the item.
     * @param item the Item
     * @return the database id of the new Item
     * @throws ObjectStoreException if an error occurs while storing
     */
    protected Integer store(Item item) throws ObjectStoreException {
        return getChadoDBConverter().store(item);
    }
    
    /**
     * Do any extra processing that is needed before the converter starts querying features
     * @param connection the Connection
     * @throws ObjectStoreException if there is a object store problem
     * @throws SQLException if there is a database problem
     */
    protected void earlyExtraProcessing(Connection connection) throws ObjectStoreException, SQLException {
        // override in subclasses as necessary
    }

    /**
     * Do any extra processing for this database, after all other processing is done
     * @param connection the Connection
     * @param featureDataMap a map from chado feature_id to data for that feature
     * @throws ObjectStoreException if there is a problem while storing
     * @throws SQLException if there is a problem
     */
    protected void extraProcessing(Connection connection, Map<Integer, FeatureData> featureDataMap)
        throws ObjectStoreException, SQLException {
        // override in subclasses as necessary
    }

    /**
     * Perform any actions needed after all processing is finished.
     * @param connection the Connection
     * @param featureDataMap a map from chado feature_id to data for that feature
     * @throws SQLException if there is a problem
     */
    protected void finishedProcessing(Connection connection, Map<Integer, FeatureData> featureDataMap) throws SQLException {
        // override in subclasses as necessary
    }

}
