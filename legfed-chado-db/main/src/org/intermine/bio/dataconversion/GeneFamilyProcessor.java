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
 * Create and store GeneFamily, Gene.geneFamily, Homologue and Gene.homologues by querying the chado feature, featureprop and phylotree tables.
 *
 * Homologous genes are defined as genes that share a gene family.
 *
 * Since this processor deals only with chado data, Items are stored in maps with Integer keys equal to
 * the chado feature.feature_id.
 *
 * project.xml parameters:
 *   organisms stored as Homologue.gene e.g. "3920_IT97K-499-35"
 *   homologue.organisms stored as Homologue.homologue e.g. "3920_IT97K-499-35 3702_Col0 3847_Williams82 3880_Mt4.0 3885_G19833 130453_V14167 130454_K30076 3827_CDCFrontier 3827_ICC4958"
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
        
	// store the full set of organisms in one place for storage since source and target organisms can overlap
	Map<Integer,Item> organismMap = new HashMap<Integer,Item>();

        // same gene multiple times
        Map<String,Item> geneMap = new HashMap<String,Item>();

        // same gene family multiple times
        Map<String,Item> geneFamilyMap = new HashMap<String,Item>();
        
        // build the Organism map from the supplied source taxon IDs
        Set<Integer> sourceOrganisms = new HashSet<Integer>();
        Map<Integer,OrganismData> chadoToOrgData = getChadoDBConverter().getChadoIdToOrgDataMap();
        for (Map.Entry<Integer,OrganismData> entry : chadoToOrgData.entrySet()) {
            Integer organismId = entry.getKey();
            OrganismData organismData = entry.getValue();
            int taxonId = organismData.getTaxonId();
            String variety = organismData.getVariety();
            Item organism = getChadoDBConverter().createItem("Organism");
            organism.setAttribute("taxonId", String.valueOf(taxonId));
            organism.setAttribute("variety", variety); // required
            store(organism);
	    organismMap.put(organismId, organism);
            sourceOrganisms.add(organismId);
        }
        LOG.info("Created "+sourceOrganisms.size()+" source organism Items.");
        if (sourceOrganisms.size()==0) {
            throw new RuntimeException("Property organisms must contain at least one taxon ID in project.xml.");
        }

        // add to the Organism map from the requested target taxon IDs
        Set<Integer> targetOrganisms = new HashSet<Integer>();
        Map<Integer,OrganismData> chadoToHomologueOrgData = getChadoDBConverter().getChadoIdToHomologueOrgDataMap();
        for (Map.Entry<Integer,OrganismData> entry : chadoToHomologueOrgData.entrySet()) {
            Integer organismId = entry.getKey();
            if (!organismMap.containsKey(organismId)) {
                OrganismData organismData = entry.getValue();
                int taxonId = organismData.getTaxonId();
                String variety = organismData.getVariety();
                Item organism = getChadoDBConverter().createItem("Organism");
                organism = getChadoDBConverter().createItem("Organism");
                organism.setAttribute("taxonId", String.valueOf(taxonId));
                organism.setAttribute("variety", variety); // required
                store(organism);
                organismMap.put(organismId, organism);
            }
            targetOrganisms.add(organismId);
        }
        LOG.info("Created "+targetOrganisms.size()+" target organism Items.");
        if (targetOrganisms.size()==0) {
            throw new RuntimeException("Property homologue.organisms must contain at least one taxon ID in project.xml.");
        }

        // CV term IDs
        int geneFamilyTypeId = 0;
        rs1 = stmt1.executeQuery("SELECT cvterm_id FROM cvterm WHERE name='gene family'");
        if (rs1.next()) geneFamilyTypeId = rs1.getInt("cvterm_id");
        rs1.close();
        if (geneFamilyTypeId==0) throw new RuntimeException("Could not determine CV term id for 'gene family'.");

        int geneTypeId = 0;
        rs1 = stmt1.executeQuery("SELECT cvterm_id FROM cvterm WHERE name='gene'");
        if (rs1.next()) geneTypeId = rs1.getInt("cvterm_id");
        rs1.close();
        if (geneTypeId==0) throw new RuntimeException("Could not determine CV term id for 'gene'.");

        int consensusRegionId = 0;
        rs1 = stmt1.executeQuery("SELECT cvterm_id FROM cvterm WHERE name='consensus_region'");
        if (rs1.next()) consensusRegionId = rs1.getInt("cvterm_id");
        rs1.close();
        if (consensusRegionId==0) throw new RuntimeException("Could not determine CV term id for 'consensus_region'.");

        // create a comma-separated list of target organism IDs for WHERE organism_id IN query clause
        String targetOrganismSQL = "(";
        boolean first = true;
        for (int organism_id : targetOrganisms) {
            if (first) {
                targetOrganismSQL += organism_id;
                first = false;
            } else {
                targetOrganismSQL += ","+organism_id;
            }
        }
        targetOrganismSQL += ")";
        LOG.info("Querying target gene records with organism_id in "+targetOrganismSQL);

        // now grab the gene families from featureprop
        rs1 = stmt1.executeQuery("SELECT DISTINCT value FROM featureprop WHERE type_id="+geneFamilyTypeId);
        while (rs1.next()) {
            String name = rs1.getString("value");
            Item geneFamily = getChadoDBConverter().createItem("GeneFamily");
            geneFamily.setAttribute("primaryIdentifier", name);
            geneFamilyMap.put(name, geneFamily);
        }
        rs1.close();
	LOG.info(geneFamilyMap.size()+" gene families retrieved from featureprop.");

        // now drill through phylotree populating the existing gene family descriptions
        rs1 = stmt1.executeQuery("SELECT name,comment FROM phylotree");
        while (rs1.next()) {
            String name = rs1.getString("name");
            String description = rs1.getString("comment");
            Item geneFamily = geneFamilyMap.get(name);
            if (geneFamily!=null) geneFamily.setAttribute("description", description);
        }
        rs1.close();
	LOG.info("Gene family descriptions added from phylotree.");

        // Get the consensus regions and their sequences and associate them with gene families
        rs1 = stmt1.executeQuery("SELECT * FROM feature WHERE type_id="+consensusRegionId);
        while (rs1.next()) {
            String uniquename = rs1.getString("uniquename");
            String residues = rs1.getString("residues");
            int seqlen = rs1.getInt("seqlen");
            String[] parts = uniquename.split("-"); // [consensus region] = [gene family]-consensus
            String geneFamilyName = parts[0];
            // only store consensus regions that match a gene family we've retrieved
            if (geneFamilyMap.containsKey(geneFamilyName)) {
                Item geneFamily = geneFamilyMap.get(geneFamilyName);
                Item consensusRegion = getChadoDBConverter().createItem("ConsensusRegion");
                consensusRegion.setAttribute("primaryIdentifier", uniquename);
                consensusRegion.setAttribute("length", String.valueOf(seqlen));
                consensusRegion.setReference("geneFamily", geneFamily);
                if (residues!=null) {
                    Item sequence = getChadoDBConverter().createItem("Sequence");
                    sequence.setAttribute("residues", residues);
                    sequence.setAttribute("length", String.valueOf(seqlen));
                    store(sequence); // store now, not going to touch it later
                    consensusRegion.setReference("sequence", sequence);
                }
                store(consensusRegion); // store now, not going to touch it later
                geneFamily.setReference("consensusRegion", consensusRegion); // reverse-reference
                store(geneFamily);      // we're done with this gene family
            }
        }
	LOG.info("Consensus region records created for gene families.");

        // Now query the source organism genes and set Gene.geneFamily; then query the target organisms for genes in the same gene family
        for (int source_organism_id : sourceOrganisms) {
            Item sourceOrganism = organismMap.get(source_organism_id);
            // query the source genes gene family
            rs1 = stmt1.executeQuery("SELECT feature.uniquename,featureprop.value" +
                                     " FROM feature,featureprop"+
                                     " WHERE feature.feature_id=featureprop.feature_id"+
                                     " AND feature.organism_id="+source_organism_id+
                                     " AND feature.type_id="+geneTypeId+
                                     " AND featureprop.type_id="+geneFamilyTypeId);
            while (rs1.next()) {
                String sourceGeneName = rs1.getString("uniquename");
                String geneFamilyName = rs1.getString("value");
                Item geneFamily = geneFamilyMap.get(geneFamilyName);
                if (geneFamily!=null) {
                    // get the source gene; could already be here from a paralogue
                    Item sourceGene;
                    if (geneMap.containsKey(sourceGeneName)) {
                        sourceGene = geneMap.get(sourceGeneName);
                    } else {
                        sourceGene = getChadoDBConverter().createItem("Gene");
                        sourceGene.setAttribute("primaryIdentifier", sourceGeneName);
                        sourceGene.setReference("organism", sourceOrganism);
                        sourceGene.setReference("geneFamily", geneFamily);
                        store(sourceGene);
                        geneMap.put(sourceGeneName, sourceGene);
                    }
                    // query the target genes for homology
                    rs2 = stmt2.executeQuery("SELECT feature.organism_id,feature.uniquename FROM feature,featureprop"+
                                             " WHERE feature.feature_id=featureprop.feature_id"+
                                             " AND feature.organism_id IN "+targetOrganismSQL+
                                             " AND feature.type_id="+geneTypeId+
                                             " AND featureprop.type_id="+geneFamilyTypeId+
                                             " AND featureprop.value='"+geneFamilyName+"'"+
                                             " AND uniquename!='"+sourceGeneName+"'");
                    while (rs2.next()) {
                        // target gene
                        int target_organism_id = rs2.getInt("organism_id");
                        String targetGeneName = rs2.getString("uniquename");
                        Item targetGene;
                        if (geneMap.containsKey(targetGeneName)) {
                            targetGene = geneMap.get(targetGeneName);
                        } else {
                            Item targetOrganism = organismMap.get(target_organism_id);
                            targetGene = getChadoDBConverter().createItem("Gene");
                            targetGene.setAttribute("primaryIdentifier", targetGeneName);
                            targetGene.setReference("organism", targetOrganism);
                            targetGene.setReference("geneFamily", geneFamily);
                            store(targetGene);
                            geneMap.put(targetGeneName, targetGene);
                        }
                        // homologue
                        Item homologue = getChadoDBConverter().createItem("Homologue");
                        String type = "orthologue";
                        if (source_organism_id==target_organism_id) type = "paralogue";
                        homologue.setAttribute("type", type);
                        homologue.setReference("geneFamily", geneFamily);
                        homologue.setReference("gene", sourceGene);
                        homologue.setReference("homologue", targetGene);
                        store(homologue); // store now, not going to touch it later
                        // add homologue to source gene's collection
                        sourceGene.addToCollection("homologues", homologue);
                    }
                    rs2.close();
                }
            }
            rs1.close();
        }

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
