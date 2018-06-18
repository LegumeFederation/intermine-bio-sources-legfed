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

import java.io.BufferedReader;
import java.io.FileReader;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.Map;
import java.util.HashMap;
import java.util.Set;
import java.util.HashSet;

import org.apache.log4j.Logger;
import org.apache.log4j.LogManager;

import org.intermine.bio.util.OrganismData;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Attribute;
import org.intermine.xml.full.Item;
import org.intermine.xml.full.Reference;

/**
 * Read plant reactome data in from a tab-delimited file, load pathways and associate them with genes.
 * We do NOT store pathways that are not associated with any of our genes.
 * The file fields are:
 * <pre>
 * identifier      name                    species           gene
 * R-CCN-1119506.1 tyrosine degradation I  Coffea canephora  Cc05_g11510
 * </pre>
 * The file is given by the parameter reactome.file in project.xml. It is likely gene_ids_by_pathway_and_species.tab from Gramene.
 *
 * NOTE: special species-specific treatment is required, as the Plant Reactome gene names do NOT match the LIS gene names by a long stretch.
 *
 * @author Sam Hokin, NCGR
 */
public class ReactomeProcessor extends ChadoProcessor {
	
    private static final Logger LOG = LogManager.getLogger(ReactomeProcessor.class);

    /**
     * Create a new ReactomeProcessor
     * @param chadoDBConverter the ChadoDBConverter that is controlling this processor
     */
    public ReactomeProcessor(ChadoDBConverter chadoDBConverter) {
        super(chadoDBConverter);
    }

    /**
     * {@inheritDoc}
     * Process by spinning through the reactome file and associating pathways for our organisms with genes from chado.
     */
    @Override
    public void process(Connection connection) throws Exception {
        
        // Items stored in maps to avoid dupe stores
	Map<Integer,Item> organismMap = new HashMap<Integer,Item>();  // restrict to desired organisms
        Map<Integer,Item> geneMap = new HashMap<Integer,Item>();  // keyed by feature_id
        Map<String,Item> pathwayMap = new HashMap<String,Item>();    // keyed by identifier
        
        // store desired species strings in a Set, hopefully spelled same way!
        Set<String> speciesSet = new HashSet<String>();              // "Phaseolus vulgaris"
        
        // initialize our DB statement
        Statement stmt = connection.createStatement();
        ResultSet rs;
        
        // build the organism map from the supplied taxon_variety IDs
        Map<Integer,OrganismData> chadoToOrgData = getChadoDBConverter().getChadoIdToOrgDataMap();
        for (Integer organismId : chadoToOrgData.keySet()) {
            OrganismData organismData = chadoToOrgData.get(organismId);
            int taxonId = organismData.getTaxonId();
            String genus = organismData.getGenus();
            String species = organismData.getSpecies();
            String variety = organismData.getVariety();
            speciesSet.add(genus+" "+species);
            Item organism = getChadoDBConverter().createItem("Organism");
            organism.setAttribute("taxonId", String.valueOf(taxonId));
            organism.setAttribute("genus", genus);
            organism.setAttribute("species", species);
            organism.setAttribute("variety", variety);
            store(organism);
	    organismMap.put(organismId, organism);
        }
        LOG.info("Species:"+speciesSet);

        // get the cv term ID for gene
        rs = stmt.executeQuery("SELECT cvterm_id FROM cvterm WHERE name='gene'");
        rs.next();
        int geneCVTermId = rs.getInt("cvterm_id");
        rs.close();

        // load the reactome file into a map
        String reactomeFilename = getChadoDBConverter().getReactomeFilename();
        LOG.info("Reading "+reactomeFilename);
        BufferedReader reader = new BufferedReader(new FileReader(reactomeFilename));
        String line = null;
        while ((line=reader.readLine())!=null) {
            if (!line.startsWith("#")) {
                String[] fields = line.split("\t");
                String identifier = fields[0];
                String pathwayName = fields[1];
                String species = fields[2];
                String geneName = fields[3];
                if (speciesSet.contains(species)) {
                    String featureName = null;
                    // HACK: concoct the correct chado feature.name
                    if (geneName.contains("_")) {
                        String[] parts = geneName.split("_");
                        if (parts[0].equals("GLYMA")) {
                            featureName = "glyma.Glyma."+parts[1];
                        }
                    }
                    if (featureName!=null) {
                        String query = "SELECT * FROM feature WHERE type_id="+geneCVTermId+" AND name='"+featureName+"'";
                        rs = stmt.executeQuery(query);
                        if (rs.next()) {
                            // get or store the gene
                            Item gene = null;
                            int feature_id = rs.getInt("feature_id");
                            if (geneMap.containsKey(feature_id)) {
                                gene = geneMap.get(feature_id);
                            } else {
                                int organism_id = rs.getInt("organism_id");
                                Item organism = organismMap.get(organism_id);
                                gene = getChadoDBConverter().createItem("Gene");
                                gene.setAttribute("primaryIdentifier", rs.getString("uniquename"));
                                gene.setAttribute("secondaryIdentifier", rs.getString("name"));
                                gene.setReference("organism", organism);
                                geneMap.put(feature_id, gene);
                            }
                            // get or store the pathway, and store species, not organism (multiple varieties have same pathways)
                            Item pathway = null;
                            if (pathwayMap.containsKey(identifier)) {
                                pathway = pathwayMap.get(identifier);
                            } else {
                                pathway = getChadoDBConverter().createItem("Pathway");
                                pathway.setAttribute("identifier", identifier);
                                pathway.setAttribute("name", pathwayName);
                                pathway.setAttribute("species", species);
                                pathwayMap.put(identifier, pathway);
                            }
                            // associate gene with pathway
                            gene.addToCollection("pathways", pathway);
                        } else {
                            LOG.error(featureName+" not found in chado database.");
                        }
                        rs.close();
                    }
                }
            }
        }

        // store the maps here so that the collection get stored
        for (Item gene : geneMap.values()) store(gene);
        for (Item pathway : pathwayMap.values()) store(pathway);
        
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

}
