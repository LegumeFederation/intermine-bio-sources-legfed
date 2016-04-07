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

import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeSet;
import java.util.Set;

import org.apache.log4j.Logger;
import org.apache.commons.lang3.StringUtils;

import org.intermine.bio.util.OrganismData;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Attribute;
import org.intermine.xml.full.Item;
import org.intermine.xml.full.Reference;

/**
 * Store the GO annotations contained in the feature.note field for genes and gene families in the LIS chado database.
 *
 * Since this processer deals only with chado data, Items are stored in maps with Integer keys equal to
 * the chado feature.feature_id.
 *
 * @author Sam Hokin, NCGR
 */
public class GOProcessor extends ChadoProcessor {
	
    private static final Logger LOG = Logger.getLogger(GeneticProcessor.class);

    /**
     * Create a new GeneticProcessor
     * @param chadoDBConverter the ChadoDBConverter that is controlling this processor
     */
    public GOProcessor(ChadoDBConverter chadoDBConverter) {
        super(chadoDBConverter);
    }

    /**
     * {@inheritDoc}
     * We process the chado database by reading the feature records for genes.
     */
    @Override
    public void process(Connection connection) throws SQLException, ObjectStoreException {

        LOG.info("Starting GOProcessor.process()");
        
        // initialize our DB statement
        Statement stmt = connection.createStatement();
        
        // build an organism map from the supplied taxon IDs
        Map<Integer,Item> organismMap = new HashMap<Integer,Item>();
        Map<Integer,OrganismData> chadoToOrgData = getChadoDBConverter().getChadoIdToOrgDataMap();
        for (Map.Entry<Integer,OrganismData> entry : chadoToOrgData.entrySet()) {
            Integer organismId = entry.getKey();
            OrganismData organismData = entry.getValue();
            int taxonId = organismData.getTaxonId();
            Item organism = getChadoDBConverter().createItem("Organism");
            organism.setAttribute("taxonId", String.valueOf(taxonId));
            store(organism);
            organismMap.put(organismId, organism);
        }
        LOG.info("Created and stored "+organismMap.size()+" organism Items.");

        // we'll store the GOTerm items in a map to avoid duplication, keyed by identifier (e.g. "GO:000037")
        Map<String,Item> goTermMap = new HashMap<String,Item>();

        // loop over the organisms to fill the GO terms
        for (Map.Entry<Integer,Item> orgEntry : organismMap.entrySet()) {
            int organism_id = orgEntry.getKey().intValue();
            Item organism = orgEntry.getValue();
        
            // load the relevant genes from the gene table, then parse out the GO identifiers
            String query = "SELECT * FROM gene WHERE organism_id="+organism_id;
            LOG.info("executing query: "+query);
            ResultSet rs = stmt.executeQuery(query);
            while (rs.next()) {
                // parse the description for GO identifiers, creating a GOAnnotation each time, and adding it to the gene's collection
                String description = rs.getString("description");
                String[] goNumbers = StringUtils.substringsBetween(description, "GO:", " ");
                if (goNumbers!=null) {
                    // create the Gene item and store the minimal stuff required for merging (and note that gene.symbol is bogus)
                    Item gene = getChadoDBConverter().createItem("Gene");
                    gene.setReference("organism", organism);
                    gene.setAttribute("primaryIdentifier", rs.getString("uniquename"));
                    // add the GO terms
                    for (int j=0; j<goNumbers.length; j++) {
                        String identifier = "GO:"+goNumbers[j];
                        // get the GO term from the map if it's there; otherwise create, store and add it to the map.
                        Item goTerm;
                        if (goTermMap.containsKey(identifier)) {
                            goTerm = goTermMap.get(identifier);
                        } else {
                            goTerm = getChadoDBConverter().createItem("GOTerm");
                            goTerm.setAttribute("identifier", identifier);
                            store(goTerm);
                            goTermMap.put(identifier, goTerm);
                        }
                        // create and store the GOAnnotation linking this gene to this GO term
                        Item goAnnotation = getChadoDBConverter().createItem("GOAnnotation");
                        goAnnotation.setReference("subject", gene);
                        goAnnotation.setReference("ontologyTerm", goTerm);
                        store(goAnnotation);
                        // have to manually set reverse reference since no reverse-reference from subject defined in OntologyAnnotation
                        gene.addToCollection("goAnnotation", goAnnotation);
                    }
                    // store the gene
                    store(gene);
                }
            }
            rs.close();
            
        } // organism

        // load the gene families from the phylotree table, since it contains the comments which may contain the GO terms
        String query = "SELECT name,comment FROM phylotree";
        LOG.info("executing query: "+query);
        ResultSet rs = stmt.executeQuery(query);
        while (rs.next()) {
            String name = rs.getString("name");
            String description = rs.getString("comment");
            String[] goNumbers = StringUtils.substringsBetween(description, "GO:", " ");
            if (goNumbers!=null) {
                // create the GeneFamily item and load the minimal attributes and GO terms
                Item geneFamily = getChadoDBConverter().createItem("GeneFamily");
                geneFamily.setAttribute("primaryIdentifier", name);
                // add the GO terms
                for (int j=0; j<goNumbers.length; j++) {
                    String identifier = "GO:"+goNumbers[j];
                    // get the GO term from the map if it's there; otherwise create, store and add it to the map.
                    Item goTerm;
                    if (goTermMap.containsKey(identifier)) {
                        goTerm = goTermMap.get(identifier);
                    } else {
                        goTerm = getChadoDBConverter().createItem("GOTerm");
                        goTerm.setAttribute("identifier", identifier);
                        store(goTerm);
                        goTermMap.put(identifier, goTerm);
                    }
                    // create and store the GOAnnotation linking this gene family to this GO term
                    Item goAnnotation = getChadoDBConverter().createItem("GOAnnotation");
                    goAnnotation.setReference("subject", geneFamily);
                    goAnnotation.setReference("ontologyTerm", goTerm);
                    store(goAnnotation);
                    // have to manually set reverse reference since no reverse-reference defined in go-annotation
                    geneFamily.addToCollection("goAnnotation", goAnnotation);
                }
                // store the gene family
                store(geneFamily);
            }
        }
        rs.close();

    } // process

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
