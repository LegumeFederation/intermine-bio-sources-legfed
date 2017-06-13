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
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import org.intermine.bio.util.OrganismData;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Attribute;
import org.intermine.xml.full.Item;
import org.intermine.xml.full.Reference;

/**
 * Read plant reactome data in from a tab-delimited file, load pathways and associate them with proteins. The file fields are:
 * <pre>
 * pathwayId pathwayName taxonId proteinName
 * </pre>
 * The file is given by the parameter reactome.file in project.xml.
 *
 * The proteinName is meant to match the feature.name of polypeptides in chado, without the .1, .2 suffix. Every polypeptide with base name matching
 * proteinName will be associated with the pathway. See the HACK below where chado names are matched with reactome names on a per-organism basis.
 *
 * @author Sam Hokin, NCGR
 */
public class ReactomeProcessor extends ChadoProcessor {
	
    private static final Logger LOG = Logger.getLogger(ReactomeProcessor.class);

    /**
     * Create a new ReactomeProcessor
     * @param chadoDBConverter the ChadoDBConverter that is controlling this processor
     */
    public ReactomeProcessor(ChadoDBConverter chadoDBConverter) {
        super(chadoDBConverter);
    }

    /**
     * {@inheritDoc}
     * We process the chado database by reading the feature, featureloc, featurepos, featuremap, feature_relationship and featureprop tables.
     */
    @Override
    public void process(Connection connection) throws Exception {
        
        // all Items are stored in maps and saved at the end
        Map<Integer,Item> organismMap = new HashMap<Integer,Item>(); // keyed by chado organism_id
        Map<String,Item> pathwayMap = new HashMap<String,Item>();    // keyed by name, which is all there is
        Map<Integer,Item> proteinMap = new HashMap<Integer,Item>();  // keyed by feature_id since stored from a chado query

        // initialize our DB statement
        Statement stmt = connection.createStatement();
        ResultSet rs;
        
        // get the cv terms for our items of interest
        rs = stmt.executeQuery("SELECT cvterm_id FROM cvterm WHERE name='polypeptide'");
        rs.next();
        int polypeptideCVTermId = rs.getInt("cvterm_id");
        rs.close();

        // build the Organism map from the supplied taxon IDs and varieties
        // NOTE: these should match the records you want in the reactome file
        Map<Integer,OrganismData> chadoToOrgData = getChadoDBConverter().getChadoIdToOrgDataMap();
        for (Integer organismId : chadoToOrgData.keySet()) {
            OrganismData organismData = chadoToOrgData.get(organismId);
            int taxonId = organismData.getTaxonId();
            String variety = organismData.getVariety();
            Item organism = getChadoDBConverter().createItem("Organism");
            organism.setAttribute("taxonId", String.valueOf(taxonId));
            organism.setAttribute("variety", variety); // required
            organismMap.put(organismId, organism);
        }
        if (organismMap.size()==0) {
            throw new RuntimeException("Property organisms must contain at least one taxon ID in project.xml.");
        }


        // read the reactome file
        String reactomeFilename = getChadoDBConverter().getReactomeFilename();
        LOG.info("Reading "+reactomeFilename);
        BufferedReader reader = new BufferedReader(new FileReader(reactomeFilename));
        String line = null;
        while ((line=reader.readLine())!=null) {
            String[] fields = line.split("\\\t");
            String pathwayId = fields[0];
            String pathwayName = fields[1];
            int taxonId = Integer.parseInt(fields[2]);
            String reactomeName = fields[3];
            // loop over matching chado organisms
            for (Integer organismId : organismMap.keySet()) {
                Item organism = organismMap.get(organismId);
                int orgTaxonId = Integer.parseInt(organism.getAttribute("taxonId").getValue());
                if (taxonId==orgTaxonId) {
                    // pathway
                    Item pathway = null;
                    if (pathwayMap.containsKey(pathwayName)) {
                        // get from pathwayMap
                        pathway = pathwayMap.get(pathwayName);
                    } else {
                        // put into pathwayMap
                        pathway = getChadoDBConverter().createItem("Pathway");
                        pathway.setAttribute("name", pathwayName);
                        pathwayMap.put(pathwayName, pathway);
                    }
                    // form the string to match to the protein name
                    String proteinName = reactomeName;
                    String[] parts = reactomeName.split("\\.");
                    if (parts.length>2) proteinName = parts[0]+"."+parts[1]; // drop .n assuming it's species.loc.n (e.g. Phvul.007G094900.1)
                    //// HACK START
                    // deal with some special mismatch cases
                    if (proteinName.startsWith("GLYMA")) proteinName = "Glyma" + proteinName.substring(5);
                    if (proteinName.startsWith("MTR_")) proteinName = "Medtr" + proteinName.substring(4);
                    //// HACK END
                    // query for the chado polypeptides
                    String query = "SELECT * FROM feature WHERE organism_id="+organismId+" AND type_id="+polypeptideCVTermId+" AND name LIKE '"+proteinName+"%'";
                    rs = stmt.executeQuery(query);
                    while (rs.next()) {
                        String primaryIdentifier = rs.getString("uniquename");
                        int feature_id = rs.getInt("feature_id");
                        Item protein = null;
                        if (proteinMap.containsKey((Integer)feature_id)) {
                            protein = proteinMap.get((Integer)feature_id);
                        } else {
                            protein = getChadoDBConverter().createItem("Protein");
                            protein.setAttribute("primaryIdentifier", primaryIdentifier);
                            proteinMap.put((Integer)feature_id, protein);
                        }
                        protein.setAttribute("reactomeName", reactomeName);
                        protein.addToCollection("pathways", pathway);
                    }
                    rs.close();
                }
            }
        }
                                                       
        // store all the stuff
        LOG.info("Storing "+organismMap.size()+" organisms...");
        for (Item item : organismMap.values()) store(item);
        
        LOG.info("Storing "+pathwayMap.size()+" pathways...");
        for (Item item : pathwayMap.values()) store(item);
        
        LOG.info("Storing "+proteinMap.size()+" proteins...");
        for (Item item : proteinMap.values()) store(item);

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
