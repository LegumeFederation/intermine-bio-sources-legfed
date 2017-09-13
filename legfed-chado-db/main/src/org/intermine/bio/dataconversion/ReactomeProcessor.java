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
 * pathwayId pathwayName proteinName
 * </pre>
 * The file is given by the parameter reactome.file in project.xml.
 *
 * It is assumed that the proteinName in the file will match the corresponding parts of uniquename from chado.
 * For example, a protein name like Phvul.009G083500.1 will match the chado uniquename Phvul.009G083500.1.v1.0 since the first two dot-separated parts match.
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
	Set<Integer> organismIDs = new HashSet<Integer>();           // store the desired chado organism_id values
        Map<String,Item> pathwayMap = new HashMap<String,Item>();    // keyed by name, which is all there is
        Map<Integer,Item> proteinMap = new HashMap<Integer,Item>();  // keyed by feature_id since stored from a chado query

        // initialize our DB statement
        Statement stmt = connection.createStatement();
        ResultSet rs;
        
        // build the organismIDs set from the supplied taxon IDs
        Map<Integer,OrganismData> chadoToOrgData = getChadoDBConverter().getChadoIdToOrgDataMap();
        for (Map.Entry<Integer,OrganismData> entry : chadoToOrgData.entrySet()) {
            Integer organismID = entry.getKey();
            organismIDs.add(organismID);
        }
        LOG.info("Will match against "+organismIDs.size()+" organisms in chado.");
        if (organismIDs.size()==0) {
            throw new RuntimeException("Property organisms must contain at least one taxon ID in project.xml.");
        }

        // form the SQL "IN" list for our organisms
        String inList = "(";
        boolean first = true;
        for (int organism_id : organismIDs) {
            if (first) {
                first = false;
            } else {
                inList += ",";
            }
            inList += organism_id;
        }
        inList += ")";

        // get the cv terms for our items of interest
        rs = stmt.executeQuery("SELECT cvterm_id FROM cvterm WHERE name='polypeptide'");
        rs.next();
        int polypeptideCVTermId = rs.getInt("cvterm_id");
        rs.close();

        // read the reactome file
        String reactomeFilename = getChadoDBConverter().getReactomeFilename();
        LOG.info("Reading "+reactomeFilename);
        BufferedReader reader = new BufferedReader(new FileReader(reactomeFilename));
        String line = null;
        while ((line=reader.readLine())!=null) {
            if (!line.startsWith("#")) {
                String[] fields = line.split("\t");
                String pathwayId = fields[0];
                String pathwayName = fields[1];
                String reactomeName = fields[2];
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
                if (parts.length>3) proteinName = parts[0]+"."+parts[1]+"."+parts[2]; // drop what comes after species.loc.n (e.g. Phvul.007G094900.1)
                // query for the chado polypeptides of the desired organisms
                String query = "SELECT * FROM feature WHERE type_id="+polypeptideCVTermId+" AND uniquename LIKE '"+proteinName+"%' AND organism_id IN "+inList;
                // DEBUG
                LOG.info(query);
                rs = stmt.executeQuery(query);
                while (rs.next()) {
                    String primaryIdentifier = rs.getString("uniquename");
                    int feature_id = rs.getInt("feature_id");
                    Item protein = null;
                    if (proteinMap.containsKey(feature_id)) {
                        protein = proteinMap.get(feature_id);
                    } else {
                        protein = getChadoDBConverter().createItem("Protein");
                        protein.setAttribute("primaryIdentifier", primaryIdentifier);
                        proteinMap.put(feature_id, protein);
                    }
                    protein.setAttribute("reactomeName", reactomeName);
                    protein.addToCollection("pathways", pathway);
                }
                rs.close();
            }
        }
        
        // store all the stuff
        LOG.info("Storing "+proteinMap.size()+" proteins...");
        for (Item item : proteinMap.values()) store(item);
        LOG.info("Storing "+pathwayMap.size()+" pathways...");
        for (Item item : pathwayMap.values()) store(item);
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
