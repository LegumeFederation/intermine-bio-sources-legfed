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

import java.util.Collection;
import java.util.HashMap;
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
 * Store the extra data from the chado.featureprop table in attributes defined in the data model.
 * It does not create any new items; these must be merged with items that have already been loaded from chado.
 *
 * Since this processer deals only with chado data, Items are stored in maps with Integer keys equal to
 * the chado feature.feature_id.
 *
 * @author Sam Hokin, NCGR
 */
public class FeaturePropProcessor extends ChadoProcessor {
	
    private static final Logger LOG = Logger.getLogger(FeaturePropProcessor.class);

    /**
     * Create a new FeaturePropProcessor
     * @param chadoDBConverter the ChadoDBConverter that is controlling this processor
     */
    public FeaturePropProcessor(ChadoDBConverter chadoDBConverter) {
        super(chadoDBConverter);
    }

    /**
     * {@inheritDoc}
     * We process the chado database by reading the featureprop table for the feature_id values of our features of interest.
     */
    @Override
    public void process(Connection connection) throws SQLException, ObjectStoreException {
        
        // initialize our DB statement
        Statement stmt = connection.createStatement();
                
        // loop over the requested organisms
        Map<Integer,OrganismData> chadoToOrgData = getChadoDBConverter().getChadoIdToOrgDataMap();
        for (Integer organismId : chadoToOrgData.keySet()) {

            int organism_id = organismId.intValue();
            OrganismData organismData = chadoToOrgData.get(organismId);
            
            // create and store organism Item
            String taxonId = organismData.getTaxonId();
            String variety = organismData.getVariety();
            Item organism = getChadoDBConverter().createItem("Organism");
            organism.setAttribute("taxonId", String.valueOf(taxonId));
            organism.setAttribute("variety", variety); // non-null variety is required
            // add organism.description if it exists
            ResultSet rs = stmt.executeQuery("SELECT * FROM organism WHERE organism_id="+organism_id);
            if (rs.next()) {
                String description = rs.getString("comment");
                if (description!=null && description.trim().length()>0) {
                    organism.setAttribute("description", description.trim());
                }
            }
            rs.close();
            LOG.info("Storing organism "+taxonId+" ("+variety+")");
            store(organism);

            // load gene attributes
            Map<Integer,Item> geneMap = generateMap(stmt, organism_id, organism, "Gene", "gene");
            geneMap = loadAttributes(geneMap, stmt, organism_id, "Gene", "description", "gene", "Note");
            store(geneMap.values());
            LOG.info("Stored "+geneMap.size()+" genes for organism "+taxonId+" ("+variety+")");


            // load linkage group attributes
            Map<Integer,Item> linkageGroupMap = generateMap(stmt, organism_id, organism, "LinkageGroup", "linkage_group");
            linkageGroupMap = loadAttributes(linkageGroupMap, stmt, organism_id, "LinkageGroup", "assignedLinkageGroup", "linkage_group", "Assigned Linkage Group");
            store(linkageGroupMap.values());
            LOG.info("Stored "+linkageGroupMap.size()+" linkage groups for organism "+taxonId+" ("+variety+")");

            // load QTL attributes
            Map<Integer,Item> qtlMap = generateMap(stmt, organism_id, organism, "QTL", "QTL");
            qtlMap = loadAttributes(qtlMap, stmt, organism_id, "QTL", "description", "QTL", "comment");
            qtlMap = loadAttributes(qtlMap, stmt, organism_id, "QTL", "traitDescription", "QTL", "Experiment Trait Description");
            qtlMap = loadAttributes(qtlMap, stmt, organism_id, "QTL", "traitName", "QTL", "Experiment Trait Name");
            qtlMap = loadAttributes(qtlMap, stmt, organism_id, "QTL", "publicationLinkageGroup", "QTL", "Publication Linkage Group");
            qtlMap = loadAttributes(qtlMap, stmt, organism_id, "QTL", "analysisMethod", "QTL", "QTL Analysis Method");
            qtlMap = loadAttributes(qtlMap, stmt, organism_id, "QTL", "traitUnit", "QTL", "Trait Unit");
            qtlMap = loadAttributes(qtlMap, stmt, organism_id, "QTL", "identifier", "QTL", "QTL Identifier");
            qtlMap = loadAttributes(qtlMap, stmt, organism_id, "QTL", "peak", "QTL", "QTL Peak");
            qtlMap = loadAttributes(qtlMap, stmt, organism_id, "QTL", "studyTreatment", "QTL", "QTL Study Treatment");
            store(qtlMap.values());
            LOG.info("Stored "+qtlMap.size()+" QTLs for organism "+taxonId+" ("+variety+")");

            // load genetic marker attributes
            Map<Integer,Item> geneticMarkerMap = generateMap(stmt, organism_id, organism, "GeneticMarker", "genetic_marker");
            geneticMarkerMap = loadAttributes(geneticMarkerMap, stmt, organism_id, "GeneticMarker", "description", "genetic_marker", "description");
            geneticMarkerMap = loadAttributes(geneticMarkerMap, stmt, organism_id, "GeneticMarker", "sourceDescription", "genetic_marker", "Source Description");
            geneticMarkerMap = loadAttributes(geneticMarkerMap, stmt, organism_id, "GeneticMarker", "canonicalMarker", "genetic_marker", "Canonical Marker");
            store(geneticMarkerMap.values());
            LOG.info("Stored "+geneticMarkerMap.size()+" genetic markers for organism "+taxonId+" ("+variety+")");

            // load protein attributes
            // NOTE: proteins are stored in chado as "polypeptide"
            Map<Integer,Item> proteinMap = generateMap(stmt, organism_id, organism, "Protein", "polypeptide");
            proteinMap = loadAttributes(proteinMap, stmt, organism_id, "Protein", "note", "polypeptide", "Note");
            store(proteinMap.values());
            LOG.info("Stored "+proteinMap.size()+" proteins for organism "+taxonId+" ("+variety+")");

            // load protein match attributes
            Map<Integer,Item> proteinMatchMap = generateMap(stmt, organism_id, organism, "ProteinMatch", "protein_match");
            proteinMatchMap = loadAttributes(proteinMatchMap, stmt, organism_id, "ProteinMatch", "signatureDesc", "protein_match", "signature_desc");
            proteinMatchMap = loadAttributes(proteinMatchMap, stmt, organism_id, "ProteinMatch", "status", "protein_match", "status");
            proteinMatchMap = loadAttributes(proteinMatchMap, stmt, organism_id, "ProteinMatch", "date", "protein_match", "date");
            store(proteinMatchMap.values());
            LOG.info("Stored "+proteinMatchMap.size()+" protein matches for organism "+taxonId+" ("+variety+")");

            // load protein HMM match attributes
            Map<Integer,Item> proteinHmmMatchMap = generateMap(stmt, organism_id, organism, "ProteinHmmMatch", "protein_hmm_match");
            proteinHmmMatchMap = loadAttributes(proteinHmmMatchMap, stmt, organism_id, "ProteinMatch", "signatureDesc", "protein_hmm_match", "signature_desc");
            proteinHmmMatchMap = loadAttributes(proteinHmmMatchMap, stmt, organism_id, "ProteinMatch", "status", "protein_hmm_match", "status");
            proteinHmmMatchMap = loadAttributes(proteinHmmMatchMap, stmt, organism_id, "ProteinMatch", "date", "protein_hmm_match", "date");
            store(proteinHmmMatchMap.values());
            LOG.info("Stored "+proteinHmmMatchMap.size()+" protein matches for organism "+taxonId+" ("+variety+")");

        } // organism loop

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
     * Store the items.
     * @param items the Item Collection
     * @throws ObjectStoreException if an error occurs while storing
     */
    protected void store(Collection<Item> items) throws ObjectStoreException {
        for (Item item : items) store(item);
        return;
    }

    /**
     * Return the cvterm.cvterm_id for the requested cvterm.name, 0 if not found.
     * Require that a non-empty definition exist - this narrows some dupes (like "description") to a single record.
     */
    int getCVTermId(Statement stmt, String name) throws SQLException {
        ResultSet rs = stmt.executeQuery("SELECT * FROM cvterm WHERE name='"+name+"' AND length(definition)>0");
        int cvTermId = 0;
        if (rs.next()) cvTermId = rs.getInt("cvterm_id");
        rs.close();
        return cvTermId;
    }

    /**
     * Create a map, keyed with chado feature_id, with records from the feature table; restrict to feature records which have corresponding records in featureprop.
     */
    Map<Integer,Item> generateMap(Statement stmt, int organism_id, Item organism, String className, String featureCVTermName) throws SQLException {
        
        int featureCVTermId = getCVTermId(stmt, featureCVTermName);

        // first load the map of items keyed by chado.feature_id and loaded with uniquename=primaryIdentifier and organism reference
        Map<Integer,Item> map = new HashMap<Integer,Item>();
        ResultSet rs = stmt.executeQuery("SELECT DISTINCT feature.feature_id,feature.uniquename " +
                                         "FROM feature,featureprop " +
                                         "WHERE feature.feature_id=featureprop.feature_id " +
                                         "AND organism_id="+organism_id+" " +
                                         "AND feature.type_id="+featureCVTermId+" " +
                                         "ORDER BY feature_id");
        while (rs.next()) {
            Integer featureId = new Integer(rs.getInt("feature_id"));
            String uniquename = rs.getString("uniquename");
            Item item = getChadoDBConverter().createItem(className);
            item.setAttribute("primaryIdentifier", uniquename);
            item.setReference("organism", organism);
            map.put(featureId, item);
        }
        rs.close();
        
        // return the map, ready to be loaded with values
        return map;
        
    }

    /**
     * Load an attribute into the items in the given map.
     */
    Map<Integer,Item> loadAttributes(Map<Integer,Item> map, Statement stmt, int organism_id, String className, String attributeName, String featureCVTermName, String propCVTermName) throws SQLException {
    
        int featureCVTermId = getCVTermId(stmt, featureCVTermName);
        int propCVTermId = getCVTermId(stmt, propCVTermName);

        // spin through the featureprop records and load the desired property value as the attribute nam
        ResultSet rs = stmt.executeQuery("SELECT featureprop.* FROM featureprop,feature " +
                                         "WHERE featureprop.feature_id=feature.feature_id " +
                                         "AND featureprop.type_id="+propCVTermId+" " +
                                         "AND feature.organism_id="+organism_id+" " +
                                         "AND feature.type_id="+featureCVTermId+" " +
                                         "ORDER BY feature_id ASC, rank DESC");
        while (rs.next()) {
            Integer featureId = new Integer(rs.getInt("feature_id"));
            if (map.containsKey(featureId)) {
                String value = rs.getString("value");
                if (value!=null && value.trim().length()>0) {
                    Item item = map.get(featureId);
                    item.setAttribute(attributeName, value.trim());
                }
            } else {
                LOG.error("No feature found for featureprop.feature_id="+featureId+".");
            }
        }
        rs.close();
        
        // return map without storing its items
        return map;
    }

}
