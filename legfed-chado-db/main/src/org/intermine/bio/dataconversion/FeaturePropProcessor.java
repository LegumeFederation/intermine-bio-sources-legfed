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
        
        // feature types of interest
        int geneCVTermId = getCVTermId(stmt, "gene");
        int geneticMarkerCVTermId = getCVTermId(stmt, "genetic_marker");
        int linkageGroupCVTermId = getCVTermId(stmt, "linkage_group");
        int polypeptideCVTermId = getCVTermId(stmt, "polypeptide");
        int proteinMatchCVTermId = getCVTermId(stmt, "protein_match");
        int proteinHmmMatchCVTermId = getCVTermId(stmt, "protein_hmm_match");
        int qtlCVTermId = getCVTermId(stmt, "QTL");

        // desired featureprop types
        int noteCVTermId = getCVTermId(stmt, "Note");
        int commentCVTermId = getCVTermId(stmt, "comment");
        int experimentTraitDescriptionCVTermId = getCVTermId(stmt, "Experiment Trait Description");
        int experimentTraitNameCVTermId = getCVTermId(stmt, "Experiment Trait Name");
        int publicationLinkageGroupCVTermId = getCVTermId(stmt, "Publication Linkage Group");
        int traitUnitCVTermId = getCVTermId(stmt, "Trait Unit");
        int qtlAnalysisMethodCVTermId = getCVTermId(stmt, "QTL Analysis Method");
        int qtlIdentifierCVTermId = getCVTermId(stmt, "QTL Identifier");
        int qtlPeakCVTermId = getCVTermId(stmt, "QTL Peak");
        int qtlStudyTreatmentCVTermId = getCVTermId(stmt, "QTL Study Treatment");
        int sourceDescriptionCVTermId = getCVTermId(stmt, "Source Description");
        int canonicalMarkerCVTermId = getCVTermId(stmt, "Canonical Marker");
        int assignedLinkageGroupCVTermId = getCVTermId(stmt, "Assigned Linkage Group");
        int signatureDescCVTermId = getCVTermId(stmt, "signature_desc");
        int statusCVTermId = getCVTermId(stmt, "status");
        int dateCVTermId = getCVTermId(stmt, "date");
        
        // loop over the requested organisms
        Map<Integer,OrganismData> chadoToOrgData = getChadoDBConverter().getChadoIdToOrgDataMap();
        for (Map.Entry<Integer,OrganismData> entry : chadoToOrgData.entrySet()) {
            
            // create and store organism Item
            Integer organismId = entry.getKey();
            int organism_id = organismId.intValue();
            OrganismData organismData = entry.getValue();
            int taxonId = organismData.getTaxonId();
            LOG.info("Initializing organism Item...");
            Item organism = getChadoDBConverter().createItem("Organism");
            organism.setAttribute("taxonId", String.valueOf(taxonId));
            // add organism.comment if it exists
            ResultSet rs = stmt.executeQuery("SELECT * FROM organism WHERE organism_id="+organism_id);
            if (rs.next()) {
                String comment = rs.getString("comment");
                if (comment!=null && comment.trim().length()>0) {
                    LOG.info("Adding organism.comment");
                    organism.setAttribute("comment", comment.trim());
                }
            }
            rs.close();
            LOG.info("Storing organism "+taxonId+"...");
            store(organism);

            // load gene attributes
            Map<Integer,Item> geneMap = generateMap(stmt, "Gene", geneCVTermId, organism_id, organism);
            if (geneMap.size()>0) {
                LOG.info("Loading attributes for "+geneMap.size()+" Gene items for organism "+taxonId);
                loadAttributes(stmt, organism_id, geneMap, geneCVTermId, "note", noteCVTermId);
                LOG.info("Storing genes for organism "+taxonId);
                for (Item item : geneMap.values()) store(item);
            }
            
            // load linkage group attributes
            Map<Integer,Item> linkageGroupMap = generateMap(stmt, "LinkageGroup", linkageGroupCVTermId, organism_id, organism);
            if (linkageGroupMap.size()>0) {
                LOG.info("Loading attributes for "+linkageGroupMap.size()+" LinkageGroup items for organism "+taxonId);
                loadAttributes(stmt, organism_id, linkageGroupMap, linkageGroupCVTermId, "assignedLinkageGroup", assignedLinkageGroupCVTermId);
                LOG.info("Storing linkage groups for organism "+taxonId);
                for (Item item : linkageGroupMap.values()) store(item);
            }

            // load QTL attributes
            Map<Integer,Item> qtlMap = generateMap(stmt, "QTL", qtlCVTermId, organism_id, organism);
            if (qtlMap.size()>0) {
                LOG.info("Loading attributes for "+qtlMap.size()+" QTL items for organism "+taxonId);
                loadAttributes(stmt, organism_id, qtlMap, qtlCVTermId, "comment", commentCVTermId);
                loadAttributes(stmt, organism_id, qtlMap, qtlCVTermId, "experimentTraitDescription", experimentTraitDescriptionCVTermId);
                loadAttributes(stmt, organism_id, qtlMap, qtlCVTermId, "experimentTraitName", experimentTraitNameCVTermId);
                loadAttributes(stmt, organism_id, qtlMap, qtlCVTermId, "publicationLinkageGroup", publicationLinkageGroupCVTermId);
                loadAttributes(stmt, organism_id, qtlMap, qtlCVTermId, "analysisMethod", qtlAnalysisMethodCVTermId);
                loadAttributes(stmt, organism_id, qtlMap, qtlCVTermId, "traitUnit", traitUnitCVTermId);
                loadAttributes(stmt, organism_id, qtlMap, qtlCVTermId, "identifier", qtlIdentifierCVTermId);
                loadAttributes(stmt, organism_id, qtlMap, qtlCVTermId, "peak", qtlPeakCVTermId);
                loadAttributes(stmt, organism_id, qtlMap, qtlCVTermId, "studyTreatment", qtlStudyTreatmentCVTermId);
                LOG.info("Storing QTLs for organism "+taxonId);
                for (Item item : qtlMap.values()) store(item);
            }

            // load genetic marker attributes
            Map<Integer,Item> geneticMarkerMap = generateMap(stmt, "GeneticMarker", geneticMarkerCVTermId, organism_id, organism);
            if (geneticMarkerMap.size()>0) {
                LOG.info("Loading attributes for "+geneticMarkerMap.size()+" GeneticMarker items for organism "+taxonId);
                loadAttributes(stmt, organism_id, geneticMarkerMap, geneticMarkerCVTermId, "comment", commentCVTermId);
                loadAttributes(stmt, organism_id, geneticMarkerMap, geneticMarkerCVTermId, "sourceDescription", sourceDescriptionCVTermId);
                loadAttributes(stmt, organism_id, geneticMarkerMap, geneticMarkerCVTermId, "canonicalMarker", canonicalMarkerCVTermId);
                LOG.info("Storing genetic markers for organism "+taxonId);
                for (Item item : geneticMarkerMap.values()) store(item);
            }

            // load polypeptide attributes
            Map<Integer,Item> polypeptideMap = generateMap(stmt, "Polypeptide", polypeptideCVTermId, organism_id, organism);
            if (polypeptideMap.size()>0) {
                LOG.info("Loading attributes for "+polypeptideMap.size()+" Polypeptide items for organism "+taxonId);
                loadAttributes(stmt, organism_id, polypeptideMap, polypeptideCVTermId, "note", noteCVTermId);
                LOG.info("Storing polypeptides for organism "+taxonId);
                for (Item item : polypeptideMap.values()) store(item);
            }

            // load protein match attributes
            Map<Integer,Item> proteinMatchMap = generateMap(stmt, "ProteinMatch", proteinMatchCVTermId, organism_id, organism);
            if (proteinMatchMap.size()>0) {
                LOG.info("Loading attributes for "+proteinMatchMap.size()+" ProteinMatch items for organism "+taxonId);
                loadAttributes(stmt, organism_id, proteinMatchMap, proteinMatchCVTermId, "signatureDesc", signatureDescCVTermId);
                loadAttributes(stmt, organism_id, proteinMatchMap, proteinMatchCVTermId, "status", statusCVTermId);
                loadAttributes(stmt, organism_id, proteinMatchMap, proteinMatchCVTermId, "date", dateCVTermId);
                LOG.info("Storing protein matches for organism "+taxonId);
                for (Item item : proteinMatchMap.values()) store(item);
            }

            // load protein HMM match attributes
            Map<Integer,Item> proteinHmmMatchMap = generateMap(stmt, "ProteinHmmMatch", proteinHmmMatchCVTermId, organism_id, organism);
            if (proteinHmmMatchMap.size()>0) {
                LOG.info("Loading attributes for "+proteinHmmMatchMap.size()+" ProteinHmmMatch items for organism "+taxonId);
                loadAttributes(stmt, organism_id, proteinHmmMatchMap, proteinHmmMatchCVTermId, "signatureDesc", signatureDescCVTermId);
                loadAttributes(stmt, organism_id, proteinHmmMatchMap, proteinHmmMatchCVTermId, "status", statusCVTermId);
                loadAttributes(stmt, organism_id, proteinHmmMatchMap, proteinHmmMatchCVTermId, "date", dateCVTermId);
                LOG.info("Storing protein HMM matches for organism "+taxonId);
                for (Item item : proteinHmmMatchMap.values()) store(item);
            }

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


    /**
     * Return the cvterm.cvterm_id for the requested cvterm.name, 0 if not found.
     * Require that a non-empty definition exist - this narrows some dupes (like "comment") to a single record.
     */
    int getCVTermId(Statement stmt, String name) throws SQLException {
        ResultSet rs = stmt.executeQuery("SELECT * FROM cvterm WHERE name='"+name+"' AND length(definition)>0");
        int cvTermId = 0;
        if (rs.next()) cvTermId = rs.getInt("cvterm_id");
        rs.close();
        return cvTermId;
    }

    /**
     * Create a map with records from the feature table; restrict to feature records which have corresponding records in featureprop.
     */
    Map<Integer,Item> generateMap(Statement stmt, String className, int cvTermId, int organism_id, Item organism) throws SQLException {
        Map<Integer,Item> map = new HashMap<Integer,Item>();
        ResultSet rs = stmt.executeQuery("SELECT DISTINCT feature.* FROM feature,featureprop WHERE feature.feature_id=featureprop.feature_id AND organism_id="+organism_id+" AND feature.type_id="+cvTermId);
        while (rs.next()) {
            Integer featureId = new Integer(rs.getInt("feature_id"));
            ChadoFeature cf = new ChadoFeature(rs);
            Item item = getChadoDBConverter().createItem(className);
            cf.populateBioEntity(item, organism);
            map.put(featureId, item);
        }
        rs.close();
        return map;
    }

    /**
     * Load the attributes from the featureprop table for a given feature type and featureprop type
     */
    void loadAttributes(Statement stmt, int organism_id, Map<Integer,Item> itemMap, int itemCVTermId, String attributeName, int propCVTermId) throws SQLException {
        LOG.info("loadAttributes: attributeName="+attributeName+", feature.type_id="+itemCVTermId+", organism_id="+organism_id);
        ResultSet rs = stmt.executeQuery("SELECT featureprop.* FROM featureprop,feature " +
                                         "WHERE featureprop.feature_id=feature.feature_id " +
                                         "AND featureprop.type_id="+propCVTermId+" " +
                                         "AND feature.organism_id="+organism_id+" " +
                                         "AND feature.type_id="+itemCVTermId+" " +
                                         "ORDER BY feature_id,type_id,rank DESC");
        while (rs.next()) {
            Integer featureId = new Integer(rs.getInt("feature_id"));
            String value = rs.getString("value");
            if (value!=null && value.trim().length()>0) {
                Item item = itemMap.get(featureId);
                if (item!=null) {
                    item.setAttribute(attributeName, value.trim());
                } else {
                    LOG.error("Item for feature_id="+featureId+" is null.");
                }
            }
        }
        rs.close();
    }

}
