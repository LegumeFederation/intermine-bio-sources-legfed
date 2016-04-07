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
        int qtlCVTermId = getCVTermId(stmt, "QTL");
        int geneticMarkerCVTermId = getCVTermId(stmt, "genetic_marker");
        int linkageGroupCVTermId = getCVTermId(stmt, "linkage_group");
        int polypeptideCVTermId = getCVTermId(stmt, "polypeptide");
        int proteinMatchCVTermId = getCVTermId(stmt, "protein_match");
        int proteinHmmMatchCVTermId = getCVTermId(stmt, "protein_hmm_match");
        int geneCVTermId = getCVTermId(stmt, "gene");

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
        
        // build the Organism map from the supplied taxon IDs
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

        // store the Items in maps keyed by feature_id
        Map<Integer,Item> geneMap = new HashMap<Integer,Item>();
        Map<Integer,Item> linkageGroupMap = new HashMap<Integer,Item>();
        Map<Integer,Item> geneticMarkerMap = new HashMap<Integer,Item>();
        Map<Integer,Item> qtlMap = new HashMap<Integer,Item>();
        Map<Integer,Item> polypeptideMap = new HashMap<Integer,Item>();
        Map<Integer,Item> proteinMatchMap = new HashMap<Integer,Item>();
        Map<Integer,Item> proteinHmmMatchMap = new HashMap<Integer,Item>();

        // loop over the requested organisms and load the maps
        for (Map.Entry<Integer,Item> orgEntry : organismMap.entrySet()) {
            int organism_id = orgEntry.getKey().intValue();
            Item organism = orgEntry.getValue();
            populateMap(stmt, geneMap, "Gene", geneCVTermId, organism_id, organism);
            populateMap(stmt, polypeptideMap, "Polypeptide", polypeptideCVTermId, organism_id, organism);
            populateMap(stmt, proteinMatchMap, "ProteinMatch", proteinMatchCVTermId, organism_id, organism);
            populateMap(stmt, proteinHmmMatchMap, "ProteinHmmMatch", proteinHmmMatchCVTermId, organism_id, organism);
            populateMap(stmt, linkageGroupMap, "LinkageGroup", linkageGroupCVTermId, organism_id, organism);
            populateMap(stmt, geneticMarkerMap, "GeneticMarker", geneticMarkerCVTermId, organism_id, organism);
            populateMap(stmt, qtlMap, "QTL", qtlCVTermId, organism_id, organism);
        }
        LOG.info("Created "+geneMap.size()+" Gene items.");
        LOG.info("Created "+polypeptideMap.size()+" Polypeptide items.");
        LOG.info("Created "+proteinMatchMap.size()+" ProteinMatch items.");
        LOG.info("Created "+proteinHmmMatchMap.size()+" ProteinHmmMatch items.");
        LOG.info("Created "+linkageGroupMap.size()+" LinkageGroup items.");
        LOG.info("Created "+geneticMarkerMap.size()+" GeneticMarker items.");
        LOG.info("Created "+qtlMap.size()+" QTL items.");

        // load gene attributes
        LOG.info("Loading gene attributes...");
        for (Map.Entry<Integer,Item> entry : geneMap.entrySet()) {
            int feature_id = entry.getKey().intValue();
            Item gene = entry.getValue();
            loadAttribute(stmt, feature_id, gene, "note", noteCVTermId);
        }

        // load QTL attributes
        LOG.info("Loading QTL attributes...");
        for (Map.Entry<Integer,Item> entry : qtlMap.entrySet()) {
            int feature_id = entry.getKey().intValue();
            Item qtl = entry.getValue();
            loadAttribute(stmt, feature_id, qtl, "comment", commentCVTermId);
            loadAttribute(stmt, feature_id, qtl, "experimentTraitDescription", experimentTraitDescriptionCVTermId);
            loadAttribute(stmt, feature_id, qtl, "experimentTraitName", experimentTraitNameCVTermId);
            loadAttribute(stmt, feature_id, qtl, "publicationLinkageGroup", publicationLinkageGroupCVTermId);
            loadAttribute(stmt, feature_id, qtl, "analysisMethod", qtlAnalysisMethodCVTermId);
            loadAttribute(stmt, feature_id, qtl, "traitUnit", traitUnitCVTermId);
            loadAttribute(stmt, feature_id, qtl, "identifier", qtlIdentifierCVTermId);
            loadAttribute(stmt, feature_id, qtl, "peak", qtlPeakCVTermId);
            loadAttribute(stmt, feature_id, qtl, "studyTreatment", qtlStudyTreatmentCVTermId);
        }

        // load genetic marker attributes
        LOG.info("Loading genetic marker attributes...");
        for (Map.Entry<Integer,Item> entry : geneticMarkerMap.entrySet()) {
            int feature_id = entry.getKey().intValue();
            Item geneticMarker = entry.getValue();
            loadAttribute(stmt, feature_id, geneticMarker, "comment", commentCVTermId);
            loadAttribute(stmt, feature_id, geneticMarker, "sourceDescription", sourceDescriptionCVTermId);
            loadAttribute(stmt, feature_id, geneticMarker, "canonicalMarker", canonicalMarkerCVTermId);
        }

        // load linkage group attributes
        LOG.info("Loading linkage group attributes...");
        for (Map.Entry<Integer,Item> entry : linkageGroupMap.entrySet()) {
            int feature_id = entry.getKey().intValue();
            Item linkageGroup = entry.getValue();
            loadAttribute(stmt, feature_id, linkageGroup, "assignedLinkageGroup", assignedLinkageGroupCVTermId);
        }

        // load polypeptide attributes
        LOG.info("Loading polypeptide attributes...");
        for (Map.Entry<Integer,Item> entry : polypeptideMap.entrySet()) {
            int feature_id = entry.getKey().intValue();
            Item polypeptide = entry.getValue();
            loadAttribute(stmt, feature_id, polypeptide, "note", noteCVTermId);
        }

        // load protein match attributes
        LOG.info("Loading protein match attributes...");
        for (Map.Entry<Integer,Item> entry : proteinMatchMap.entrySet()) {
            int feature_id = entry.getKey().intValue();
            Item proteinMatch = entry.getValue();
            loadAttribute(stmt, feature_id, proteinMatch, "signatureDesc", signatureDescCVTermId);
            loadAttribute(stmt, feature_id, proteinMatch, "status", statusCVTermId);
            loadAttribute(stmt, feature_id, proteinMatch, "date", dateCVTermId);
        }

        // load protein HMM match attributes
        LOG.info("Loading protein HMM match attributes...");
        for (Map.Entry<Integer,Item> entry : proteinHmmMatchMap.entrySet()) {
            int feature_id = entry.getKey().intValue();
            Item proteinHmmMatch = entry.getValue();
            loadAttribute(stmt, feature_id, proteinHmmMatch, "signatureDesc", signatureDescCVTermId);
            loadAttribute(stmt, feature_id, proteinHmmMatch, "status", statusCVTermId);
            loadAttribute(stmt, feature_id, proteinHmmMatch, "date", dateCVTermId);
        }

        // store our Items
        LOG.info("Storing genes...");
        for (Item item : geneMap.values()) store(item);

        LOG.info("Storing linkage groups...");
        for (Item item : linkageGroupMap.values()) store(item);

        LOG.info("Storing genetic markers...");
        for (Item item : geneticMarkerMap.values()) store(item);

        LOG.info("Storing QTLs...");
        for (Item item : qtlMap.values()) store(item);

        LOG.info("Storing polypeptides...");
        for (Item item : polypeptideMap.values()) store(item);

        LOG.info("Storing protein matches...");
        for (Item item : proteinMatchMap.values()) store(item);

        LOG.info("Storing protein HMM matches...");
        for (Item item : proteinHmmMatchMap.values()) store(item);

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
     * Populate a map with records from the feature table
     */
    void populateMap(Statement stmt, Map<Integer,Item> map, String className, int cvTermId, int organism_id, Item organism) throws SQLException {
        ResultSet rs = stmt.executeQuery("SELECT * FROM feature WHERE organism_id="+organism_id+" AND type_id="+cvTermId);
        while (rs.next()) {
            ChadoFeature cf = new ChadoFeature(rs);
            Item item = getChadoDBConverter().createItem(className);
            cf.populateBioEntity(item, organism);
            map.put(new Integer(cf.feature_id), item);
        }
        rs.close();
    }

    /**
     * Load an attribute from the featureprop table
     */
    void loadAttribute(Statement stmt, int feature_id, Item item, String attributeName, int cvTermId) throws SQLException {
        ResultSet rs = stmt.executeQuery("SELECT * FROM featureprop WHERE feature_id="+feature_id+" AND type_id="+cvTermId+" ORDER BY rank");
        if (rs.next() && rs.getString("value")!=null && rs.getString("value").length()>0) {
            item.setAttribute(attributeName, rs.getString("value"));
        }
        rs.close();
    }

}
