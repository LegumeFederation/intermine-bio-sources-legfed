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
import java.io.File;
import java.io.FileReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import org.intermine.bio.util.OrganismData;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Attribute;
import org.intermine.xml.full.Item;
import org.intermine.xml.full.Reference;
import org.intermine.xml.full.ReferenceList;

/**
 * Store the genetic marker and QTL data from a CMap tab-delimited file. Fields are:
 *
 * map_acc map_name map_start map_stop feature_acc feature_name feature_aliases feature_start feature_stop feature_type_acc is_landmark
 *
 * @author Sam Hokin, NCGR
 */
public class CMapFileProcessor extends ChadoProcessor {
	
    private static final Logger LOG = Logger.getLogger(CMapFileProcessor.class);

    // global scope for use in methods
    Item organism;

    // store items in maps for ultimate storage, may need access in methods, so global
    ItemMap chromosomeMap = new ItemMap();

    // these are keyed by acc in file
    Map<String,Item> linkageGroupMap = new HashMap<String,Item>();
    Map<String,Item> geneticMarkerMap = new HashMap<String,Item>();
    Map<String,Item> qtlMap = new HashMap<String,Item>();

    // these are just collected in sets for final storage
    Set<Item> linkageGroupPositionSet = new HashSet<Item>();
    Set<Item> linkageGroupRangeSet = new HashSet<Item>();
    
    /**
     * Create a new CMapFileProcessor
     * @param chadoDBConverter the ChadoDBConverter that is controlling this processor
     */
    public CMapFileProcessor(ChadoDBConverter chadoDBConverter) {
        super(chadoDBConverter);
    }

    /**
     * {@inheritDoc}
     * We process the chado database by reading the feature, featureloc, featurepos, featuremap, feature_relationship and featureprop tables.
     */
    @Override
    public void process(Connection connection) throws SQLException, ObjectStoreException {

        // ---------------------------------------------------------
        // ---------------- INITIAL DATA LOADING -------------------
        // ---------------------------------------------------------

        // get the CMap file, do nothing if it's null
        File cMapFile = getChadoDBConverter().getCmapFile();
        if (cMapFile==null) {
            LOG.error("CMap file has not been set in project.xml.");
            System.exit(1);
        }

        // get chado organism_id from supplied taxon ID - enforce processing a single organism!
        Map<Integer,OrganismData> chadoToOrgData = getChadoDBConverter().getChadoIdToOrgDataMap();
        if (chadoToOrgData.size()>1) {
            System.err.println("ERROR - multiple chado organisms specified in data source; CMapFileProcessor can only process one organism at a time.");
            System.exit(1);
        }
        Integer organismId = 0;
        for (Integer key : chadoToOrgData.keySet()) {
            organismId = key.intValue();
        }
        OrganismData org = chadoToOrgData.get(organismId);
        int organism_id = organismId.intValue(); // chado organism.organism_id
        int taxonId = org.getTaxonId();

        // create organism Item - global so can be used in populate routines
        organism = getChadoDBConverter().createItem("Organism");
        BioStoreHook.setSOTerm(getChadoDBConverter(), organism, "organism", getChadoDBConverter().getSequenceOntologyRefId());
        organism.setAttribute("taxonId", String.valueOf(taxonId));

        // ---------------------------------------------------------------------------------------------------------------------------
        // Load linkage groups, genetic markers and QTLs from the CMap file and associate QTLs and genetic markers with linkage groups
        // ---------------------------------------------------------------------------------------------------------------------------

        try {
            
            BufferedReader cmapReader = new BufferedReader(new FileReader(cMapFile));
            String cmapLine = cmapReader.readLine(); // first line is header
            LOG.info("Reading CMap file:"+cMapFile.getName()+" with header:"+cmapLine);
            
            while ((cmapLine=cmapReader.readLine())!=null) {
                
                CMapRecord cmap = new CMapRecord(cmapLine);
                if (cmap.map_acc!=null) {
                
                    // add this linkage group to its map if not already in
                    Item linkageGroup = null;
                    if (linkageGroupMap.containsKey(cmap.map_acc)) {
                        linkageGroup = linkageGroupMap.get(cmap.map_acc);
                    } else {
                        linkageGroup = getChadoDBConverter().createItem("LinkageGroup");
                        populateLinkageGroup(linkageGroup, cmap);
                        linkageGroupMap.put(cmap.map_acc, linkageGroup);
                    }
                    
                    // add this QTL to this linkage group if appropriate
                    if (cmap.isQTL()) {
                        if (!qtlMap.containsKey(cmap.feature_acc)) {
                            Item qtl = getChadoDBConverter().createItem("QTL");
                            populateQTL(qtl, linkageGroup, cmap);
                            qtlMap.put(cmap.feature_acc, qtl);
                        }
                    }
                    
                    // add this genetic marker to this linkage group if appropriate
                    if (cmap.isMarker()) {
                        if (!geneticMarkerMap.containsKey(cmap.feature_acc)) {
                            Item geneticMarker = getChadoDBConverter().createItem("GeneticMarker");
                            populateGeneticMarker(geneticMarker, linkageGroup, cmap);
                            geneticMarkerMap.put(cmap.feature_acc, geneticMarker);
                        }
                    }
                        
                }
                
            }

            cmapReader.close();
            LOG.info("Created "+linkageGroupMap.size()+" LinkageGroup items.");
            LOG.info("Created "+qtlMap.size()+" QTL items.");
            LOG.info("Created "+geneticMarkerMap.size()+" GeneticMarker items.");

        } catch (Exception ex) {

            LOG.error(ex.getMessage());

        }

        // ----------------------------------------------------------------
        // we're done, so store everything
        // ----------------------------------------------------------------
        LOG.info("Storing organism for taxonId="+taxonId+".");
        store(organism);

        LOG.info("Storing linkage groups...");
        for (Item item : linkageGroupMap.values()) store(item);

        LOG.info("Storing linkage group positions...");
        for (Item item : linkageGroupPositionSet) store(item);

        LOG.info("Storing linkage group ranges...");
        for (Item item : linkageGroupRangeSet) store(item);

        LOG.info("Storing QTLs...");
        for (Item item : qtlMap.values()) store(item);

        LOG.info("Storing genetic markers...");
        for (Item item : geneticMarkerMap.values()) store(item);
 
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
     * Populate a LinkageGroup Item from a CMap record; use map_name as primary!
     */
    void populateLinkageGroup(Item linkageGroup, CMapRecord cmap) {
        BioStoreHook.setSOTerm(getChadoDBConverter(), linkageGroup, "linkage_group", getChadoDBConverter().getSequenceOntologyRefId());
        linkageGroup.setReference("organism", organism);
        linkageGroup.setAttribute("secondaryIdentifier", cmap.map_acc);
        linkageGroup.setAttribute("primaryIdentifier", cmap.map_name);
        linkageGroup.setAttribute("length", String.valueOf(cmap.map_stop));
    }

    /**
     * Populate a QTL Item from a related LinkageGroup Item and CMap record; use map_name as primary!
     */
    void populateQTL(Item qtl, Item linkageGroup, CMapRecord cmap) throws ObjectStoreException {
        BioStoreHook.setSOTerm(getChadoDBConverter(), qtl, "QTL", getChadoDBConverter().getSequenceOntologyRefId());
        qtl.setReference("organism", organism);
        qtl.setAttribute("secondaryIdentifier", cmap.feature_acc);
        qtl.setAttribute("primaryIdentifier", cmap.feature_name);
        // create and store linkage group range; place it in map as well for future processing
        Item linkageGroupRange = getChadoDBConverter().createItem("LinkageGroupRange");
        linkageGroupRange.setAttribute("begin", String.valueOf(cmap.feature_start));
        linkageGroupRange.setAttribute("end", String.valueOf(cmap.feature_stop));
        linkageGroupRange.setAttribute("length", String.valueOf(round(cmap.feature_stop-cmap.feature_start,2)));
        linkageGroupRange.setReference("linkageGroup", linkageGroup);
        linkageGroupRangeSet.add(linkageGroupRange);
        // add to QTL collection
        qtl.addToCollection("linkageGroupRanges", linkageGroupRange);
        // add to linkage group collection
        linkageGroup.addToCollection("QTLs", qtl);
    }

    /**
     * Populate a GeneticMarker Item from a LinkageGroup Item and CMap record; use map_name as primary!
     */
    void populateGeneticMarker(Item geneticMarker, Item linkageGroup, CMapRecord cmap) throws ObjectStoreException {
        BioStoreHook.setSOTerm(getChadoDBConverter(), geneticMarker, "genetic_marker", getChadoDBConverter().getSequenceOntologyRefId());
        geneticMarker.setReference("organism", organism);
        geneticMarker.setAttribute("secondaryIdentifier", cmap.feature_acc);
        geneticMarker.setAttribute("primaryIdentifier", cmap.feature_name);
        geneticMarker.setAttribute("type", cmap.feature_type_acc);
        // create and store linkage group position; place it in a map as well for future processing
        Item linkageGroupPosition = getChadoDBConverter().createItem("LinkageGroupPosition");
        linkageGroupPosition.setAttribute("position", String.valueOf(cmap.feature_start));
        linkageGroupPosition.setReference("linkageGroup", linkageGroup);
        linkageGroupPositionSet.add(linkageGroupPosition);
        geneticMarker.addToCollection("linkageGroupPositions", linkageGroupPosition);
        // add to linkage group collection
        linkageGroup.addToCollection("geneticMarkers", geneticMarker);
    }

    /**
     * Round a double to the given number of places
     */
    public static double round(double value, int places) {
        if (places < 0) throw new IllegalArgumentException();
        BigDecimal bd = new BigDecimal(value);
        bd = bd.setScale(places, RoundingMode.HALF_UP);
        return bd.doubleValue();
    }


}
