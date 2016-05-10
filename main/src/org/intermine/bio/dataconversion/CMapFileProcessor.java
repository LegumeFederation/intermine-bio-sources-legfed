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
import java.io.Reader;

import java.math.BigDecimal;
import java.math.RoundingMode;

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
public class CMapFileProcessor extends LegfedFileProcessor {
	
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
     * @param legfedFileConverter the LegfedFileConverter that is controlling this processor
     */
    public CMapFileProcessor(LegfedFileConverter legfedFileConverter) {
        super(legfedFileConverter);
    }

    /**
     * {@inheritDoc}
     * We process the chado database by reading the feature, featureloc, featurepos, featuremap, feature_relationship and featureprop tables.
     */
    @Override
    public void process(Reader reader) throws ObjectStoreException {

        // ---------------------------------------------------------
        // ---------------- INITIAL DATA LOADING -------------------
        // ---------------------------------------------------------

        // get the organism taxon ID; enforce only one
        OrganismData[] organisms = getLegfedFileConverter().getOrganismsToProcess().toArray(new OrganismData[0]);
        if (organisms.length>1) {
            String error = "Multiple organisms specified in data source; GeneticMarkerGFFProcessor can only process one organism at a time.";
            LOG.error(error);
            throw new RuntimeException(error);
        }
        int taxonId = organisms[0].getTaxonId();

        // create organism Item - global so can be used in populate routines
        organism = getLegfedFileConverter().createItem("Organism");
        BioStoreHook.setSOTerm(getLegfedFileConverter(), organism, "organism", getLegfedFileConverter().getSequenceOntologyRefId());
        organism.setAttribute("taxonId", String.valueOf(taxonId));

        // ---------------------------------------------------------------------------------------------------------------------------
        // Load linkage groups, genetic markers and QTLs from the CMap file and associate QTLs and genetic markers with linkage groups
        // ---------------------------------------------------------------------------------------------------------------------------

        try {
            
            BufferedReader cmapReader = new BufferedReader(reader);
            String cmapLine = cmapReader.readLine(); // first line is header
            LOG.info("Reading CMap file with header:"+cmapLine);
            
            while ((cmapLine=cmapReader.readLine())!=null) {
                
                CMapRecord cmap = new CMapRecord(cmapLine);
                if (cmap.map_acc!=null) {
                
                    // add this linkage group to its map if not already in
                    Item linkageGroup = null;
                    if (linkageGroupMap.containsKey(cmap.map_acc)) {
                        linkageGroup = linkageGroupMap.get(cmap.map_acc);
                    } else {
                        linkageGroup = getLegfedFileConverter().createItem("LinkageGroup");
                        populateLinkageGroup(linkageGroup, cmap);
                        linkageGroupMap.put(cmap.map_acc, linkageGroup);
                    }
                    
                    // add this QTL to this linkage group if appropriate
                    if (cmap.isQTL()) {
                        if (!qtlMap.containsKey(cmap.feature_acc)) {
                            Item qtl = getLegfedFileConverter().createItem("QTL");
                            populateQTL(qtl, linkageGroup, cmap);
                            qtlMap.put(cmap.feature_acc, qtl);
                        }
                    }
                    
                    // add this genetic marker to this linkage group if appropriate
                    if (cmap.isMarker()) {
                        if (!geneticMarkerMap.containsKey(cmap.feature_acc)) {
                            Item geneticMarker = getLegfedFileConverter().createItem("GeneticMarker");
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
     * Populate a LinkageGroup Item from a CMap record; use map_name as primary!
     */
    void populateLinkageGroup(Item linkageGroup, CMapRecord cmap) {
        BioStoreHook.setSOTerm(getLegfedFileConverter(), linkageGroup, "linkage_group", getLegfedFileConverter().getSequenceOntologyRefId());
        linkageGroup.setReference("organism", organism);
        linkageGroup.setAttribute("secondaryIdentifier", cmap.map_acc);
        linkageGroup.setAttribute("primaryIdentifier", cmap.map_name);
        linkageGroup.setAttribute("length", String.valueOf(cmap.map_stop));
    }

    /**
     * Populate a QTL Item from a related LinkageGroup Item and CMap record; use map_name as primary!
     */
    void populateQTL(Item qtl, Item linkageGroup, CMapRecord cmap) throws ObjectStoreException {
        BioStoreHook.setSOTerm(getLegfedFileConverter(), qtl, "QTL", getLegfedFileConverter().getSequenceOntologyRefId());
        qtl.setReference("organism", organism);
        qtl.setAttribute("secondaryIdentifier", cmap.feature_acc);
        qtl.setAttribute("primaryIdentifier", cmap.feature_name);
        // create and store linkage group range; place it in map as well for future processing
        Item linkageGroupRange = getLegfedFileConverter().createItem("LinkageGroupRange");
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
        BioStoreHook.setSOTerm(getLegfedFileConverter(), geneticMarker, "genetic_marker", getLegfedFileConverter().getSequenceOntologyRefId());
        geneticMarker.setReference("organism", organism);
        geneticMarker.setAttribute("secondaryIdentifier", cmap.feature_acc);
        geneticMarker.setAttribute("primaryIdentifier", cmap.feature_name);
        geneticMarker.setAttribute("type", cmap.feature_type_acc);
        // create and store linkage group position; place it in a map as well for future processing
        Item linkageGroupPosition = getLegfedFileConverter().createItem("LinkageGroupPosition");
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
