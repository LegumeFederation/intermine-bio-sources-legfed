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

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

/**
 * Store the genetic marker and QTL data from a CMap tab-delimited file. Fields are:
 * <pre>
 * map_acc map_name map_start map_stop feature_acc feature_name feature_aliases feature_start feature_stop feature_type_acc is_landmark
 * </pre>
 * The file name gives the taxon ID of the organism. This allows you to have several files in a single src.data.dir 
 * that are run in a single invocation of this converter. The format is: anything-other-than-underscore_3847.txt.
 *
 * @author Sam Hokin
 */
public class CMapFileConverter extends BioFileConverter {
	
    private static final Logger LOG = Logger.getLogger(CMapFileConverter.class);
    
    /**
     * Create a new CMapFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public CMapFileConverter(ItemWriter writer, Model model) {
        super(writer, model);
    }

    /**
     * Get the taxon ID from the current file name, e.g. soybean_3847.txt
     */
    public int getTaxonId() {
        try {
            String fileName = getCurrentFile().getName();
            String[] chunks = fileName.split("_");
            String[] parts = chunks[1].split("\\.");
            int taxonId = Integer.parseInt(parts[0]);
            return taxonId;
        } catch (Exception ex) {
            throw new RuntimeException("Could not parse GFF filename; format should be: no-underscores_12345.gff3 where taxonID=12345. "+ex.getMessage());
        }
    }

    /**
     * {@inheritDoc}
     * Read in the CMap file and store the linkage groups, QTLs and genetic markers along with their ranges and positions
     */
    @Override
    public void process(Reader reader) throws Exception {

        LOG.info("Processing CMap file "+getCurrentFile().getName()+"...");
        
        // these are keyed by acc in file
        Map<String,Item> linkageGroupMap = new HashMap<String,Item>();
        Map<String,Item> geneticMarkerMap = new HashMap<String,Item>();
        Map<String,Item> qtlMap = new HashMap<String,Item>();

        // these are just collected in sets for final storage
        Set<Item> linkageGroupPositionSet = new HashSet<Item>();
        Set<Item> linkageGroupRangeSet = new HashSet<Item>();

        // create the organism Item - global so can be used in methods
        int taxonId = getTaxonId();
        Item organism = createItem("Organism");
        BioStoreHook.setSOTerm(this, organism, "organism", getSequenceOntologyRefId());
        organism.setAttribute("taxonId", String.valueOf(taxonId));

        // ---------------------------------------------------------------------------------------------------------------------------
        // Load linkage groups, genetic markers and QTLs from the CMap file and associate QTLs and genetic markers with linkage groups
        // ---------------------------------------------------------------------------------------------------------------------------

        BufferedReader cmapReader = new BufferedReader(reader);
        String cmapLine = cmapReader.readLine(); // first line is header
        while ((cmapLine=cmapReader.readLine())!=null) {
            
            CMapRecord cmap = new CMapRecord(cmapLine);
            if (cmap.map_acc!=null) {
                
                // add this linkage group to its map if not already in
                Item linkageGroup = null;
                if (linkageGroupMap.containsKey(cmap.map_acc)) {
                    linkageGroup = linkageGroupMap.get(cmap.map_acc);
                } else {
                    linkageGroup = createItem("LinkageGroup");
                    BioStoreHook.setSOTerm(this, linkageGroup, "linkage_group", getSequenceOntologyRefId());
                    linkageGroup.setReference("organism", organism);
                    linkageGroup.setAttribute("secondaryIdentifier", cmap.map_acc);
                    linkageGroup.setAttribute("primaryIdentifier", cmap.map_name);
                    linkageGroup.setAttribute("length", String.valueOf(cmap.map_stop));
                    linkageGroupMap.put(cmap.map_acc, linkageGroup);
                }
                
                // add this QTL to this linkage group if appropriate
                if (cmap.isQTL()) {
                    if (!qtlMap.containsKey(cmap.feature_acc)) {
                        Item qtl = createItem("QTL");
                        BioStoreHook.setSOTerm(this, qtl, "QTL", getSequenceOntologyRefId());
                        qtl.setReference("organism", organism);
                        qtl.setAttribute("secondaryIdentifier", cmap.feature_acc);
                        qtl.setAttribute("primaryIdentifier", cmap.feature_name);
                        // create and store linkage group range; place it in map as well for future processing
                        Item linkageGroupRange = createItem("LinkageGroupRange");
                        linkageGroupRange.setAttribute("begin", String.valueOf(cmap.feature_start));
                        linkageGroupRange.setAttribute("end", String.valueOf(cmap.feature_stop));
                        linkageGroupRange.setAttribute("length", String.valueOf(round(cmap.feature_stop-cmap.feature_start,2)));
                        linkageGroupRange.setReference("linkageGroup", linkageGroup);
                        linkageGroupRangeSet.add(linkageGroupRange);
                        // add to QTL collection
                        qtl.addToCollection("linkageGroupRanges", linkageGroupRange);
                        // add to linkage group collection
                        linkageGroup.addToCollection("QTLs", qtl);
                        qtlMap.put(cmap.feature_acc, qtl);
                    }
                }
                
                // add this genetic marker to this linkage group if appropriate
                if (cmap.isMarker()) {
                    if (!geneticMarkerMap.containsKey(cmap.feature_acc)) {
                        Item geneticMarker = createItem("GeneticMarker");
                        BioStoreHook.setSOTerm(this, geneticMarker, "genetic_marker", getSequenceOntologyRefId());
                        geneticMarker.setReference("organism", organism);
                        geneticMarker.setAttribute("secondaryIdentifier", cmap.feature_acc);
                        geneticMarker.setAttribute("primaryIdentifier", cmap.feature_name);
                        geneticMarker.setAttribute("type", cmap.feature_type_acc);
                        // create and store linkage group position; place it in a map as well for future processing
                        Item linkageGroupPosition = createItem("LinkageGroupPosition");
                        linkageGroupPosition.setAttribute("position", String.valueOf(cmap.feature_start));
                        linkageGroupPosition.setReference("linkageGroup", linkageGroup);
                        linkageGroupPositionSet.add(linkageGroupPosition);
                        geneticMarker.addToCollection("linkageGroupPositions", linkageGroupPosition);
                        // add to linkage group collection
                        linkageGroup.addToCollection("geneticMarkers", geneticMarker);
                        geneticMarkerMap.put(cmap.feature_acc, geneticMarker);
                    }
                }
                
            }
            
        }
        
        cmapReader.close();
        LOG.info("Created "+linkageGroupMap.size()+" LinkageGroup items.");
        LOG.info("Created "+qtlMap.size()+" QTL items.");
        LOG.info("Created "+geneticMarkerMap.size()+" GeneticMarker items.");


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
     * Round a double to the given number of places
     */
    public static double round(double value, int places) {
        if (places < 0) throw new IllegalArgumentException();
        BigDecimal bd = new BigDecimal(value);
        bd = bd.setScale(places, RoundingMode.HALF_UP);
        return bd.doubleValue();
    }

}
