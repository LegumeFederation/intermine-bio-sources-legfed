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

    // store items at end in close() method to avoid dupes
    Set<Item> organismSet = new HashSet<Item>();

    Map<String,Item> linkageGroupMap = new HashMap<String,Item>();    // keyed by acc in file
    Map<String,Item> geneticMarkerMap = new HashMap<String,Item>();   // keyed by acc in file
    Map<String,Item> qtlMap = new HashMap<String,Item>();             // keyed by acc in file

    Set<Item> linkageGroupPositionSet = new HashSet<Item>();
    Set<Item> linkageGroupRangeSet = new HashSet<Item>();
    
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
    public String getTaxonId() {
        try {
            String fileName = getCurrentFile().getName();
            String[] chunks = fileName.split("_");
            String[] parts = chunks[1].split("\\.");
            return parts[0];
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
        
        // create and add the organism Item to its map
        Item organism = createItem("Organism");
        BioStoreHook.setSOTerm(this, organism, "organism", getSequenceOntologyRefId());
        organism.setAttribute("taxonId", getTaxonId());
        organismSet.add(organism);

        // ---------------------------------------------------------------------------------------------------------------------------
        // Load linkage groups, genetic markers and QTLs from the CMap file and associate QTLs and genetic markers with linkage groups
        // ---------------------------------------------------------------------------------------------------------------------------

        BufferedReader cmapReader = new BufferedReader(reader);
        String cmapLine = cmapReader.readLine(); // first line is header
        while ((cmapLine=cmapReader.readLine())!=null) {
            
            CMapRecord cmap = new CMapRecord(cmapLine);
            if (cmap.map_acc!=null) {
                
                // add this linkage group to its map if not already in
                // we'll use map_acc for the primaryIdentifier since map_name is often something minimal like "3" which is likely not unique
                Item linkageGroup = null;
                if (linkageGroupMap.containsKey(cmap.map_acc)) {
                    linkageGroup = linkageGroupMap.get(cmap.map_acc);
                } else {
                    linkageGroup = createItem("LinkageGroup");
                    BioStoreHook.setSOTerm(this, linkageGroup, "linkage_group", getSequenceOntologyRefId());
                    linkageGroup.setReference("organism", organism);
                    linkageGroup.setAttribute("primaryIdentifier", cmap.map_acc);
                    linkageGroup.setAttribute("secondaryIdentifier", cmap.map_name);
                    linkageGroup.setAttribute("length", String.valueOf(cmap.map_stop));
                    linkageGroupMap.put(cmap.map_acc, linkageGroup);
                }
                
                // add this QTL to this linkage group if appropriate
                // we'll use map_name as primaryIdentifier since it's hopefully unique and concise
                if (cmap.isQTL()) {
                    if (!qtlMap.containsKey(cmap.feature_acc)) {
                        Item qtl = createItem("QTL");
                        BioStoreHook.setSOTerm(this, qtl, "QTL", getSequenceOntologyRefId());
                        qtl.setReference("organism", organism);
                        qtl.setAttribute("primaryIdentifier", cmap.feature_name);
                        qtl.setAttribute("secondaryIdentifier", cmap.feature_acc);
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
                // we'll use map_name as primaryIdentifier since it's hopefully unique and concise
                if (cmap.isMarker()) {
                    if (!geneticMarkerMap.containsKey(cmap.feature_acc)) {
                        Item geneticMarker = createItem("GeneticMarker");
                        BioStoreHook.setSOTerm(this, geneticMarker, "genetic_marker", getSequenceOntologyRefId());
                        geneticMarker.setReference("organism", organism);
                        geneticMarker.setAttribute("primaryIdentifier", cmap.feature_name);
                        geneticMarker.setAttribute("secondaryIdentifier", cmap.feature_acc);
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
 
    }


    /**
     * Store the items we've collected from the CMap files
     */
    @Override
    public void close() throws Exception {
    
        LOG.info("Storing "+organismSet.size()+" organisms...");
        for (Item item : organismSet) store(item);

        LOG.info("Storing "+linkageGroupMap.size()+" linkage groups...");
        for (Item item : linkageGroupMap.values()) store(item);

        LOG.info("Storing "+linkageGroupPositionSet.size()+" linkage group positions...");
        for (Item item : linkageGroupPositionSet) store(item);

        LOG.info("Storing "+linkageGroupRangeSet.size()+" linkage group ranges...");
        for (Item item : linkageGroupRangeSet) store(item);

        LOG.info("Storing "+qtlMap.size()+" QTLs...");
        for (Item item : qtlMap.values()) store(item);

        LOG.info("Storing "+geneticMarkerMap.size()+" genetic markers...");
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
