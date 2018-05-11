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
 * DEPRECATED - DOES NOT LOAD ORG and OTHER DATA FROM HEADER LINES; DOES NOT SUPPORT Organism.variety.
 *
 * Store the genetic marker and QTL data from a CMap tab-delimited file. Fields are:
 *
 * map_acc map_name map_start map_stop feature_acc feature_name feature_aliases feature_start feature_stop feature_type_acc is_landmark
 *
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
    Map<String,Item> markerMap = new HashMap<String,Item>();   // keyed by acc in file
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
     * {@inheritDoc}
     * Read in the CMap file and store the linkage groups, QTLs and genetic markers along with their ranges and positions
     */
    @Override
    public void process(Reader reader) throws Exception {

        // don't process README files
        if (getCurrentFile().getName().contains("README")) return;

        LOG.info("Processing CMap file "+getCurrentFile().getName()+"...");
        
        // create and add the organism Item to its map
        Item organism = createItem("Organism");
        organism.setAttribute("taxonId", getTaxonId());
        organismSet.add(organism);

        // create and store the genetic map
        Item geneticMap = createItem("GeneticMap");
        geneticMap.setAttribute("primaryIdentifier", getGeneticMapName());
        geneticMap.setReference("organism", organism);
        store(geneticMap);

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
                // change length if this record shows a larger map_stop (update to handle ADF-generated files)
                Item linkageGroup = null;
                if (linkageGroupMap.containsKey(cmap.map_acc)) {
                    linkageGroup = linkageGroupMap.get(cmap.map_acc);
                    double currentLength = Double.parseDouble(linkageGroup.getAttribute("length").getValue());
                    if (currentLength<cmap.map_stop) linkageGroup.setAttribute("length", String.valueOf(cmap.map_stop));
                } else {
                    linkageGroup = createItem("LinkageGroup");
                    linkageGroup.setAttribute("primaryIdentifier", cmap.map_acc);
                    linkageGroup.setAttribute("secondaryIdentifier", cmap.map_name);
                    linkageGroup.setAttribute("length", String.valueOf(cmap.map_stop));
                    linkageGroup.setReference("organism", organism);
                    linkageGroup.setReference("geneticMap", geneticMap);
                    linkageGroupMap.put(cmap.map_acc, linkageGroup);
                }
                
                // add this QTL to this linkage group if appropriate
                // we'll use map_name as primaryIdentifier since it's hopefully unique and concise
                if (cmap.isQTL()) {
                    if (!qtlMap.containsKey(cmap.feature_acc)) {
                        Item qtl = createItem("QTL");
                        qtl.setReference("organism", organism);
                        if (cmap.feature_name.contains(":")) {
                            // use the part before colon for primary identifier
                            String parts[] = cmap.feature_name.split(":");
                            qtl.setAttribute("primaryIdentifier", parts[0]);
                        } else {
                            qtl.setAttribute("primaryIdentifier", cmap.feature_name);
                        }
                        if (cmap.feature_acc.contains(":")) {
                            // use the part after colon for secondary identifier
                            String parts[] = cmap.feature_acc.split(":");
                            qtl.setAttribute("secondaryIdentifier", parts[1]);
                        } else {
                            qtl.setAttribute("secondaryIdentifier", cmap.feature_acc);
                        }
                        // create and store linkage group range; place it in map as well for future processing
                        Item linkageGroupRange = createItem("LinkageGroupRange");
                        linkageGroupRange.setAttribute("begin", String.valueOf(cmap.feature_start));
                        linkageGroupRange.setAttribute("end", String.valueOf(cmap.feature_stop));
                        linkageGroupRange.setAttribute("length", String.valueOf(round(cmap.feature_stop-cmap.feature_start,2)));
                        linkageGroupRange.setReference("linkageGroup", linkageGroup);
                        linkageGroupRangeSet.add(linkageGroupRange);
                        // SOYBASE: add a comment if the length is exactly 2.0 cM - artificially imposed
                        if (getTaxonId().equals("3847") && (cmap.feature_stop-cmap.feature_start)==2.0) {
                            qtl.setAttribute("description", "Length on linkage group arbitrarily set to 2.0 cM.");
                        }
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
                    if (!markerMap.containsKey(cmap.feature_acc)) {
                        Item marker = createItem("GeneticMarker");
                        marker.setReference("organism", organism);
                        marker.setAttribute("primaryIdentifier", cmap.feature_name);
                        marker.setAttribute("secondaryIdentifier", cmap.feature_acc);
                        marker.setAttribute("type", cmap.feature_type_acc);
                        // create and store linkage group position; place it in a map as well for future processing
                        Item linkageGroupPosition = createItem("LinkageGroupPosition");
                        linkageGroupPosition.setAttribute("position", String.valueOf(cmap.feature_start));
                        linkageGroupPosition.setReference("linkageGroup", linkageGroup);
                        linkageGroupPositionSet.add(linkageGroupPosition);
                        marker.addToCollection("linkageGroupPositions", linkageGroupPosition);
                        // add to linkage group collection
                        linkageGroup.addToCollection("markers", marker);
                        markerMap.put(cmap.feature_acc, marker);
                    }
                }
                
            }
            
        }
        
        cmapReader.close();
        LOG.info("Created "+linkageGroupMap.size()+" LinkageGroup items.");
        LOG.info("Created "+qtlMap.size()+" QTL items.");
        LOG.info("Created "+markerMap.size()+" GeneticMarker items.");
 
    }


    /**
     * Store the items we've collected from the CMap files
     */
    @Override
    public void close() throws Exception {
    
        LOG.info("Storing "+organismSet.size()+" organisms...");
        store(organismSet);

        LOG.info("Storing "+linkageGroupMap.size()+" linkage groups...");
        store(linkageGroupMap.values());

        LOG.info("Storing "+linkageGroupPositionSet.size()+" linkage group positions...");
        store(linkageGroupPositionSet);

        LOG.info("Storing "+linkageGroupRangeSet.size()+" linkage group ranges...");
        store(linkageGroupRangeSet);

        LOG.info("Storing "+qtlMap.size()+" QTLs...");
        store(qtlMap.values());

        LOG.info("Storing "+markerMap.size()+" genetic markers...");
        store(markerMap.values());

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

    /**
     * Get the genetic map name from the file name, e.g. GeneticMapFoo_3917_24659904.gt
     */
    public String getGeneticMapName() {
        String fileName = getCurrentFile().getName();
        String[] chunks = fileName.split("_");
        return chunks[0];
    }

    /**
     * Get the taxon ID from the current file name, e.g. GeneticMapFoo_3917_24659904.gt
     */
    public String getTaxonId() {
        String fileName = getCurrentFile().getName();
        String[] chunks = fileName.split("_");
        return chunks[1];
    }

}
