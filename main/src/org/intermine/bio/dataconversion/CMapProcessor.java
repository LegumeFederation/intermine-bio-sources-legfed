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
import java.io.FileNotFoundException;
import java.io.IOException;

import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;
import java.util.Set;

import org.apache.log4j.Logger;

import org.intermine.bio.util.OrganismData;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Attribute;
import org.intermine.xml.full.Item;
import org.intermine.xml.full.Reference;
import org.intermine.xml.full.ReferenceList;

/**
 * Store the genetic marker and QTL data from a CMap tab-delimited file.
 * NOTE: files are currently HARDCODED
 *
 * @author Sam Hokin, NCGR
 */
public class CMapProcessor extends ChadoProcessor {
	
    private static final Logger LOG = Logger.getLogger(CMapProcessor.class);

    // set true for extra verbose debug lines
    private static final boolean DEBUG = true;
    
    // count specific debug lines
    private HashMap<String,Integer> debugLineMap = new HashMap<String,Integer>();

    // global scope for use in methods
    Item organism;

    // NOTE: HARDCODED GFF FILE
    String gffFilename = "/home/intermine/data/soybean_map_version_4_0.gff3";

    // NOTE: HARDCODED CMAP FILE
    String cMapFilename = "/home/intermine/data/Soybean-GmComposite2003.txt";
        
    // global for use in methods
    Map<String,Item> chromosomeMap = new HashMap<String,Item>();
    
    // global storage maps for relating QTLs to genetic markers, populated in methods
    Map<String,Item> linkageGroupRangeMap = new HashMap<String,Item>();
    Map<String,Item> linkageGroupPositionMap = new HashMap<String,Item>();
    
    /**
     * Create a new CMapProcessor
     * @param chadoDBConverter the ChadoDBConverter that is controlling this processor
     */
    public CMapProcessor(ChadoDBConverter chadoDBConverter) {
        super(chadoDBConverter);
    }

    /**
     * {@inheritDoc}
     * We process the chado database by reading the feature, featureloc, featurepos, featuremap, feature_relationship and featureprop tables.
     */
    @Override
    public void process(Connection connection) throws SQLException, ObjectStoreException {

        // initialize our DB statement and stuff for chado queries
        Statement stmt = connection.createStatement();
        String query;
        ResultSet rs;
        
        // ---------------------------------------------------------
        // ---------------- INITIAL DATA LOADING -------------------
        // ---------------------------------------------------------

        // get chado organism_id from supplied taxon ID - enforce processing a single organism!
        Map<Integer,OrganismData> chadoToOrgData = getChadoDBConverter().getChadoIdToOrgDataMap();
        if (chadoToOrgData.size()>1) {
            System.err.println("ERROR - multiple chado organisms specified in data source; CMapProcessor can only process one organism at a time.");
            System.exit(1);
        }
        Integer organismId = 0;
        for (Integer key : chadoToOrgData.keySet()) {
            organismId = key.intValue();
        }
        OrganismData org = chadoToOrgData.get(organismId);
        int organism_id = organismId.intValue();
        int taxonId = org.getTaxonId();

        // create organism Item - global so can be used in populate routines
        organism = getChadoDBConverter().createItem("Organism");
        organism.setAttribute("taxonId", String.valueOf(taxonId));
        store(organism);
        LOG.info("Created and stored organism Item for taxonId="+taxonId+".");

        // ----------------------------------------------------------------------------------------------
        // Load the GFF data into a map for easy retrieval when item found
        // Also create chromosome Items from GFF, but updating primaryIdentifier to match chado values
        // ----------------------------------------------------------------------------------------------

        Map<String,GFFRecord> gffMap = new HashMap<String,GFFRecord>();

        // wrapped in try block due to I/O exceptions
        try {
            
            BufferedReader gffReader = new BufferedReader(new FileReader(gffFilename));
            LOG.info("Reading GFF file:"+gffFilename);
            String gffLine = null;
            while ((gffLine=gffReader.readLine()) != null) {
                GFFRecord gff = new GFFRecord(gffLine);
                if (gff.seqid!=null && gff.attributeName!=null) {
                    gffMap.put(gff.attributeName, gff);
                    if (!chromosomeMap.containsKey(gff.seqid)) {
                        Item chromosome = getChadoDBConverter().createItem("Chromosome");
                        String uniquename = gff.seqid.replace("Gm","Chr"); // chado uniquename
                        chromosome.setReference("organism", organism);
                        chromosome.setAttribute("primaryIdentifier", uniquename); // so we match up in SequenceProcessor
                        chromosomeMap.put(gff.seqid, chromosome); // use gff.seqid since we'll get chromosome from gff record
                    }
                }
            }
            gffReader.close();
            LOG.info("Read "+gffMap.size()+" GFF records from "+gffFilename+".");
            LOG.info("Created "+chromosomeMap.size()+" Chromosome items from "+gffFilename+".");

        } catch (Exception ex) {

            LOG.error(ex.getMessage());

        }


        // -----------------------------------------------------------------------------------------------------------
        // Load genes and gene families from the chado database, putting them in maps using the Item identifier as key
        // -----------------------------------------------------------------------------------------------------------

        // GENES
        // load the relevant genes from the feature table
        // key is feature_id
        Map<Integer,Item> geneMap = new HashMap<Integer,Item>();
        int geneCVTermId = 0;
        rs = stmt.executeQuery("SELECT * FROM cvterm WHERE name='gene'");
        if (rs.next()) geneCVTermId = rs.getInt("cvterm_id");
        rs.close();
        query = "SELECT * FROM feature WHERE organism_id="+organism_id+" AND type_id="+geneCVTermId;
        LOG.info("executing query: "+query);
        rs = stmt.executeQuery(query);
        while (rs.next()) {
            int feature_id = rs.getInt("feature_id");
            ChadoFeature cf = new ChadoFeature(rs);
            Item gene = getChadoDBConverter().createItem("Gene");
            cf.populateBioEntity(gene, organism);
            geneMap.put(new Integer(feature_id), gene);
        }
        rs.close();
        LOG.info("Created "+geneMap.size()+" Gene items.");

        // GENE FAMILIES
        // load the gene families from the featureprop table for our organism - distinct values with gene family type_id
        // key is value
        Map<String,Item> geneFamilyMap = new HashMap<String,Item>();
        int geneFamilyCVTermId = 0;
        rs = stmt.executeQuery("SELECT * FROM cvterm WHERE name='gene family'");
        if (rs.next()) geneFamilyCVTermId = rs.getInt("cvterm_id");
        rs.close();
        query = "SELECT DISTINCT featureprop.value FROM featureprop,feature " +
            "WHERE featureprop.type_id="+geneFamilyCVTermId+" AND featureprop.feature_id=feature.feature_id AND feature.organism_id="+organism_id;
        LOG.info("executing query: "+query);
        rs = stmt.executeQuery(query);
        while (rs.next()) {
            String value = rs.getString("value");
            Item geneFamily = getChadoDBConverter().createItem("GeneFamily");
            geneFamily.setReference("organism", organism);
            geneFamily.setAttribute("primaryIdentifier", value);
            geneFamilyMap.put(value, geneFamily);
        }
        rs.close();
        LOG.info("Created "+geneFamilyMap.size()+" GeneFamily items.");

        // ----------------------------------------------------------------------------------------------
        // Load items from the CMap file and create associations
        // ----------------------------------------------------------------------------------------------

        Map<String,Item> linkageGroupMap = new HashMap<String,Item>();
        Map<String,Item> geneticMarkerMap = new HashMap<String,Item>();
        Map<String,Item> qtlMap = new HashMap<String,Item>();

        // wrap in try block due to I/O exceptions
        try {
            
            BufferedReader cmapReader = new BufferedReader(new FileReader(cMapFilename));
            String cmapLine = cmapReader.readLine(); // header
            LOG.info("Reading CMap file:"+cMapFilename+" with header:"+cmapLine);
            
            while ((cmapLine=cmapReader.readLine())!=null) {
                
                CMapRecord cmap = new CMapRecord(cmapLine);
                if (cmap.map_acc!=null) {
                
                    // add this linkage group if not already in, checking for primaryIdentifier==cmap.map_acc
                    Item linkageGroup = getItemByPrimaryIdentifier(linkageGroupMap, cmap.map_acc);
                    if (linkageGroup==null) {
                        linkageGroup = getChadoDBConverter().createItem("LinkageGroup");
                        populateLinkageGroup(linkageGroup, cmap);
                        linkageGroupMap.put(linkageGroup.getIdentifier(), linkageGroup);
                    }
                    
                    // add a QTL if appropriate, using Item identifier as key
                    if (cmap.isQTL()) {
                        Item qtl = getChadoDBConverter().createItem("QTL");
                        populateQTL(qtl, linkageGroup, cmap);
                        qtlMap.put(qtl.getIdentifier(), qtl);
                    }
                    
                    // add a genetic marker if appropriate, using Item identifier as key
                    if (cmap.isMarker()) {
                        GFFRecord gff = gffMap.get(cmap.feature_name);
                        Item geneticMarker = getChadoDBConverter().createItem("GeneticMarker");
                        populateGeneticMarker(geneticMarker, linkageGroup, cmap, gff);
                        geneticMarkerMap.put(geneticMarker.getIdentifier(), geneticMarker);
                    }

                }
                
            }
            cmapReader.close();

        } catch (Exception ex) {

            LOG.error(ex.getMessage());

        }

        LOG.info("Created "+linkageGroupMap.size()+" LinkageGroup items.");
        LOG.info("Created "+qtlMap.size()+" QTL items.");
        LOG.info("Created "+geneticMarkerMap.size()+" GeneticMarker items.");


        // ----------------------------------------------------------------------------------------------
        // Loop through the linkage groups and associate QTLs with genetic markers within their range
        // ----------------------------------------------------------------------------------------------

        LOG.info("Associating genetic markers with QTLs...");
        for (Map.Entry<String,Item> entry : linkageGroupMap.entrySet()) {
            String linkageGroupRefId = entry.getKey();
            Item linkageGroup = entry.getValue();
            if (DEBUG) LOG.info("Linkage group:"+linkageGroup.getAttribute("primaryIdentifier").getValue());
            // QTL collection
            ReferenceList qtlRefList = linkageGroup.getCollection("QTLs");
            // GeneticMarker collection
            ReferenceList geneticMarkerRefList = linkageGroup.getCollection("geneticMarkers");
            // spin through the QTLs and add the genetic markers which fall within their range on the linkage group
            if (qtlRefList!=null && geneticMarkerRefList!=null) {
                List<String> qtlRefIdList = qtlRefList.getRefIds();
                List<String> geneticMarkerRefIdList = geneticMarkerRefList.getRefIds();
                for (String qtlRefId : qtlRefIdList) {
                    Item qtl = qtlMap.get(qtlRefId);
                    double begin = 0.0;
                    double end = 0.0;
                    // there may be multiple linkage groups for this QTL so we have to scan them
                    List<String> linkageGroupRangeRefIdList = qtl.getCollection("linkageGroupRanges").getRefIds();
                    for (String linkageGroupRangeRefId : linkageGroupRangeRefIdList) {
                        Item linkageGroupRange = linkageGroupRangeMap.get(linkageGroupRangeRefId);
                        String thisLinkageGroupRefId = linkageGroupRange.getReference("linkageGroup").getRefId();
                        if (thisLinkageGroupRefId.equals(linkageGroupRefId)) {
                            begin = Double.parseDouble(linkageGroupRange.getAttribute("begin").getValue());
                            end = Double.parseDouble(linkageGroupRange.getAttribute("end").getValue());
                            if (DEBUG) LOG.info("QTL:"+qtl.getAttribute("primaryIdentifier").getValue()+" begin="+begin+" end="+end);
                        }
                    }
                    if (end>0.0) {
                        // spin through the genetic markers and add them to qtl's collection if within range
                        // also keep track of lowest and highest genomic locations of markers for overlapping gene retrieval
                        int startMin = 100000000;
                        int endMax = 0;
                        Item chromosome = null;
                        for (String geneticMarkerRefId : geneticMarkerRefIdList) {
                            Item geneticMarker = geneticMarkerMap.get(geneticMarkerRefId);
                            List<String> linkageGroupPositionRefIdList = geneticMarker.getCollection("linkageGroupPositions").getRefIds();
                            for (String linkageGroupPositionRefId : linkageGroupPositionRefIdList) {
                                Item linkageGroupPosition = linkageGroupPositionMap.get(linkageGroupPositionRefId);
                                String thisLinkageGroupRefId = linkageGroupPosition.getReference("linkageGroup").getRefId();
                                if (thisLinkageGroupRefId.equals(linkageGroupRefId)) {
                                    // this position is on this linkage group, now check for QTL overlap
                                    double position = Double.parseDouble(linkageGroupPosition.getAttribute("position").getValue());
                                    if (position>=begin && position<=end) {
                                        if (DEBUG) LOG.info("..."+geneticMarker.getAttribute("primaryIdentifier").getValue()+": position="+position);
                                        // marker is in QTL span, add it QTL's geneticMarkers collection
                                        qtl.addToCollection("geneticMarkers", geneticMarker);
                                        // update startMin and endMax from this marker IF present in the GFF file
                                        String feature_name = geneticMarker.getAttribute("secondaryIdentifier").getValue();
                                        GFFRecord gff = gffMap.get(feature_name);
                                        if (gff!=null) {
                                            if (gff.start<startMin) startMin = gff.start;
                                            if (gff.end>endMax) endMax = gff.end;
                                            chromosome = chromosomeMap.get(gff.seqid); // chromosome should be same for all markers here
                                        }
                                    }
                                }
                            }
                        }
                        if (endMax>0) {
                            // get the chado feature_id of the chromosome
                            String uniquename = chromosome.getAttribute("primaryIdentifier").getValue();
                            query = "SELECT feature_id FROM feature WHERE organism_id="+organism_id+" AND uniquename='"+uniquename+"'";
                            rs = stmt.executeQuery(query);
                            int srcfeature_id = 0;
                            if (rs.next()) srcfeature_id = rs.getInt("feature_id");
                            rs.close();
                            if (srcfeature_id>0) {
                                // find the genes that overlap [startMin,endMax] on this chromosome
                                query = "SELECT featureloc.feature_id FROM featureloc,feature WHERE " +
                                    "featureloc.feature_id=feature.feature_id AND type_id="+geneCVTermId+" AND srcfeature_id="+srcfeature_id+" AND fmin<"+endMax+" AND fmax>"+startMin;
                                if (DEBUG) LOG.info("executing query: "+query);
                                rs = stmt.executeQuery(query);
                                while (rs.next()) {
                                    int feature_id = rs.getInt("feature_id");
                                    Item gene = geneMap.get(new Integer(feature_id));
                                    if (gene!=null) {
                                        qtl.addToCollection("overlappingGenes", gene);
                                    }
                                }
                                rs.close();
                            }
                        }
                    } // end>0
                }
            }
        }
        LOG.info("...done.");

        
        // ----------------------------------------------------------------------------------------------
        // Chado data associations
        // ----------------------------------------------------------------------------------------------

        // query the featureprop table to associate our genes with gene families
        LOG.info("Associating genes with gene families from featureprop....");
        for (Map.Entry<Integer,Item> entry : geneMap.entrySet()) {
            int feature_id = (int) entry.getKey();
            Item gene = entry.getValue();
            query = "SELECT * FROM featureprop WHERE feature_id="+feature_id;
            rs = stmt.executeQuery(query);
            while (rs.next()) {
                int featureprop_id = rs.getInt("featureprop_id"); // not used
                int type_id = rs.getInt("type_id");               // one of three possibilities
                String value = rs.getString("value");             // geneFamilyMap key (description)
                // branch on type_id
                if (type_id==2125) {
                    gene.setAttribute("note", value);
                } else if (type_id==78562) {
                    gene.setAttribute("familyRepresentative", value);
                } else if (type_id==51206) {
                    // add reference to the GeneFamily to this Gene
                    Item geneFamily = geneFamilyMap.get(value);
                    if (geneFamily!=null) {
                        gene.setReference("geneFamily", geneFamily);
                    }
                }
            }
            rs.close();
        }
        LOG.info("...done.");


        // ----------------------------------------------------------------
        // we're done, so store everything
        // ----------------------------------------------------------------

        LOG.info("Storing chromosomes...");
        for (Item item : chromosomeMap.values()) {
            store(item);
        }
        
        LOG.info("Storing linkage groups...");
        for (Item item : linkageGroupMap.values()) {
            store(item);
        }

        LOG.info("Storing genetic markers...");
        for (Item item : geneticMarkerMap.values()) {
            store(item);
        }
 
        LOG.info("Storing QTLs...");
        for (Item item : qtlMap.values()) {
            store(item);
        }

        LOG.info("Storing genes...");
        for (Item item : geneMap.values()) {
            store(item);
        }

        LOG.info("Storing gene families...");
        for (Item item : geneFamilyMap.values()) {
            store(item);
        }


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
     * Store and tally debug lines in debugLineMap.
     * @param a string to be output
     */
    private void debugToLog(String output) {
	if (debugLineMap.containsKey(output)) {
	    int count = debugLineMap.get(output).intValue();
	    count++;
	    debugLineMap.put(output, new Integer(count));
	} else {
	    debugLineMap.put(output, new Integer(1));
	}
    }

    /**
     * Dump the debugLineMap and clear it.
     */
    private void logDebugLines() {
	for (String debugLine : debugLineMap.keySet()) {
	    LOG.info(debugLine+" {"+debugLineMap.get(debugLine).intValue()+"}");
	}
	debugLineMap.clear();
    }

    /**
     * Populate a LinkageGroup Item from a CMap record
     */
    void populateLinkageGroup(Item linkageGroup, CMapRecord cmap) {
        linkageGroup.setReference("organism", organism);
        linkageGroup.setAttribute("primaryIdentifier", cmap.map_acc);
        if (cmap.map_name!=null && cmap.map_name.length()>0) {
            linkageGroup.setAttribute("secondaryIdentifier", cmap.map_name);
        }
        linkageGroup.setAttribute("length", String.valueOf(cmap.map_stop));
    }

    /**
     * Populate a QTL Item from a related LinkageGroup Item and CMap record
     */
    void populateQTL(Item qtl, Item linkageGroup, CMapRecord cmap) throws ObjectStoreException {
        qtl.setReference("organism", organism);
        qtl.setAttribute("primaryIdentifier", cmap.feature_acc);
        if (cmap.feature_name!=null && cmap.feature_name.length()>0) {
            qtl.setAttribute("secondaryIdentifier", cmap.feature_name);
        }
        // create and store linkage group range; place it in map as well for future processing
        Item linkageGroupRange = getChadoDBConverter().createItem("LinkageGroupRange");
        linkageGroupRange.setAttribute("begin", String.valueOf(cmap.feature_start));
        linkageGroupRange.setAttribute("end", String.valueOf(cmap.feature_stop));
        linkageGroupRange.setReference("linkageGroup", linkageGroup);
        store(linkageGroupRange);
        qtl.addToCollection("linkageGroupRanges", linkageGroupRange);
        // add to map
        linkageGroupRangeMap.put(linkageGroupRange.getIdentifier(), linkageGroupRange);
        // add to linkage group collection
        linkageGroup.addToCollection("QTLs", qtl);
    }

    /**
     * Populate a GeneticMarker Item from a LinkageGroup Item, CMap record and GFF record
     */
    void populateGeneticMarker(Item geneticMarker, Item linkageGroup, CMapRecord cmap, GFFRecord gff) throws ObjectStoreException {
        geneticMarker.setReference("organism", organism);
        geneticMarker.setAttribute("primaryIdentifier", cmap.feature_acc);
        if (cmap.feature_name!=null && cmap.feature_name.length()>0) {
            geneticMarker.setAttribute("secondaryIdentifier", cmap.feature_name);
        }
        geneticMarker.setAttribute("type", cmap.feature_type_acc);
        if (gff!=null) {
            // get chromosome Item from chromosomeMap
            Item chromosome = chromosomeMap.get(gff.seqid);
            if (chromosome!=null) {
                geneticMarker.setReference("chromosome", chromosome);
                geneticMarker.setAttribute("length", String.valueOf(gff.end-gff.start+1));
                Item location = getChadoDBConverter().createItem("Location");
                location.setAttribute("start", String.valueOf(gff.start));
                location.setAttribute("end", String.valueOf(gff.end));
                location.setReference("feature", geneticMarker);
                location.setReference("locatedOn", chromosome);
                store(location);
                geneticMarker.setReference("chromosomeLocation", location);
            }
        }
        // create and store linkage group position; place it in a map as well for future processing
        Item linkageGroupPosition = getChadoDBConverter().createItem("LinkageGroupPosition");
        linkageGroupPosition.setAttribute("position", String.valueOf(cmap.feature_start));
        linkageGroupPosition.setReference("linkageGroup", linkageGroup);
        store(linkageGroupPosition);
        geneticMarker.addToCollection("linkageGroupPositions", linkageGroupPosition);
        // add to map
        linkageGroupPositionMap.put(linkageGroupPosition.getIdentifier(), linkageGroupPosition);
        // add to linkage group collection
        linkageGroup.addToCollection("geneticMarkers", geneticMarker);
    }

    /**
     * Get an Item by its primaryIdentifier; return null if not found
     */
    Item getItemByPrimaryIdentifier(Map<String,Item> map, String primaryIdentifier) {
        for (Item item : map.values()) {
            String primaryId = item.getAttribute("primaryIdentifier").getValue();
            if (primaryId.equals(primaryIdentifier)) return item;
        }
        return null;
    }

    /**
     * Get an Item by its chadoFeatureId; return null if not found
     */
    Item getItemByChadoFeatureId(Map<String,Item> map, int chadoFeatureId) {
        for (Item item : map.values()) {
            int cfId = Integer.parseInt(item.getAttribute("chadoFeatureId").getValue());
            if (cfId==chadoFeatureId) return item;
        }
        return null;
    }

}
