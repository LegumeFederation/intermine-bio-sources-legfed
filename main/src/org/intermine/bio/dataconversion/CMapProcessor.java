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
 * Store the genetic marker and QTL data from a CMap tab-delimited file, along with a GFF file containing 
 * genetic markers and a QTL-markers file containing QTL-marker relations.
 *
 * NOTE: input text files are currently HARDCODED (Soybase)
 *
 * @author Sam Hokin, NCGR
 */
public class CMapProcessor extends ChadoProcessor {
	
    private static final Logger LOG = Logger.getLogger(CMapProcessor.class);

    // global scope for use in methods
    Item organism;

    // NOTE: HARDCODED GFF FILE
    String gffFilename = "/home/intermine/data/soybean_map_version_4_0.gff3";

    // NOTE: HARDCODED CMAP FILE
    String cMapFilename = "/home/intermine/data/Soybean-GmComposite2003.txt";

    // NOTE: HARDCODED QTL-MARKERS FILE
    String qtlMarkersFilename = "/home/intermine/data/SoyBase-QTL-markers.txt";
        
    // global for use in methods
    ItemMap chromosomeMap = new ItemMap();
    
    // global storage maps for relating QTLs to genetic markers, populated in methods
    ItemMap linkageGroupRangeMap = new ItemMap();
    ItemMap linkageGroupPositionMap = new ItemMap();
    
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
        int organism_id = organismId.intValue(); // chado organism.organism_id
        int taxonId = org.getTaxonId();

        // create organism Item - global so can be used in populate routines
        organism = getChadoDBConverter().createItem("Organism");
        organism.setAttribute("taxonId", String.valueOf(taxonId));
        store(organism);
        LOG.info("Created and stored organism Item for taxonId="+taxonId+".");

        // -------------------------------------------------------------------------------------------------------------------
        // Load the GFF data into a map for easy retrieval when item found. Also create chromosome Items from GFF.
        //
        // NOTE: this is hardcoded for Soybase chromosomes: Soybase GFF and LIS chado don't use same chromosome IDs,
        //       so we update primaryIdentifier to match chado values.
        // -------------------------------------------------------------------------------------------------------------------

        Map<String,GFFRecord> gffMap = new HashMap<String,GFFRecord>();

        try {
            
            BufferedReader gffReader = new BufferedReader(new FileReader(gffFilename));
            LOG.info("Reading GFF file:"+gffFilename);
            String gffLine = null;
            while ((gffLine=gffReader.readLine()) != null) {
                GFFRecord gff = new GFFRecord(gffLine);
                if (gff.seqid!=null && gff.attributeName!=null) {
                    // NOTE: key is attributeName
                    gffMap.put(gff.attributeName, gff);
                    // create chromosome and add to chromosomeMap if needed; avoid scaffolds
                    String name = gff.seqid.replace("Gm","glyma.Chr"); // chado chromosome name HARDCODED
                    if (!chromosomeMap.containsSecondaryIdentifier(name)) {
                        // create the chromosome Item from the chado record
                        query = "SELECT * FROM feature WHERE type_id=43403 AND organism_id="+organism_id+" AND name='"+name+"'";
                        rs = stmt.executeQuery(query);
                        if (rs.next()) {
                            ChadoFeature cf = new ChadoFeature(rs);
                            Item chromosome = getChadoDBConverter().createItem("Chromosome");
                            cf.populateBioEntity(chromosome, organism);
                            chromosomeMap.put(chromosome);
                        }
                        rs.close();
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
            // NOTE: feature_id is key
            geneMap.put(new Integer(feature_id), gene);
        }
        rs.close();
        LOG.info("Created "+geneMap.size()+" Gene items.");

        // GENE FAMILIES
        // load the gene families from the featureprop table for our organism - distinct values with gene family type_id
        ItemMap geneFamilyMap = new ItemMap();
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
            geneFamily.setAttribute("primaryIdentifier", value);
            // NOTE: value is key
            geneFamilyMap.put(value, geneFamily);
        }
        rs.close();
        LOG.info("Created "+geneFamilyMap.size()+" GeneFamily items.");

        // ----------------------------------------------------------------------------------------------
        // Load items from the CMap file and associate with linkage groups
        // ----------------------------------------------------------------------------------------------

        ItemMap linkageGroupMap = new ItemMap();
        ItemMap geneticMarkerMap = new ItemMap();
        ItemMap qtlMap = new ItemMap();

        try {
            
            BufferedReader cmapReader = new BufferedReader(new FileReader(cMapFilename));
            String cmapLine = cmapReader.readLine(); // header
            LOG.info("Reading CMap file:"+cMapFilename+" with header:"+cmapLine);
            
            while ((cmapLine=cmapReader.readLine())!=null) {
                
                CMapRecord cmap = new CMapRecord(cmapLine);
                if (cmap.map_acc!=null) {
                
                    // add this linkage group if not already in, checking for primaryIdentifier==cmap.map_acc
                    Item linkageGroup = linkageGroupMap.getByPrimaryIdentifier(cmap.map_acc);
                    if (linkageGroup==null) {
                        linkageGroup = getChadoDBConverter().createItem("LinkageGroup");
                        populateLinkageGroup(linkageGroup, cmap);
                        linkageGroupMap.put(linkageGroup);
                    }
                    
                    // add a QTL if appropriate, using Item identifier as key
                    if (cmap.isQTL()) {
                        Item qtl = getChadoDBConverter().createItem("QTL");
                        populateQTL(qtl, linkageGroup, cmap);
                        qtlMap.put(qtl);
                    }
                    
                    // add a genetic marker if appropriate, using Item identifier as key
                    if (cmap.isMarker()) {
                        GFFRecord gff = gffMap.get(cmap.feature_name);
                        Item geneticMarker = getChadoDBConverter().createItem("GeneticMarker");
                        populateGeneticMarker(geneticMarker, linkageGroup, cmap, gff);
                        geneticMarkerMap.put(geneticMarker);
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

        // ----------------------------------------------------------------------------------------------
        // Run through the QTL-Markers file and add the associated markers to the given QTLs
        // NOTE: marker ZZ is a placeholder, not a real marker
        // ----------------------------------------------------------------------------------------------
        
        try {
            
            BufferedReader qtlMarkersReader = new BufferedReader(new FileReader(qtlMarkersFilename));
            String qtlMarkersLine = qtlMarkersReader.readLine(); // header
            LOG.info("Reading QTL-Markers file:"+qtlMarkersFilename+" with header:"+qtlMarkersLine);
            
            while ((qtlMarkersLine=qtlMarkersReader.readLine())!=null) {
                QTLMarkersRecord rec = new QTLMarkersRecord(qtlMarkersLine);
                if (rec.qtlName!=null && rec.markerName!=null && !rec.markerName.equals("ZZ")) {
                    // find the QTL given that qtlName=secondaryIdentifier
                    Item qtl = qtlMap.getBySecondaryIdentifier(rec.qtlName);
                    if (qtl!=null) {
                        // now find the genetic marker given that markerName=secondaryIdentifier
                        Item geneticMarker = geneticMarkerMap.getBySecondaryIdentifier(rec.markerName);
                        if (geneticMarker!=null) {
                            // add this genetic marker to this QTL's collection
                            qtl.addToCollection("associatedGeneticMarkers", geneticMarker);
                        } else {
                            LOG.info("Genetic marker "+rec.markerName+" is missing from geneticMarkerMap.");
                        }
                    } else {
                        LOG.info("QTL "+rec.qtlName+" is missing from qtlMap.");
                    }
                }
            }
            
            qtlMarkersReader.close();
            LOG.info("Done associating genetic markers with QTLs.");

        } catch (Exception ex) {

            LOG.error(ex.getMessage());

        }

        // ----------------------------------------------------------------------------------------------------------
        // Loop through the QTLs and find the genes which overlap the region spanned by the QTL's associated markers.
        // NOTE: only include markers which are on the same linkage group as the QTL. Sometimes a QTL gets
        // extra markers from other linkage groups. ASSUMPTION: QTL has only one linkage group. True for Soybase.
        // ---------------------------------------------------------------------------------------------------------

        LOG.info("Associating genes with QTLs via associated marker(s)...");
        for (Item qtl : qtlMap.values()) {

            List<Item> geneticMarkers = geneticMarkerMap.getForCollection(qtl, "associatedGeneticMarkers");
            List<Item> linkageGroupRanges = linkageGroupRangeMap.getForCollection(qtl, "linkageGroupRanges");

            if (geneticMarkers.size()>0 && linkageGroupRanges.size()>0) {

                // get this QTL's linkage group for marker check below, assumption is there's only one
                Item qtlLinkageGroup = null;
                for (Item linkageGroupRange : linkageGroupRanges) {
                    qtlLinkageGroup = linkageGroupMap.getForReference(linkageGroupRange, "linkageGroup");
                }
                // spin through this QTL's markers and find the full genomic range spanned: chromosome[minStart,maxEnd]
                int minStart = 100000000;
                int maxEnd = 0;
                int srcfeature_id = 0;
                for (Item geneticMarker : geneticMarkers) {
                    // check that we're on the same linkage group as QTL, again assume single linkage group
                    Item gmLinkageGroup= null;
                    List<Item> linkageGroupPositions = linkageGroupPositionMap.getForCollection(geneticMarker, "linkageGroupPositions");
                    for (Item linkageGroupPosition : linkageGroupPositions) {
                        gmLinkageGroup = linkageGroupMap.getForReference(linkageGroupPosition, "linkageGroup");
                    }
                    if (gmLinkageGroup!=null && gmLinkageGroup.equals(qtlLinkageGroup)) {
                        Item chromosome = chromosomeMap.getForReference(geneticMarker, "chromosome");
                        if (chromosome!=null) {
                            srcfeature_id = Integer.parseInt(chromosome.getAttribute("chadoFeatureId").getValue());
                            // get this marker's range on chromosome from GFF record
                            String attributeName = geneticMarker.getAttribute("secondaryIdentifier").getValue();
                            GFFRecord gff = gffMap.get(attributeName);
                            if (gff!=null) {
                                if (gff.start<minStart) minStart = gff.start;
                                if (gff.end>maxEnd) maxEnd = gff.end;
                            }
                        }
                    }
                }

                // now query the genes within chromosome[minStart,maxEnd] and associate them with this QTL
                if (srcfeature_id>0 && maxEnd>0) {
                    query = "SELECT featureloc.feature_id FROM featureloc,feature WHERE " +
                        "featureloc.feature_id=feature.feature_id AND type_id="+geneCVTermId+" AND srcfeature_id="+srcfeature_id+" AND fmin<="+maxEnd+" AND fmax>="+minStart;
                    rs = stmt.executeQuery(query);
                    while (rs.next()) {
                        int feature_id = rs.getInt("feature_id");
                        Item gene = geneMap.get(new Integer(feature_id));
                        if (gene!=null) {
                            qtl.addToCollection("overlappedGenes", gene);
                        }
                    }
                    rs.close();
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
        linkageGroupRangeMap.put(linkageGroupRange);
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
            String uniquename = gff.seqid.replace("Gm","Chr"); // chado uniquename
            Item chromosome = chromosomeMap.getByPrimaryIdentifier(uniquename);
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
        linkageGroupPositionMap.put(linkageGroupPosition);
        // add to linkage group collection
        linkageGroup.addToCollection("geneticMarkers", geneticMarker);
    }

}
