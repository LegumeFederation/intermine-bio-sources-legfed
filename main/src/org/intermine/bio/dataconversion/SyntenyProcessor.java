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
 * Read synteny blocks for two organisms from a DAGchainer synteny GFF file and store them as SyntenyBlock items.
 * This is designed to use the GFF annotation produced by DAGchainer.
 *
 * NOTE: input text files are currently HARDCODED
 *
 * @author Sam Hokin, NCGR
 */
public class SyntenyProcessor extends ChadoProcessor {
	
    private static final Logger LOG = Logger.getLogger(SyntenyProcessor.class);

    // DAGchainer synteny gff file location is passed in as parameter dagChainerFile in project.xml
    File dagChainerFile;

    // InterMine organism IDs
    int sourceTaxonId;
    int targetTaxonId;

    // chado organism_ids
    int source_organism_id;
    int target_organism_id;

    /**
     * Create a new SyntenyProcessor
     * @param chadoDBConverter the ChadoDBConverter that is controlling this processor
     */
    public SyntenyProcessor(ChadoDBConverter chadoDBConverter) {
        super(chadoDBConverter);
    }

    /**
     * {@inheritDoc}
     * We process the GFF file by creating SytenicRegion items and storing them.
     * In order to do so, we need to access the two organisms' chromosomes from the chado database.
     */
    @Override
    public void process(Connection connection) throws SQLException, ObjectStoreException {

        // ---------------------------------------------------------
        // ---------------- INITIAL DATA LOADING -------------------
        // ---------------------------------------------------------

        // initialize our DB statement and stuff for chado queries
        Statement stmt = connection.createStatement();
        String query;
        ResultSet rs;

        // get the taxonomy IDs corresponding to the source and target
        try {
            sourceTaxonId = getChadoDBConverter().getSourceOrganism();
            targetTaxonId = getChadoDBConverter().getTargetOrganism();
        } catch (Exception ex) {
            throw new RuntimeException("Error retrieving source and target taxonomy IDs from project.xml: "+ex.toString());
        }
        LOG.info("sourceTaxonId="+sourceTaxonId);
        LOG.info("targetTaxonId="+targetTaxonId);
        if (sourceTaxonId==0 || targetTaxonId==0) {
            LOG.error("Error retrieving source and target taxonomy IDs from project.xml.");
            // HACK - until I figure out the bug in this
            sourceTaxonId = 3885;
            targetTaxonId = 3847;
        }

        // get the chado organism_ids, which are keys of the ChadoIdToOrgDataMap;
        try {
            for (Map.Entry<Integer,OrganismData> entry : getChadoDBConverter().getChadoIdToOrgDataMap().entrySet()) {
                if (entry.getValue().getTaxonId()==sourceTaxonId) source_organism_id = (int)entry.getKey();
                if (entry.getValue().getTaxonId()==targetTaxonId) target_organism_id = (int)entry.getKey();
            }
        } catch (Exception ex) {
            throw new RuntimeException("Error determining source and target chado organism_ids: "+ex.toString());
        }
        LOG.info("source_organism_id="+source_organism_id);
        LOG.info("target_organism_id="+target_organism_id);
        if (source_organism_id==0 || target_organism_id==0) {
            LOG.error("Error retrieving source and target chado organism_ids from project.xml.");
            throw new RuntimeException("Error retrieving source and target chado organism_ids from project.xml.");
        }

        // create and store organism Items
        Item sourceOrganism = getChadoDBConverter().createItem("Organism");
        sourceOrganism.setAttribute("taxonId", String.valueOf(sourceTaxonId));
        store(sourceOrganism);
        Item targetOrganism = getChadoDBConverter().createItem("Organism");
        targetOrganism.setAttribute("taxonId", String.valueOf(targetTaxonId));
        store(targetOrganism);
        LOG.info("Created and stored sourceOrganism, taxonId="+sourceTaxonId+" and targetOrganism, taxonId="+targetTaxonId);

        // get the all-important DAGchainer GFF file
        dagChainerFile = getChadoDBConverter().getDagChainerFile();
        LOG.info("dagChainerFile="+dagChainerFile.getAbsolutePath());
        
        // -------------------------------------------------------------------------------------------------------------------
        // Load the GFF data into a map. Create source and target chromosome Items.
        //
        // NOTE: this is hardcoded for Pv source and Gm target chromosomes.
        // -------------------------------------------------------------------------------------------------------------------

        Map<String,GFFRecord> gffMap = new HashMap<String,GFFRecord>();
        ItemMap sourceChromosomeMap = new ItemMap();
        ItemMap targetChromosomeMap = new ItemMap();
        
        try {
            
            BufferedReader gffReader = new BufferedReader(new FileReader(dagChainerFile));
            LOG.info("Reading GFF file:"+dagChainerFile.getName());
            String gffLine = null;
            while ((gffLine=gffReader.readLine()) != null) {
                GFFRecord gff = new GFFRecord(gffLine);
                if (gff.seqid!=null && gff.attributeName!=null) {
                    // create source chromosome, store and add to sourceChromosomeMap if needed
                    String sourceName = gff.seqid.replace("Pv","phavu.Chr"); // chado chromosome name HARDCODED
                    if (!sourceChromosomeMap.containsSecondaryIdentifier(sourceName)) {
                        // create the chromosome Item from the chado record
                        query = "SELECT * FROM feature WHERE type_id=43403 AND organism_id="+source_organism_id+" AND name='"+sourceName+"'";
                        rs = stmt.executeQuery(query);
                        if (rs.next()) {
                            ChadoFeature cf = new ChadoFeature(rs);
                            Item chromosome = getChadoDBConverter().createItem("Chromosome");
                            cf.populateBioEntity(chromosome, sourceOrganism);
                            store(chromosome);
                            sourceChromosomeMap.put(chromosome);
                            LOG.info("Stored new source chromosome "+chromosome.getAttribute("primaryIdentifier").getValue());
                        }
                        rs.close();
                    }
                    // create target chromosome, store and add to targetChromosomeMap if needed
                    String targetChromosome = gff.getTargetChromosome();
                    if (targetChromosome!=null) {
                        String targetName = targetChromosome.replace("Gm","glyma.Chr"); // chado chromosome name HARDCODED
                        if (!targetChromosomeMap.containsSecondaryIdentifier(targetName)) {
                            // create the chromosome Item from the chado record
                            query = "SELECT * FROM feature WHERE type_id=43403 AND organism_id="+target_organism_id+" AND name='"+targetName+"'";
                            rs = stmt.executeQuery(query);
                            if (rs.next()) {
                                ChadoFeature cf = new ChadoFeature(rs);
                                Item chromosome = getChadoDBConverter().createItem("Chromosome");
                                cf.populateBioEntity(chromosome, targetOrganism);
                                store(chromosome);
                                targetChromosomeMap.put(chromosome);
                                LOG.info("Stored new target chromosome "+chromosome.getAttribute("primaryIdentifier").getValue());
                            }
                            rs.close();
                        }
                    }
                    // store in map if it's a syntenic_region
                    if (gff.type.equals("syntenic_region")) {
                        gffMap.put(gff.attributeName, gff);
                        LOG.info("Synteny region: "+gff.attributeName);
                    }
                }
            }
            gffReader.close();
            
            LOG.info("Read "+gffMap.size()+" GFF records from "+dagChainerFile.getName()+".");
            LOG.info("Created "+sourceChromosomeMap.size()+" source Chromosome items from "+dagChainerFile.getName()+".");
            LOG.info("Created "+targetChromosomeMap.size()+" target Chromosome items from "+dagChainerFile.getName()+".");

        } catch (Exception ex) {

            LOG.error(ex.toString());
            ex.printStackTrace();
            System.exit(1);

        }

        // ----------------------------------------------------------------------
        // Now spin through the gffMap records and store the SyntenyBlock Items
        // ----------------------------------------------------------------------

        LOG.info("Creating, linking and storing synteny blocks...");

        for (GFFRecord gff : gffMap.values())  {

            // populate the source region and its location
            Item sourceChromosome = sourceChromosomeMap.getBySecondaryIdentifier(gff.seqid.replace("Pv","phavu.Chr")); // HARDCODED SOURCE GENOME
            Item sourceRegion = getChadoDBConverter().createItem("SyntenyRegion");
            Item sourceChromosomeLocation = getChadoDBConverter().createItem("Location");
            gff.populateSourceRegion(sourceRegion, sourceOrganism, sourceChromosome, sourceChromosomeLocation);

            // populate the target region and its location
            Item targetChromosome = targetChromosomeMap.getBySecondaryIdentifier(gff.getTargetChromosome().replace("Gm","glyma.Chr")); // HARDCODED TARGET GENOME
            Item targetRegion = getChadoDBConverter().createItem("SyntenyRegion");
            Item targetChromosomeLocation = getChadoDBConverter().createItem("Location");
            gff.populateTargetRegion(targetRegion, targetOrganism, targetChromosome, targetChromosomeLocation);

            // associate the two regions with a synteny block
            Item syntenyBlock = getChadoDBConverter().createItem("SyntenyBlock");
            syntenyBlock.setAttribute("primaryIdentifier", gff.attributeName);
            syntenyBlock.setAttribute("medianKs", String.valueOf(gff.attributeMedianKs));
            syntenyBlock.setReference("sourceRegion", sourceRegion);
            syntenyBlock.setReference("targetRegion", targetRegion);
            store(syntenyBlock);

	    sourceRegion.setReference("syntenyBlock", syntenyBlock);
	    targetRegion.setReference("syntenyBlock", syntenyBlock);

            store(sourceRegion);
            store(sourceChromosomeLocation);

            store(targetRegion);
            store(targetChromosomeLocation);
            
        }

        LOG.info("...done.");    

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

}
