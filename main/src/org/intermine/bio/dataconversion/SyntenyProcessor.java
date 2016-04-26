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
 * Read synteny blocks for two organisms from a DAGchainer synteny GFF file and store them as SyntenyBlock items,
 * each related to a source and target SyntenicRegion. This is designed to use the GFF annotation produced by DAGchainer.
 *
 * The source and target organism taxIDs are taken from the DAGchainer GFF file name.
 * For example, if source=Phaseolus vulgaris and target=Glycine max, then the file will be, e.g. synteny.3885_with_3847.gff3, or
 * any other file name which contains ".3885_chars_3847.gff" and no other dots or underscores.
 *
 * @author Sam Hokin, NCGR
 */
public class SyntenyProcessor extends ChadoProcessor {
	
    private static final Logger LOG = Logger.getLogger(SyntenyProcessor.class);

    /**
     * Create a new SyntenyProcessor
     * @param chadoDBConverter the ChadoDBConverter that is controlling this processor
     */
    public SyntenyProcessor(ChadoDBConverter chadoDBConverter) {
        super(chadoDBConverter);
    }

    /**
     * {@inheritDoc}
     * We process the GFF file by creating SyntenicRegion items and storing them.
     * In order to do so, we need to access the two organisms' chromosomes from the chado database.
     */
    @Override
    public void process(Connection connection) throws SQLException, ObjectStoreException {

        // ---------------------------------------------------------
        // ---------------- INITIAL DATA LOADING -------------------
        // ---------------------------------------------------------

        // initialize our DB statement and stuff for chado queries
        Statement stmt = connection.createStatement();

        // InterMine organism IDs
        int sourceTaxonId = 0;
        int targetTaxonId = 0;

        // chado organism_ids
        int source_organism_id = 0;
        int target_organism_id = 0;

        // DAGchainer synteny gff file location is passed in as parameter dagChainerFile in project.xml
        File dagChainerFile = getChadoDBConverter().getDagChainerFile();
        LOG.info("dagChainerFile="+dagChainerFile.getAbsolutePath());

        // parse the taxonomy IDs corresponding to the source and target from the DAGchainer file name
        LOG.info("Attempting to get dagChainerFile...");
        String dagChainerFilename = dagChainerFile.getName();
        LOG.info("dagChainerFilename="+dagChainerFilename);
        try {
            String[] bothChunks = dagChainerFilename.split("_");
            String[] sourceChunks = bothChunks[0].split("\\.");
            String[] targetChunks = bothChunks[2].split("\\.");
            String sourceString = sourceChunks[1];
            String targetString = targetChunks[0];
            LOG.info("source taxID="+sourceString+" target taxID="+targetString);
            sourceTaxonId = Integer.parseInt(sourceString);
            targetTaxonId = Integer.parseInt(targetString);
        } catch (Exception ex) {
            LOG.error(ex.toString());
            throw new RuntimeException(ex);
        }

        // get the chado organism_ids, which are keys of the ChadoIdToOrgDataMap;
        try {
            for (Map.Entry<Integer,OrganismData> entry : getChadoDBConverter().getChadoIdToOrgDataMap().entrySet()) {
                if (entry.getValue().getTaxonId()==sourceTaxonId) source_organism_id = (int) entry.getKey();
                if (entry.getValue().getTaxonId()==targetTaxonId) target_organism_id = (int) entry.getKey();
            }
        } catch (Exception ex) {
            LOG.error("Error determining source and target chado organism_ids: "+ex.toString());
            throw new RuntimeException(ex);
        }
        LOG.info("source_organism_id="+source_organism_id);
        LOG.info("target_organism_id="+target_organism_id);
        if (source_organism_id==0 || target_organism_id==0) {
            LOG.error("Error retrieving source and target chado organism_ids from project.xml.");
            throw new RuntimeException("Error retrieving source and target chado organism_ids from project.xml.");
        }

        // create and store organism Items
        Item sourceOrganism = getChadoDBConverter().createItem("Organism");
        BioStoreHook.setSOTerm(getChadoDBConverter(), sourceOrganism, "organism", getChadoDBConverter().getSequenceOntologyRefId());
        sourceOrganism.setAttribute("taxonId", String.valueOf(sourceTaxonId));
        store(sourceOrganism);
        Item targetOrganism = getChadoDBConverter().createItem("Organism");
        BioStoreHook.setSOTerm(getChadoDBConverter(), targetOrganism, "organism", getChadoDBConverter().getSequenceOntologyRefId());
        targetOrganism.setAttribute("taxonId", String.valueOf(targetTaxonId));
        store(targetOrganism);
        LOG.info("Created and stored sourceOrganism, taxonId="+sourceTaxonId+" and targetOrganism, taxonId="+targetTaxonId);
        
        // -------------------------------------------------------------------------------------------------------------------
        // Load the GFF data into a map. Create source and target chromosome Items.
        // -------------------------------------------------------------------------------------------------------------------

        Map<String,GFFRecord> gffMap = new HashMap<String,GFFRecord>();
        ItemMap sourceChromosomeMap = new ItemMap();
        ItemMap targetChromosomeMap = new ItemMap();
        
        try {
            
            // get chromosome cvterm_id
            ResultSet rs0 = stmt.executeQuery("SELECT cvterm_id FROM cvterm WHERE name='chromosome' AND definition LIKE 'Structural unit%'");
            rs0.next();
            int chrCVTermId = rs0.getInt("cvterm_id");
            rs0.close();
            
            BufferedReader gffReader = new BufferedReader(new FileReader(dagChainerFile));
            LOG.info("Reading GFF file:"+dagChainerFile.getName());
            String gffLine = null;
            while ((gffLine=gffReader.readLine()) != null) {
                GFFRecord gff = new GFFRecord(gffLine);
                if (gff.seqid!=null && gff.attributeName!=null) {
                    String sourceName = gff.seqid; // can tweak sourceName to be different from seqid to match chado
                    // add to the map
                    if (!sourceChromosomeMap.containsPrimaryIdentifier(sourceName)) {
                        // create the chromosome Item from the chado record
                        ResultSet rs = stmt.executeQuery("SELECT * FROM feature WHERE type_id="+chrCVTermId+" AND organism_id="+source_organism_id+" AND name='"+sourceName+"'");
                        if (rs.next()) {
                            ChadoFeature cf = new ChadoFeature(rs);
                            Item chromosome = getChadoDBConverter().createItem("Chromosome");
                            BioStoreHook.setSOTerm(getChadoDBConverter(), chromosome, "chromosome", getChadoDBConverter().getSequenceOntologyRefId());
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
                        String targetName = targetChromosome; // can tweak targetName to be different from targetChromosome to match chado
                        if (!targetChromosomeMap.containsPrimaryIdentifier(sourceName)) {
                            // create the chromosome Item from the chado record
                            ResultSet rs = stmt.executeQuery("SELECT * FROM feature WHERE type_id="+chrCVTermId+" AND organism_id="+target_organism_id+" AND name='"+targetName+"'");
                            if (rs.next()) {
                                ChadoFeature cf = new ChadoFeature(rs);
                                Item chromosome = getChadoDBConverter().createItem("Chromosome");
                                BioStoreHook.setSOTerm(getChadoDBConverter(), chromosome, "chromosome", getChadoDBConverter().getSequenceOntologyRefId());
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
                        LOG.info("Syntenic region: "+gff.attributeName);
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
            
            String srcChrName = gff.seqid;
            String tgtChrName = gff.getTargetChromosome();

            // populate the source region and its location
            Item sourceChromosome = sourceChromosomeMap.getByPrimaryIdentifier(srcChrName);
            Item sourceRegion = getChadoDBConverter().createItem("SyntenicRegion");
            BioStoreHook.setSOTerm(getChadoDBConverter(), sourceRegion, "syntenic_region", getChadoDBConverter().getSequenceOntologyRefId());
            Item sourceChromosomeLocation = getChadoDBConverter().createItem("Location");
            gff.populateSourceRegion(sourceRegion, sourceOrganism, sourceChromosome, sourceChromosomeLocation);
            
            // populate the target region and its location
            Item targetChromosome = targetChromosomeMap.getByPrimaryIdentifier(tgtChrName);
            Item targetRegion = getChadoDBConverter().createItem("SyntenicRegion");
            BioStoreHook.setSOTerm(getChadoDBConverter(), targetRegion, "syntenic_region", getChadoDBConverter().getSequenceOntologyRefId());
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
