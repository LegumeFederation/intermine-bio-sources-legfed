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
public class SyntenyFileProcessor extends LegfedFileProcessor {
	
    private static final Logger LOG = Logger.getLogger(SyntenyFileProcessor.class);

    /**
     * Create a new SyntenyFileProcessor
     * @param legfedFileConverter the LegfedFileConverter that is controlling this processor
     */
    public SyntenyFileProcessor(LegfedFileConverter legfedFileConverter) {
        super(legfedFileConverter);
    }

    /**
     * {@inheritDoc}
     * We process the GFF file by creating SyntenyBlock and SyntenicRegion items and storing them.
     */
    @Override
    public void process(Reader reader) throws ObjectStoreException {

        // ---------------------------------------------------------
        // -- Get the source and target organisms from the file name
        // ---------------------------------------------------------

        // InterMine organism IDs
        int sourceTaxonId = 0;
        int targetTaxonId = 0;

        // parse the DAGchainer file name for the taxonomy IDs corresponding to the source and target
        try {
            String dagChainerFileName = getLegfedFileConverter().getCurrentFileName();
            String[] bothChunks = dagChainerFileName.split("_");
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

        // create and store organism Items
        Item sourceOrganism = getLegfedFileConverter().createItem("Organism");
        BioStoreHook.setSOTerm(getLegfedFileConverter(), sourceOrganism, "organism", getLegfedFileConverter().getSequenceOntologyRefId());
        sourceOrganism.setAttribute("taxonId", String.valueOf(sourceTaxonId));
        store(sourceOrganism);
        Item targetOrganism = getLegfedFileConverter().createItem("Organism");
        BioStoreHook.setSOTerm(getLegfedFileConverter(), targetOrganism, "organism", getLegfedFileConverter().getSequenceOntologyRefId());
        targetOrganism.setAttribute("taxonId", String.valueOf(targetTaxonId));
        store(targetOrganism);
        LOG.info("Created and stored sourceOrganism, taxonId="+sourceTaxonId+" and targetOrganism, taxonId="+targetTaxonId);
        
        // -------------------------------------------------------------------------------------------------------
        // Load the GFF data into a map. Create source and target chromosome Maps keyed by name=primaryIdentifier.
        // -------------------------------------------------------------------------------------------------------

        Map<String,GFFRecord> gffMap = new HashMap<String,GFFRecord>();
        Map<String,Item> sourceChromosomeMap = new HashMap<String,Item>();
        Map<String,Item> targetChromosomeMap = new HashMap<String,Item>();
        
        try {
            
            BufferedReader gffReader = new BufferedReader(reader);
            LOG.info("Reading DAGChainer GFF file...");
            String gffLine = null;
            while ((gffLine=gffReader.readLine()) != null) {
                GFFRecord gff = new GFFRecord(gffLine);
                String sourceName = gff.seqid; // could potentially alter to match the chado values here
                String targetName = gff.getTargetChromosome(); // could potentially alter to match the chado values here
                if (gff.attributeName!=null) {
                    if (sourceName!=null && !sourceChromosomeMap.containsKey(sourceName)) {
                        // create the chromosome Item and add to the source map
                        Item chromosome = getLegfedFileConverter().createItem("Chromosome");
                        BioStoreHook.setSOTerm(getLegfedFileConverter(), chromosome, "chromosome", getLegfedFileConverter().getSequenceOntologyRefId());
                        chromosome.setAttribute("primaryIdentifier", sourceName);
                        store(chromosome);
                        sourceChromosomeMap.put(sourceName, chromosome);
                        LOG.info("Stored new source chromosome:"+sourceName);
                    }
                    if (targetName!=null && !targetChromosomeMap.containsKey(targetName)) {
                        // create the chromosome Item and add to the target map
                        Item chromosome = getLegfedFileConverter().createItem("Chromosome");
                        BioStoreHook.setSOTerm(getLegfedFileConverter(), chromosome, "chromosome", getLegfedFileConverter().getSequenceOntologyRefId());
                        chromosome.setAttribute("primaryIdentifier", targetName);
                        store(chromosome);
                        targetChromosomeMap.put(targetName, chromosome);
                        LOG.info("Stored new target chromosome:"+targetName);
                    }
                    // store synteny blocks in gff map; form unique (we hope) key
                    if (gff.type.equals("syntenic_region")) {
                        String syntenyBlockID = gff.attributeName+"."+gff.attributeID;
                        gffMap.put(syntenyBlockID, gff);
                    }
                }
            }
            gffReader.close();
            
            LOG.info("Read "+gffMap.size()+" syntenic_region GFF records.");
            LOG.info("Created "+sourceChromosomeMap.size()+" source Chromosome items.");
            LOG.info("Created "+targetChromosomeMap.size()+" target Chromosome items.");

        } catch (Exception ex) {

            LOG.error(ex.toString());
            ex.printStackTrace();
            System.exit(1);

        }

        // ----------------------------------------------------------------------
        // Now spin through the gffMap records and store the SyntenyBlock Items
        // ----------------------------------------------------------------------

        LOG.info("Creating, linking and storing synteny blocks...");

        for (Map.Entry<String,GFFRecord> entry : gffMap.entrySet()) {
            
            String syntenyBlockID = entry.getKey();
            GFFRecord gff = entry.getValue();
            
            String sourceName = gff.seqid;
            String targetName = gff.getTargetChromosome();

            // populate the source region and its location
            Item sourceChromosome = sourceChromosomeMap.get(sourceName);
            Item sourceRegion = getLegfedFileConverter().createItem("SyntenicRegion");
            BioStoreHook.setSOTerm(getLegfedFileConverter(), sourceRegion, "syntenic_region", getLegfedFileConverter().getSequenceOntologyRefId());
            Item sourceChromosomeLocation = getLegfedFileConverter().createItem("Location");
            gff.populateSourceRegion(sourceRegion, sourceOrganism, sourceChromosome, sourceChromosomeLocation);
            
            // populate the target region and its location
            Item targetChromosome = targetChromosomeMap.get(targetName);
            Item targetRegion = getLegfedFileConverter().createItem("SyntenicRegion");
            BioStoreHook.setSOTerm(getLegfedFileConverter(), targetRegion, "syntenic_region", getLegfedFileConverter().getSequenceOntologyRefId());
            Item targetChromosomeLocation = getLegfedFileConverter().createItem("Location");
            gff.populateTargetRegion(targetRegion, targetOrganism, targetChromosome, targetChromosomeLocation);

            // associate the two regions with a synteny block
            Item syntenyBlock = getLegfedFileConverter().createItem("SyntenyBlock");
            LOG.info("Storing SyntenyBlock:"+syntenyBlockID);
            syntenyBlock.setAttribute("primaryIdentifier", syntenyBlockID);
            syntenyBlock.setAttribute("medianKs", String.valueOf(gff.attributeMedianKs));
            syntenyBlock.setReference("sourceRegion", sourceRegion);
            syntenyBlock.setReference("targetRegion", targetRegion);
            store(syntenyBlock);

	    sourceRegion.setReference("syntenyBlock", syntenyBlock);
            store(sourceRegion);
            store(sourceChromosomeLocation);

	    targetRegion.setReference("syntenyBlock", syntenyBlock);
            store(targetRegion);
            store(targetChromosomeLocation);
            
        }

        LOG.info("...done.");    

    }

}
