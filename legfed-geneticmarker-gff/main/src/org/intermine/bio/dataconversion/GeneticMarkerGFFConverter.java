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
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import org.intermine.bio.io.gff3.GFF3Record;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.xml.full.Item;

/**
 * Store the genetic marker and genomic data from a GFF file.
 *
 * The GFF file may contain duplicates - the last one read is the one that is loaded.
 *
 * The GFF filename gives the taxon ID of the organism. This allows you to have several GFF files in a single src.data.dir 
 * that are run in a single invocation of this converter. The format is: anything-other-than-underscore_3845.gff3.
 *
 * @author Sam Hokin
 */
public class GeneticMarkerGFFConverter extends BioFileConverter {
	
    private static final Logger LOG = Logger.getLogger(GeneticMarkerGFFConverter.class);

    // the items will all be stored in the close() method, to avoid dupes
    Set<Item> organismSet = new HashSet<Item>();
    Map<String,Item> chromosomeMap = new HashMap<String,Item>();
    Map<String,Item> supercontigMap = new HashMap<String,Item>();
    Map<String,Item> markerMap = new HashMap<String,Item>();
    Map<String,Item> locationMap = new HashMap<String,Item>();
    
    /**
     * Create a new GeneticMarkerGFFConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public GeneticMarkerGFFConverter(ItemWriter writer, Model model) {
        super(writer, model);
    }

    /**
     * Get the taxon ID from the current file name, e.g. soybean_3847.gff3
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
     * Process the supplied GFF file to create GeneticMarker Items, along with Chromosome and Supercontig Items from the seqid column.
     */
    @Override
    public void process(Reader reader) throws Exception {

        LOG.info("Processing GFF file "+getCurrentFile().getName()+"...");

        // create the organism Item and add it to its map
        Item organism = createItem("Organism");
        BioStoreHook.setSOTerm(this, organism, "organism", getSequenceOntologyRefId());
        organism.setAttribute("taxonId", getTaxonId());
        organismSet.add(organism);

        // -------------------------------------------------------------------------------------------------------------------
        // Load the chromosomes, supercontigs and genetic markers from the GFF file
        // -------------------------------------------------------------------------------------------------------------------
        
        String line;
        BufferedReader gffReader = new BufferedReader(reader);
        while ((line=gffReader.readLine()) != null) {
            if (!line.startsWith("#")) {
                GFF3Record gff = new GFF3Record(line);
                Item marker = createItem("GeneticMarker");
                Item location = createItem("Location");
                // standard chromosome naming convention, e.g. "glyma.Chr01"
                if (gff.getSequenceID().contains("Chr")) {
                    // create and store this chromosome if not already in map; else get it from the map
                    String chrName = gff.getSequenceID();
                    Item chromosome = null;
                    if (chromosomeMap.containsKey(chrName)) {
                        chromosome = chromosomeMap.get(chrName);
                    } else {
                        chromosome = createItem("Chromosome");
                        BioStoreHook.setSOTerm(this, chromosome, "chromosome", getSequenceOntologyRefId());
                        chromosome.setAttribute("primaryIdentifier", chrName);
                        chromosomeMap.put(chrName, chromosome);
                    }
                    // populate the chromosome location
                    location.setAttribute("start", String.valueOf(gff.getStart()));
                    location.setAttribute("end", String.valueOf(gff.getEnd()));
                    location.setReference("feature", marker);
                    location.setReference("locatedOn", chromosome);
                    // associate the genetic marker with this chromosome/location
                    marker.setReference("chromosome", chromosome);
                    marker.setReference("chromosomeLocation", location);
                } else {
                    // create and store this supercontig if not already in map; else get it from the map
                    String supercontigName = gff.getSequenceID();
                    Item supercontig = null;
                    if (supercontigMap.containsKey(supercontigName)) {
                        supercontig = supercontigMap.get(supercontigName);
                    } else {
                        supercontig = createItem("Supercontig");
                        BioStoreHook.setSOTerm(this, supercontig, "supercontig", getSequenceOntologyRefId());
                        supercontig.setAttribute("primaryIdentifier", supercontigName);
                        supercontigMap.put(supercontigName, supercontig);
                    }
                    // populate the supercontig location
                    location.setAttribute("start", String.valueOf(gff.getStart()));
                    location.setAttribute("end", String.valueOf(gff.getEnd()));
                    location.setReference("feature", marker);
                    location.setReference("locatedOn", supercontig);
                    // associate the genetic marker with this supercontig/location
                    marker.setReference("supercontig", supercontig);
                    marker.setReference("supercontigLocation", location);
                }
                // set other marker attributes
                BioStoreHook.setSOTerm(this, marker, "genetic_marker", getSequenceOntologyRefId());
                marker.setReference("organism", organism);
                List<String> names = gff.getNames();
                String name = names.get(0);
                marker.setAttribute("primaryIdentifier", name);
                marker.setAttribute("type", gff.getType());
                marker.setAttribute("length", String.valueOf(gff.getEnd()-gff.getStart()+1));
                // add this marker and location to the maps; overwrites previous if same name
                markerMap.put(name, marker);
                locationMap.put(name, location);
            }
        }

        LOG.info("Created "+markerMap.size()+" distinct GeneticMarker items.");

    }

    /**
     * Store the items we've collected from the GFF files
     */
    @Override
    public void close() throws Exception {

        LOG.info("Storing "+organismSet.size()+" organism items...");
	for (Item organism : organismSet) store(organism);
        
	LOG.info("Storing "+chromosomeMap.size()+" chromosome items...");
	for (Item chromosome : chromosomeMap.values()) store(chromosome);

        LOG.info("Storing "+supercontigMap.size()+" supercontig items...");
        for (Item supercontig : supercontigMap.values()) store(supercontig);

        LOG.info("Storing "+markerMap.size()+" marker and location items...");
        for (Item item : markerMap.values()) store(item);
        for (Item item : locationMap.values()) store(item);

    }
    
}
