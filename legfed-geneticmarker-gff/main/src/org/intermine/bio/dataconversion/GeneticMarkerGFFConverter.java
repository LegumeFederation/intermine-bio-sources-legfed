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
     * Process the supplied GFF file to create GeneticMarker Items, along with Chromosome and Supercontig Items from the seqid column.
     */
    public void process(Reader reader) throws Exception {

        LOG.info("Processing GFF file "+getCurrentFile().getName()+"...");

        // store everything in maps for uniqueness, store at the end; all keyed by name
        Map<String,Item> chromosomeMap = new HashMap<String,Item>();
        Map<String,Item> supercontigMap = new HashMap<String,Item>();
        Map<String,Item> markerMap = new HashMap<String,Item>();
        Map<String,Item> locationMap = new HashMap<String,Item>();

        // create the organism Item
        int taxonId = getTaxonId();
        Item organism = createItem("Organism");
        BioStoreHook.setSOTerm(this, organism, "organism", getSequenceOntologyRefId());
        organism.setAttribute("taxonId", String.valueOf(taxonId));

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
        LOG.info("Created "+chromosomeMap.size()+" Chromosome items.");
        LOG.info("Created "+supercontigMap.size()+" Supercontig items.");

        // now store all the items
        store(organism);
        for (Item item : chromosomeMap.values()) store(item);
        for (Item item : supercontigMap.values()) store(item);
        for (Item item : locationMap.values()) store(item);
        for (Item item : markerMap.values()) store(item);

    }
    
}
