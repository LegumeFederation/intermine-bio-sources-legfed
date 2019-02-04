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

import java.util.Map;
import java.util.HashMap;

import org.apache.log4j.Logger;

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.xml.full.Item;

/**
 * Store genetic marker genomic positions, type and optional motif from a tab-delimited file. TaxonID/Strain Name are entered as shown.
 * Motif is optional.
 *
 * TaxonID 3847
 * Strain Williams82
 * #primaryIdentifier secondaryIdentifier Type Chromosome Start   End     Motif
 * Sat_332            BARCSOYSSR_01_0019  SSR  Gm01       355597  355646  (AT)25
 *
 * @author Sam Hokin, NCGR
 */
public class MarkerChromosomeFileConverter extends BioFileConverter {
	
    private static final Logger LOG = Logger.getLogger(MarkerChromosomeFileConverter.class);

    // maps contain Items that are repeated across files
    Map<String,Item> chromosomeMap = new HashMap<String,Item>();
    Map<String,Item> organismMap = new HashMap<String,Item>();
    Map<String,Item> strainMap = new HashMap<String,Item>();

    /**
     * Create a new MarkerChromosomeFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public MarkerChromosomeFileConverter(ItemWriter writer, Model model) {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     * Process the marker-chromosome relationships by reading in from a tab-delimited file.
     */
    @Override
    public void process(Reader reader) throws Exception {

        // don't process README files
        if (getCurrentFile().getName().contains("README")) return;

        LOG.info("Processing file "+getCurrentFile().getName()+"...");

        // organism/strain for this file
        Item organism = null;
	Item strain = null;

        BufferedReader markerReader = new BufferedReader(reader);
	String line;
        while ((line=markerReader.readLine())!=null) {
            if (!line.startsWith("#") && line.trim().length()>0) {
                String[] parts = line.split("\t");
                String key = parts[0];
                String value = parts[1];
		
                // organism and strain lines at top
                if (key.toLowerCase().equals("taxonid")) {
                    String taxonId = parts[1];
                    if (organismMap.containsKey(taxonId)) {
                        organism = organismMap.get(taxonId);
                    } else {
                        organism = createItem("Organism");
                        organism.setAttribute("taxonId", taxonId);
                        store(organism);
                        LOG.info("Stored organism:"+taxonId);
                        organismMap.put(taxonId, organism);
                    }
                    
                } else if (key.toLowerCase().equals("strain")) {
                    String strainName = parts[1];
                    if (strainMap.containsKey(strainName)) {
                        strain = strainMap.get(strainName);
                    } else {
                        strain = createItem("Strain");
                        strain.setAttribute("primaryIdentifier", strainName);
                        strain.setReference("organism", organism);
                        store(strain);
                        LOG.info("Stored strain:"+strainName);
                        strainMap.put(strainName, strain);
                    }

                } else {
		    
                    // check that we've got an organism and strain
                    if (organism==null) {
                        LOG.error("Organism has not been formed for marker/chromosome import in file "+getCurrentFile().getName());
                        throw new RuntimeException("Organism has not been formed for marker/chromosome import in file "+getCurrentFile().getName());
                    }
                    if (strain==null) {
                        LOG.error("Strain has not been formed for marker/chromosome import in file "+getCurrentFile().getName());
                        throw new RuntimeException("Strain has not been formed for marker/chromosome import in file "+getCurrentFile().getName());
                    }
		    
                    MarkerChromosomeRecord rec = new MarkerChromosomeRecord(line);
                    
                    // create the marker
                    Item marker = createItem("GeneticMarker");
                    marker.setAttribute("primaryIdentifier", rec.primaryIdentifier);
                    marker.setReference("organism", organism);
		    marker.setReference("strain", strain);
                    if (rec.secondaryIdentifier.length()>0) marker.setAttribute("secondaryIdentifier", rec.secondaryIdentifier);
                    if (rec.type.length()>0) marker.setAttribute("type", rec.type);
                    if (rec.motif!=null && rec.motif.length()>0) marker.setAttribute("motif", rec.motif);
                    // set the chromosome or supercontig reference
                    // HACK: assume supercontig has "scaffold" or "contig" in the name
                    boolean isSupercontig = (rec.chromosome.toLowerCase().contains("scaffold") || rec.chromosome.toLowerCase().contains("contig"));
                    Item chromosome = null;
                    if (chromosomeMap.containsKey(rec.chromosome)) {
                        chromosome = chromosomeMap.get(rec.chromosome);
                    } else {
                        // create and store this chromosome/supercontig
                        if (isSupercontig) {
                            chromosome = createItem("Supercontig");
                        } else {
                            chromosome = createItem("Chromosome");
                        }
                        chromosome.setAttribute("primaryIdentifier", rec.chromosome);
                        chromosome.setReference("organism", organism);
			chromosome.setReference("strain", strain);
                        store(chromosome);
                        chromosomeMap.put(rec.chromosome, chromosome);
                        LOG.info("Created and stored chromosome/supercontig: "+rec.chromosome);
                    }
                    if (isSupercontig) {
                        marker.setReference("supercontig", chromosome);
                    } else {
                        marker.setReference("chromosome", chromosome);
                    }
                    
                    // create and store the Location object
                    Item location = createItem("Location");
                    location.setReference("locatedOn", chromosome);
                    location.setReference("feature", marker);
                    location.setAttribute("start", String.valueOf(rec.start));
                    location.setAttribute("end", String.valueOf(rec.end));
                    location.setAttribute("strand", String.valueOf(+1));
                    store(location);
                    
                    // set the chromosomeLocation/supercontigLocation reference and store the marker
                    if (isSupercontig) {
                        marker.setReference("supercontigLocation", location);
                    } else {
                        marker.setReference("chromosomeLocation", location);
                    }                    
                    store(marker);

                }
	    }
	}
        
	markerReader.close();
    }
}
