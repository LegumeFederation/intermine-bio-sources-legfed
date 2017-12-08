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

import org.intermine.model.bio.Organism;

/**
 * Store the genetic marker and genomic data from a GFF file.
 *
 * The GFF file may contain duplicates - the last one read is the one that is loaded.
 *
 * Above the GFF data is a couple comment lines specifying the organism:
 *
 * #TaxonID     3920
 * #Variety     IT97K-499-35
 *
 * NOTE: the GFF file is presumed to NOT contain supercontigs.
 *
 * @author Sam Hokin
 */
public class GeneticMarkerGFFConverter extends BioFileConverter {
	
    private static final Logger LOG = Logger.getLogger(GeneticMarkerGFFConverter.class);

    // the items will all be stored in the close() method, to avoid dupes
    Map<String,Item> organismMap = new HashMap<String,Item>();
    Map<String,Item> chromosomeMap = new HashMap<String,Item>();

    // to check for dups
    Map<String,String> markerMap = new HashMap<String,String>();
    
    /**
     * Create a new GeneticMarkerGFFConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public GeneticMarkerGFFConverter(ItemWriter writer, Model model) {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     * Process the supplied GFF file to create GeneticMarker Items, along with Chromosome and Supercontig Items from the seqid column.
     */
    @Override
    public void process(Reader reader) throws Exception {

	if (getCurrentFile().getName().equals("README")) return;
	
        LOG.info("Processing GFF file "+getCurrentFile().getName()+"...");

        // persistent items per file
	int taxonId = 0;
	String variety = null;
	Item organism = null;
        
        String line;
        BufferedReader gffReader = new BufferedReader(reader);
        while ((line=gffReader.readLine()) != null) {

	    // create and store organism if we're ready
	    if (organism==null && taxonId>0 && variety!=null) {
		String key = taxonId+"_"+variety;
		if (organismMap.containsKey(key)) {
		    organism = organismMap.get(key);
		} else {
		    organism = createItem("Organism");
		    organism.setAttribute("taxonId", String.valueOf(taxonId));
		    organism.setAttribute("variety", variety);
		    store(organism);
		    organismMap.put(key, organism);
		    LOG.info("Created organism: "+taxonId+" ("+variety+")");
		}
	    }

	    if (line.startsWith("#TaxonID")) {

		String[] parts = line.split("\t");
		taxonId = Integer.parseInt(parts[1]);

	    } else if (line.startsWith("#Variety")) {

		String[] parts = line.split("\t");
		variety = parts[1];

	    } else if (line.startsWith("#") || line.trim().length()==0) {

		// do nothing, comment or blank line

	    } else {

		// bail if organism is not set
		if (organism==null) {
		    LOG.error("Organism not set: taxonId="+taxonId+", variety="+variety);
		    throw new RuntimeException("Organism not set: taxonId="+taxonId+", variety="+variety);
		}

                GFF3Record gff = new GFF3Record(line);

		// set marker name to be first of those listed in GFF record
		List<String> names = gff.getNames();
		String name = names.get(0);
		
		// check for dupe marke
		String key = taxonId+"_"+variety;
		if (markerMap.containsKey(name)) {
		    String otherKey = markerMap.get(name);
		    if (key.equals(otherKey)) {
			LOG.info("Ignoring duplicate marker: "+name+" for organism "+key);
			continue;
		    } else {
			markerMap.put(name, key);
		    }
		} else {
		    markerMap.put(name, key);
		}
		
		String chrName = gff.getSequenceID();
		Item chromosome;
		if (chromosomeMap.containsKey(chrName)) {
		    chromosome = chromosomeMap.get(chrName);
		} else {
		    chromosome = createItem("Chromosome");
		    chromosome.setAttribute("primaryIdentifier", chrName);
		    chromosome.setReference("organism", organism);
		    store(chromosome);
		    chromosomeMap.put(chrName, chromosome);
		    LOG.info("Created chromosome: "+chrName);
		}
		
                Item location = createItem("Location");
		Item marker = createItem("GeneticMarker");

		// populate and store the chromosome location
		location.setAttribute("start", String.valueOf(gff.getStart()));
		location.setAttribute("end", String.valueOf(gff.getEnd()));
		location.setReference("locatedOn", chromosome);
		location.setReference("feature", marker);
		store(location);

		// populate and store the genetic marker
                marker.setAttribute("primaryIdentifier", name);
                marker.setReference("organism", organism);
                marker.setAttribute("type", gff.getType());
                marker.setAttribute("length", String.valueOf(gff.getEnd()-gff.getStart()+1));
		marker.setReference("chromosome", chromosome);
		marker.setReference("chromosomeLocation", location);
		store(marker);

            }
        }
	
    }

}
