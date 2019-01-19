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

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.HashMap;
import java.util.HashSet;

import org.apache.log4j.Logger;

import org.intermine.bio.io.gff3.GFF3Record;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.xml.full.Item;

import org.intermine.model.bio.Organism;

/**
 * Store the genetic marker and genomic data from a GFF file.
 *
 * ONE ORGANISM+STRAIN per FILE.
 *
 *  <source name="phavu-blair-gff" type="legfed-geneticmarker-file">
 *    <property name="src.data.dir" value="/home/legfed/datastore/mixed.map1.blair.7PMp"/>
 *    <property name="gffFilename" value="phavu.mixed.map1.blair.7PMp.map.gff3"/>
 *  </source>
 *
 * Put the taxon ID and strain at the top in file comments:
 *
 * #TaxonID        130453
 * #Strain         V14167
 *
 * HACK: this contains a hack to associate markers with Supercontig if the "chromosome" name contains "scaffold" (case-insensitively).
 *
 * @author Sam Hokin
 */
public class GeneticMarkerGFFConverter extends BioFileConverter {
	
    private static final Logger LOG = Logger.getLogger(GeneticMarkerGFFConverter.class);

    // optional file name for single-file operation
    String gffFilename;

    // stored to avoid dupes
    Map<String,Item> organismMap = new HashMap<>();
    Map<String,Item> strainMap = new HashMap<>();
    Map<String,Item> sequenceMap = new HashMap<>();
    Set<String> markerSet = new HashSet<>();
    
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

        // this file's organism and strain
        Item organism = null;
	Item strain = null;
        
        // run through the lines
        LOG.info("Processing GFF file "+getCurrentFile().getName()+"...");
        String line;
        BufferedReader gffReader = new BufferedReader(reader);
        while ((line=gffReader.readLine())!=null) {

	    String[] parts = line.split("\t");

            // header stuff
            if (parts[0].toLowerCase().equals("#taxonid")) {
		String taxonId = parts[1];
		if (organismMap.containsKey(taxonId)) {
		    organism = organismMap.get(taxonId);
		} else {
		    organism = createItem("Organism");
		    organism.setAttribute("taxonId", taxonId);
		    store(organism);
		    organismMap.put(taxonId, organism);
		    LOG.info("Stored organism "+taxonId);
		}

            } else if (parts[0].toLowerCase().equals("#strain")) {
                String strainName = parts[1];
		if (strainMap.containsKey(strainName)) {
		    strain = strainMap.get(strainName);
		} else {
		    strain = createItem("Strain");
		    strain.setAttribute("primaryIdentifier", strainName);
		    strain.setReference("organism", organism);
		    store(strain);
		    LOG.info("Stored strain "+strainName);
		}
                
            } else if (line.startsWith("#")) {

		continue; // comment

	    } else {

		// check organism and strain given
		if (organism==null) {
		    LOG.error("Organism not set.");
		    throw new RuntimeException("Organism not set.");
		}
		if (strain==null) {
		    LOG.error("Strain not set.");
		    throw new RuntimeException("Strain not set.");
                }

                GFF3Record gff = new GFF3Record(line);

                // set marker name to be first of those listed in GFF record
                List<String> names = gff.getNames();
                String name = names.get(0);
		
                // check for dupe marker, bail if dupe for this organism
                if (markerSet.contains(name)) {
                    LOG.info("Ignoring duplicate marker: "+name);
                    continue;
                } else {
                    markerSet.add(name);
                }
		
                Item sequence;
                String sequenceName = gff.getSequenceID();
                boolean isSupercontig = sequenceName.toLowerCase().contains("scaffold");
                if (isSupercontig) sequenceName = sequenceName.toLowerCase();
                if (sequenceMap.containsKey(sequenceName)) {
                    sequence = sequenceMap.get(sequenceName);
                } else {
                    if (isSupercontig) {
                        sequence = createItem("Supercontig");
                        sequence.setAttribute("primaryIdentifier", sequenceName); // we insist on it being "scaffold"
		    } else {
                        sequence = createItem("Chromosome");
                        sequence.setAttribute("primaryIdentifier", sequenceName);
                    }
                    sequence.setReference("organism", organism);
		    sequence.setReference("strain", strain);
                    store(sequence);
                    sequenceMap.put(sequenceName, sequence);
                    LOG.info("Stored sequence "+sequenceName);
                }
		
                Item location = createItem("Location");
                Item marker = createItem("GeneticMarker");
                
                // populate and store the location
                location.setAttribute("start", String.valueOf(gff.getStart()));
                location.setAttribute("end", String.valueOf(gff.getEnd()));
                location.setReference("locatedOn", sequence);
                location.setReference("feature", marker);
                store(location);
                
                // populate and store the genetic marker
                marker.setAttribute("primaryIdentifier", name);
                marker.setReference("organism", organism);
		marker.setReference("strain", strain);
                marker.setAttribute("type", gff.getType());
                marker.setAttribute("length", String.valueOf(gff.getEnd()-gff.getStart()+1));
                if (isSupercontig) {
                    marker.setReference("supercontig", sequence);
                    marker.setReference("supercontigLocation", location);
                } else {
                    marker.setReference("chromosome", sequence);
                    marker.setReference("chromosomeLocation", location);
                }
                store(marker);
            }
        }
        gffReader.close();
    }
}
