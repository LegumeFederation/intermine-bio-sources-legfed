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
 * ONE ORGANISM per FILE.
 *
 * The taxon ID of the organism are optionally given as project.xml parameters along with file locations:
 *
 *  <source name="phavu-blair-gff" type="legfed-geneticmarker-file">
 *    <property name="taxonId" value="3885"/>
 *    <property name="variety" value="McDonald123"/>
 *    <property name="src.data.dir" value="/home/legfed/datastore/mixed.map1.blair.7PMp"/>
 *    <property name="gffFilename" value="phavu.mixed.map1.blair.7PMp.map.gff3"/>
 *  </source>
 *
 * OR you can put the TaxonID and Variety in file comments:
 *
 * #TaxonID        130453
 * #Variety        V14167
 *
 * NOTE: this contains a hack to associate markers with Supercontig if the "chromosome" name contains "scaffold" (case-insensitively).
 *
 * @author Sam Hokin
 */
public class GeneticMarkerGFFConverter extends BioFileConverter {
	
    private static final Logger LOG = Logger.getLogger(GeneticMarkerGFFConverter.class);

    // organism parameters
    String taxonId;
    String variety;
    
    // optional file name for single-file operation
    String gffFilename;

    // stored to avoid dupes
    Map<String,Item> organismMap = new HashMap<>();
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
     * Set the taxon ID from a project.xml parameter.
     */
    public void setTaxonId(String input) {
        taxonId = input;
        LOG.info("Setting taxonId="+taxonId+" from project.xml.");
    }

    /**
     * Set the variety from a project.xml parameter.
     */
    public void setVariety(String input) {
        variety = input;
        LOG.info("Setting variety="+variety+" from project.xml.");
    }

    /**
     * Set the GFF filename from a project.xml parameter.
     */
    public void setGffFilename(String filename) {
        gffFilename = filename;
        LOG.info("Setting gffFilename="+gffFilename);
    }

    /**
     * {@inheritDoc}
     * Process the supplied GFF file to create GeneticMarker Items, along with Chromosome and Supercontig Items from the seqid column.
     */
    @Override
    public void process(Reader reader) throws Exception {

        // this file's organism
        Item organism = null;
        
        // bail if the specified file isn't found
        if (gffFilename!=null && !getCurrentFile().getName().equals(gffFilename)) throw new RuntimeException("Specified single GFF file "+gffFilename+" not found.");

        // run through the lines
        LOG.info("Processing GFF file "+getCurrentFile().getName()+"...");
        String line;
        BufferedReader gffReader = new BufferedReader(reader);
        while ((line=gffReader.readLine())!=null) {

            // header stuff
            if (line.toLowerCase().startsWith("#taxonid")) {

                String[] parts = line.split("\t");
                taxonId = parts[1];
                LOG.info("Setting taxonId="+taxonId+" from GFF file header.");

            } else if (line.toLowerCase().startsWith("#variety")) {

                String[] parts = line.split("\t");
                variety = parts[1];
                LOG.info("Setting variety="+variety+" from GFF file header.");
                
            } else if (!line.startsWith("#")) {

                if (taxonId==null) throw new RuntimeException("Taxon ID not set, not reading GFF data.");

                if (organism==null) {
                    // create the organism, add to the map
                    String key = taxonId;
                    if (variety!=null) key += "_"+variety;
                    if (organismMap.containsKey(key)) {
                        organism = organismMap.get(key);
                    } else {
                        organism = createItem("Organism");
                        organism.setAttribute("taxonId", taxonId);
                        if (variety!=null) organism.setAttribute("variety", String.valueOf(variety));
                        store(organism);
                        organismMap.put(key, organism);
                    }
                }

                GFF3Record gff = new GFF3Record(line);

                // set marker name to be first of those listed in GFF record
                List<String> names = gff.getNames();
                String name = names.get(0);
		
                // check for dupe marker, bail if dupe for this organism
                String markerKey = taxonId+"_"+variety+":"+name;
                if (markerSet.contains(markerKey)) {
                    LOG.info("Ignoring duplicate marker: "+markerKey);
                    continue;
                } else {
                    markerSet.add(markerKey);
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
                    store(sequence);
                    sequenceMap.put(sequenceName, sequence);
                    LOG.info("Created and stored "+sequenceName);
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
