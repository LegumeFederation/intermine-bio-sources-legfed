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
 * ONLY ONE ORGANISM/FILE PER INVOCATION.
 *
 * The taxon ID of the organism are given as project.xml parameters along with file locations:
 *
 *  <source name="phavu-blair-gff" type="legfed-geneticmarker-file">
 *    <property name="taxonId" value="3885"/>
 *    <property name="src.data.dir" value="/home/legfed/datastore/mixed.map1.blair.7PMp"/>
 *    <property name="gffFilename" value="phavu.mixed.map1.blair.7PMp.map.gff3"/>
 *  </source>
 *
 * @author Sam Hokin
 */
public class GeneticMarkerGFFConverter extends BioFileConverter {
	
    private static final Logger LOG = Logger.getLogger(GeneticMarkerGFFConverter.class);

    // project.xml parameters
    int taxonId = 0;
    String gffFilename = null;

    // stored in close()
    Map<String,Item> chromosomeMap = new HashMap<String,Item>();
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
     * Set the taxon ID from a project.xml parameter.
     */
    public void setTaxonId(String input) {
        taxonId = Integer.parseInt(input);
        LOG.info("Setting taxonId="+taxonId);
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

        // validation
        if (taxonId==0) throw new BuildException("taxonId property not set");
        if (gffFilename==null) throw new BuildException("gffFilename property not set");

        // only run for given Gff file
        if (getCurrentFile().getName().equals(gffFilename)) {
	
            LOG.info("Processing GFF file "+getCurrentFile().getName()+"...");

            // create and store the organism
            Item organism = createItem("Organism");
            organism.setAttribute("taxonId", String.valueOf(taxonId));
            store(organism);

            // run through the lines
            String line;
            BufferedReader gffReader = new BufferedReader(reader);
            while ((line=gffReader.readLine())!=null) {

                if (!line.startsWith("#")) {

                    GFF3Record gff = new GFF3Record(line);

                    // set marker name to be first of those listed in GFF record
                    List<String> names = gff.getNames();
                    String name = names.get(0);
		
                    // // check for dupe marker
                    // String key = taxonId+"_"+variety;
                    // if (markerMap.containsKey(name)) {
                    //     String otherKey = markerMap.get(name);
                    //     if (key.equals(otherKey)) {
		    //     LOG.info("Ignoring duplicate marker: "+name+" for organism "+key);
		    //     continue;
                    //     } else {
                    //         markerMap.put(name, key);
                    //     }
                    // } else {
                    //     markerMap.put(name, key);
                    // }
		
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
                        LOG.info("Created and stored chromosome: "+chrName);
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
            gffReader.close();
            
        }
	
    }

}
