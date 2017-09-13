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
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

/**
 * Store SNP marker genomic positions from a VCF file.
 *
 * <pre>
 * CHROM  POS	ID	REF	ALT	QUAL	FILTER	INFO
 * Vu01	  74363	2_37329	A	G	999	.	DP=378
 * </pre>
 *
 * Markers do not get secondary identifiers here. They are merged strictly on primaryIdentifier.
 *
 * @author Sam Hokin
 */
public class SNPVCFFileConverter extends BioFileConverter {
	
    private static final Logger LOG = Logger.getLogger(SNPVCFFileConverter.class);

    Map<String,Item> chromMap;

    /**
     * Create a new SNPVCFFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public SNPVCFFileConverter(ItemWriter writer, Model model) {
        super(writer, model);
	chromMap = new HashMap<String,Item>();
    }

    /**
     * {@inheritDoc}
     * Read in the VCF file and store the SNP markers with chromosome (or supercontig) and position.
     */
    @Override
    public void process(Reader reader) throws Exception {

        // don't process README files
        if (getCurrentFile().getName().equals("README")) return;

        LOG.info("Processing VCF file "+getCurrentFile().getName()+"...");

        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {

	    SNPVCFRecord record = new SNPVCFRecord(line);

	    // will return null chrom if a comment
	    if (record.chrom==null) continue;
	    
	    // create this marker, length 1
	    Item marker = createItem("GeneticMarker");
	    marker.setAttribute("primaryIdentifier", record.id);
	    marker.setAttribute("type", "SNP");
	    marker.setAttribute("length", "1");

	    // get (or create and store) the chromosome (or supercontig) and references
	    boolean isSupercontig = (record.chrom.toLowerCase().contains("contig"));
	    Item chrom;
	    if (chromMap.containsKey(record.chrom)) {
		chrom = chromMap.get(record.chrom);
	    } else {
		if (isSupercontig) {
		    chrom = createItem("Supercontig");
		} else {
		    chrom = createItem("Chromosome");
		}
		chrom.setAttribute("primaryIdentifier", record.chrom);
		store(chrom);
		chromMap.put(record.chrom, chrom);
	    }
	    if (isSupercontig) {
		marker.setReference("supercontig", chrom);
	    } else {
		marker.setReference("chromosome", chrom);
	    }

	    // create the location on this chromosome/supercontig and references
	    Item location = createItem("Location");
	    location.setAttribute("start", String.valueOf(record.pos));
	    location.setAttribute("end", String.valueOf(record.pos));
	    location.setReference("locatedOn", chrom);
	    location.setReference("feature", marker);
	    if (isSupercontig) {
		marker.setReference("supercontigLocation", location);
	    } else {
		marker.setReference("chromosomeLocation", location);
	    }

	    // store the marker and location
	    store(marker);
	    store(location);

	}
        
        br.close();
        
    }


    /**
     * Store the items held in maps (there are none)
     */
    @Override
    public void close() throws Exception {
    }

}
