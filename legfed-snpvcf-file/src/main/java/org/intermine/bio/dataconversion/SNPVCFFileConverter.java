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
 * #CHROM  POS	ID	REF	ALT	QUAL	FILTER	INFO
 * Vu01	  74363	2_37329	A	G	999	.	DP=378
 *
 * Markers do not get secondary identifiers here; they are merged on primary identifier.
 *
 * NOTE: it is assumed that only chromosomes, not supercontigs, are listed in the VCF file.
 *
 * @author Sam Hokin
 */
public class SNPVCFFileConverter extends BioFileConverter {
	
    private static final Logger LOG = Logger.getLogger(SNPVCFFileConverter.class);

    Map<String,Item> chromosomeMap;

    /**
     * Create a new SNPVCFFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public SNPVCFFileConverter(ItemWriter writer, Model model) {
        super(writer, model);
	chromosomeMap = new HashMap<String,Item>();
    }

    /**
     * {@inheritDoc}
     * Read in the VCF file and store the SNP markers with chromosome and position.
     */
    @Override
    public void process(Reader reader) throws Exception {

        // don't process README files
        if (getCurrentFile().getName().contains("README")) return;

        LOG.info("Processing VCF file "+getCurrentFile().getName()+"...");
        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {

            if (!line.startsWith("#") && line.trim().length()>0) {

		// parse the line
		SNPVCFRecord rec = new SNPVCFRecord(line);
		
		// create this marker, length 1
		Item marker = createItem("GeneticMarker");
		marker.setAttribute("primaryIdentifier", rec.id);
		marker.setAttribute("type", "SNP");
		marker.setAttribute("length", "1");
		
		// set the chromosome reference
		Item chromosome;
		if (chromosomeMap.containsKey(rec.chromosome)) {
		    chromosome = chromosomeMap.get(rec.chromosome);
		} else {
		    chromosome = createItem("Chromosome");
		    chromosome.setAttribute("primaryIdentifier", rec.chromosome);
		    store(chromosome);
		    chromosomeMap.put(rec.chromosome, chromosome);
		    LOG.info("Created and stored chromosome: "+rec.chromosome);
		}
		marker.setReference("chromosome", chromosome);
		
		// create and store the location on this chromosome
		Item location = createItem("Location");
		location.setReference("locatedOn", chromosome);
		location.setReference("feature", marker);
		location.setAttribute("start", String.valueOf(rec.pos));
		location.setAttribute("end", String.valueOf(rec.pos));
		store(location);

		// set the chromosomeLocation reference and store the marker
		marker.setReference("chromosomeLocation", location);
		
		// store the marker
		store(marker);

            }

        }
        
        br.close();
        
    }
}
