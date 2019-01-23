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
 * Store GWAS QTL/marker data from a SoyBase (or other) dump. See GWASFileRecord for individual record format.
 *
 * TaxonID	         3847
 * PrimaryIdentifier     KGK20170714.1 [required; must be first GWAS experiment datum]
 * PlatformName          SoySNP50k [optional]
 * PlatformDetails       Illumina Infinium BeadChip [optional]
 * NumberLociTested      52041 [optional]
 * NumberGermplasmTested 12116 [optional]
 * PubMedId              26224783 [optional]
 * DOI                   10.1534/g3.115.019000 [optional]
 *
 * @author Sam Hokin
 */
public class GWASFileConverter extends BioFileConverter {
	
    private static final Logger LOG = Logger.getLogger(GWASFileConverter.class);

    // store items in maps to avoid duplicates
    Map<String,Item> organismMap = new HashMap<String,Item>();
    Map<String,Item> chromosomeMap = new HashMap<String,Item>();
    Map<String,Item> markerMap = new HashMap<String,Item>();
    Map<String,Item> qtlMap = new HashMap<String,Item>();
    
    /**
     * Create a new GWASFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public GWASFileConverter(ItemWriter writer, Model model) {
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

        // persistent items
        Item organism = null;
        Item gwasExperiment = null;

        BufferedReader bufferedReader = new BufferedReader(reader);
	String line;
        while ((line=bufferedReader.readLine())!=null) {

            if (line.startsWith("#") || line.trim().length()==0) {
                continue;
            }
            
            String[] parts = line.split("\t");
            String key = parts[0];
            String value = parts[1];

            if (key.toLowerCase().equals("taxonid")) {
                String taxonId = value;
                if (organismMap.containsKey(taxonId)) {
                    organism = organismMap.get(taxonId);
                } else {
                    organism = createItem("Organism");
                    organism.setAttribute("taxonId", taxonId);
                    store(organism);
                    LOG.info("Created and stored Organism:"+taxonId);
                    organismMap.put(taxonId, organism);
                }

            } else if (key.toLowerCase().equals("primaryidentifier")) {
                // primaryIdentifier must be FIRST GWAS experiment record!
                gwasExperiment = createItem("GWASExperiment");
                gwasExperiment.setAttribute("primaryIdentifier", value);

            } else if (key.toLowerCase().equals("platformname")) {
                gwasExperiment.setAttribute("platformName", value);

            } else if (key.toLowerCase().equals("numberlocitested")) {
                gwasExperiment.setAttribute("numberLociTested", value);

            } else if (key.toLowerCase().equals("numbergermplasmtested")) {
                gwasExperiment.setAttribute("numberGermplasmTested", value);

            } else if (key.toLowerCase().equals("platformdetails")) {
                gwasExperiment.setAttribute("platformDetails", value);

            } else if (key.toLowerCase().equals("pmid")) {
                int pmid = Integer.parseInt(value);
                Item publication = createItem("Publication");
                publication.setAttribute("pubMedId", String.valueOf(pmid));
                store(publication);
                LOG.info("Stored publication PMID="+pmid);
                gwasExperiment.addToCollection("publications", publication);

            } else if (key.toLowerCase().equals("doi")) {
                String doi = value;
                Item publication = createItem("Publication");
                publication.setAttribute("doi", doi);
                store(publication);
                LOG.info("Stored publication DOI="+doi);
                gwasExperiment.addToCollection("publications", publication);

            } else {
                
                // check that we've got an organism - fatal exit if not
                if (organism==null) {
                    String errorMsg = "Organism has not been formed for GWAS record import in file "+getCurrentFile().getName();
                    LOG.error(errorMsg);
                    throw new RuntimeException(errorMsg);
                }
                
                // check that we've got a GWAS experiment - fatal exit if not
                if (gwasExperiment==null) {
                    String errorMsg = "GWAS experiment has not been created";
                    LOG.error(errorMsg);
                    throw new RuntimeException(errorMsg);
                }
                
                // process a record
                GWASFileRecord rec = new GWASFileRecord(line);
                
                //
                // Genetic Marker, Chromosome
                //
                Item marker = null;
                String markerKey = rec.locusName;
                if (markerMap.containsKey(markerKey)) {
                    marker = markerMap.get(markerKey);
                } else {
                    marker = createItem("GeneticMarker");
                    marker.setReference("organism", organism);
                    marker.setAttribute("primaryIdentifier", rec.locusName);
                    marker.setAttribute("type", rec.type);
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
                        chromosome.setReference("organism", organism);
                        chromosome.setAttribute("primaryIdentifier", rec.chromosome);
                        chromosomeMap.put(rec.chromosome, chromosome);
                        store(chromosome);
                        LOG.info("Created and stored chromosome/supercontig: "+rec.chromosome);
                    }
                    if (isSupercontig) {
                        marker.setReference("supercontig", chromosome);
                    } else {
                        marker.setReference("chromosome", chromosome);
                    }
                    // create and store the Location object, use + strand since it's not defined
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
                    markerMap.put(markerKey, marker);
                }

                //
                // QTL
                //
                Item qtl = null;
                String qtlKey = rec.gwasName;
                if (qtlMap.containsKey(qtlKey)) {
                    qtl = qtlMap.get(qtlKey);
                    LOG.info("Duplicate QTL:"+qtlKey+" with marker:"+markerKey);
                } else {
                    qtl = createItem("QTL");
                    qtl.setReference("organism", organism);
                    qtl.setAttribute("primaryIdentifier", rec.gwasName);
                    qtl.setAttribute("traitName", rec.gwasClass);
                    qtl.setAttribute("pValue", String.valueOf(rec.pValue));
                    qtl.setReference("gwasExperiment", gwasExperiment);
                    qtlMap.put(qtlKey, qtl);
                    LOG.info("Added QTL:"+qtlKey+" with marker:"+markerKey);
                }
                    
                // relate the QTL and marker
                marker.addToCollection("QTLs", qtl);
            }
        }

        // finally store this GWAS experiment
        store(gwasExperiment);
        
        bufferedReader.close();
    }

    /**
     * Store markers and QTLs.
     */
    @Override
    public void close() throws Exception {
        for (Item qtl : qtlMap.values()) store(qtl);
        for (Item marker : markerMap.values()) store(marker);
    }
}
