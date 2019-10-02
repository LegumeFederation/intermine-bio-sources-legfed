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
 * Store GWAS/marker data from a tab-delimited file.
 *
 * TaxonID                 3847
 * Strain                  Williams82
 * Name                    KGK20170714.1
 * PlatformName            SoySNP50k
 * PlatformDetails         Illumina Infinium BeadChip
 * NumberLociTested        52,041
 * NumberGermplasmTested   12,116
 * Assembly                Wm82.a2.v1 <-- CAREFUL about start/end!
 * DOI                     10.3835/plantgenome2015.04.0024
 * PMID                    123456
 *
 * Record fields loaded by GWASFileRecord:
 *
 * phenotype  ontology_identifier marker      p_value  chromosome start    end
 * Seed oil   SOY:0001668         ss715591649 1.12E-09 Gm05       41780982 41780982
 *
 * @author Sam Hokin
 */
public class GWASFileConverter extends BioFileConverter {
	
    private static final Logger LOG = Logger.getLogger(GWASFileConverter.class);

    // store items in maps to avoid duplicates
    Map<String,Item> organismMap = new HashMap<String,Item>();
    Map<String,Item> strainMap = new HashMap<String,Item>();
    Map<String,Item> chromosomeMap = new HashMap<String,Item>();
    Map<String,Item> phenotypeMap = new HashMap<String,Item>();
    Map<String,Item> ontologyTermMap = new HashMap<String,Item>();
    Map<String,Item> markerMap = new HashMap<String,Item>();
    
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
        Item strain = null;
        Item gwas = null;
        Item publication = null;

        BufferedReader bufferedReader = new BufferedReader(reader);
	String line;
        while ((line=bufferedReader.readLine())!=null) {

            String[] parts = line.split("\t");
            if (line.startsWith("#") || line.trim().length()==0 || parts.length<2) {
                continue;
            }
            
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
                    LOG.info("Stored organism:"+taxonId);
                    organismMap.put(taxonId, organism);
                }

            } else if (key.toLowerCase().equals("strain")) {
                String strainName = value;
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

            } else if (key.toLowerCase().equals("name")) {
                gwas = createItem("GWAS");
                gwas.setAttribute("primaryIdentifier", value);

            } else if (key.toLowerCase().equals("platformname")) {
                gwas.setAttribute("platformName", value);

            } else if (key.toLowerCase().equals("platformdetails")) {
                gwas.setAttribute("platformDetails", value);

            } else if (key.toLowerCase().equals("numberlocitested")) {
                gwas.setAttribute("numberLociTested", value);

            } else if (key.toLowerCase().equals("numbergermplasmtested")) {
                gwas.setAttribute("numberGermplasmTested", value);

            } else if (key.toLowerCase().equals("assembly")) {
                // do nothing with assembly

            } else if (key.toLowerCase().equals("pmid")) {
                int pmid = Integer.parseInt(value);
                publication = createItem("Publication");
                publication.setAttribute("pubMedId", String.valueOf(pmid));
                store(publication);
                LOG.info("Stored publication PMID="+pmid);
                gwas.addToCollection("publications", publication);

            } else if (key.toLowerCase().equals("doi")) {
                String doi = value;
                publication = createItem("Publication");
                publication.setAttribute("doi", doi);
                store(publication);
                LOG.info("Stored publication DOI="+doi);
                gwas.addToCollection("publications", publication);

            } else {
                
                // check that we've got an organism - fatal exit if not
                if (organism==null) {
                    String errorMsg = "Organism has not been set for GWAS record import in file "+getCurrentFile().getName();
                    LOG.error(errorMsg);
                    throw new RuntimeException(errorMsg);
                }
                
                // check that we've got a strain - fatal exit if not
                if (strain==null) {
                    String errorMsg = "Strain has not been set for GWAS record import in file "+getCurrentFile().getName();
                    LOG.error(errorMsg);
                    throw new RuntimeException(errorMsg);
                }
                
                // check that we've got a GWAS experiment - fatal exit if not
                if (gwas==null) {
                    String errorMsg = "GWAS experiment has not been created";
                    LOG.error(errorMsg);
                    throw new RuntimeException(errorMsg);
                }
                
                /////////////////////////////////
                // process a GWASResult record //
                /////////////////////////////////
                GWASFileRecord rec = new GWASFileRecord(line);

                Item marker = null;
                if (markerMap.containsKey(rec.marker)) {
                    marker = markerMap.get(rec.marker);
                } else {
                    marker = createItem("GeneticMarker");
                    marker.setReference("organism", organism);
                    marker.setReference("strain", strain);
                    marker.setAttribute("primaryIdentifier", rec.marker);
                    marker.setAttribute("type", rec.type);
                    // set the chromosome or supercontig reference
                    boolean isSupercontig = DatastoreUtils.isSupercontig(rec.chromosome);
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
                        store(chromosome);
                        LOG.info("Stored chromosome/supercontig: "+rec.chromosome);
                        chromosomeMap.put(rec.chromosome, chromosome);
                    }
                    if (isSupercontig) {
                        marker.setReference("supercontig", chromosome);
                    } else {
                        marker.setReference("chromosome", chromosome);
                    }

                    // COORDINATES MAY BE ON A DIFFERENT ASSEMBLY FROM THE MINE
                    // LOAD THE MARKER LOCATIONS FROM A DIFFERENT SOURCE
                    // // create and store the Location object, use + strand since it's not defined
                    // Item location = createItem("Location");
                    // location.setReference("locatedOn", chromosome);
                    // location.setReference("feature", marker);
                    // location.setAttribute("start", String.valueOf(rec.start));
                    // location.setAttribute("end", String.valueOf(rec.end));
                    // location.setAttribute("strand", String.valueOf(+1));
                    // store(location);
                    // // set the chromosomeLocation/supercontigLocation reference and store the marker
                    // if (isSupercontig) {
                    //     marker.setReference("supercontigLocation", location);
                    // } else {
                    //     marker.setReference("chromosomeLocation", location);
                    // }

                    store(marker);
                    markerMap.put(rec.marker, marker);
                }

                Item phenotype = null;
                if (phenotypeMap.containsKey(rec.phenotype)) {
                    phenotype = phenotypeMap.get(rec.phenotype);
                } else {
                    phenotype = createItem("Phenotype");
                    phenotype.setAttribute("primaryIdentifier", rec.phenotype);
                    phenotypeMap.put(rec.phenotype, phenotype);
                }
                if (publication!=null) phenotype.addToCollection("publications", publication);

                Item ontologyAnnotation = null;
                if (rec.ontologyIdentifier!=null) {
                    Item ontologyTerm = null;
                    if (ontologyTermMap.containsKey(rec.ontologyIdentifier)) {
                        ontologyTerm = ontologyTermMap.get(rec.ontologyIdentifier);
                    } else {
                        ontologyTerm = createItem("OntologyTerm");
                        ontologyTerm.setAttribute("identifier", rec.ontologyIdentifier);
                        store(ontologyTerm);
                        ontologyTermMap.put(rec.ontologyIdentifier, ontologyTerm);
                    }
                    ontologyAnnotation = createItem("OntologyAnnotation");
                    ontologyAnnotation.setReference("ontologyTerm", ontologyTerm);
                    ontologyAnnotation.setReference("subject", phenotype);
                    store(ontologyAnnotation);
                }

                Item gwasResult = createItem("GWASResult");
                if (rec.pvalue>0) gwasResult.setAttribute("pValue", String.valueOf(rec.pvalue));
                gwasResult.setReference("study", gwas);
                gwasResult.setReference("phenotype", phenotype);
                gwasResult.setReference("marker", marker);
                store(gwasResult);
            }
        }

        // finally store this GWAS experiment
        store(gwas);
        
        bufferedReader.close();
    }

    /**
     * Store some remaining Items.
     * @override
     */
    public void close() throws ObjectStoreException {
        store(phenotypeMap.values());
    }

}
