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
import java.util.List;
import java.util.ArrayList;

import org.apache.log4j.Logger;

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

/**
 * Store QTL data and relationships from tab-delimited files organized by publication.
 * Mapping population name may be given, e.g. Sanzi_x_Vita7, OR genotyping study, if appropriate.
 *
 * TaxonID     3817
 * PMID        123456 [optional]
 * DOI         10.1101/202044 [optional]
 * Journal     Intl. J. Genetics [optional]
 * Year        2007 [optional]
 * Volume      23 [optional]
 * Pages       345-357 [optional]
 * Title       A cowpea genetic map for the ages [optional]
 * MappingPopulation Sanzi_x_Vita7 [optional, can have more than one; parent1_x_parent2 results in strains being created and associated as parents]
 * GenotypingStudy   MAGIC-2007 [optional]
 * Description A description of the purpose of this mapping population/genotyping study and why the parents were chosen. [optional]
 * #QTLName          Phenotype
 * Seed weight 3-2   Seed weight
 *
 * @author Sam Hokin, NCGR
 */
public class QTLFileConverter extends BioFileConverter {
	
    private static final Logger LOG = Logger.getLogger(QTLFileConverter.class);

    // maps for Items that occur in multiple files
    Map<String,Item> organismMap = new HashMap<>();
    Map<String,Item> strainMap = new HashMap<>();
    Map<String,Item> phenotypeMap = new HashMap<>();
    Map<String,Item> qtlMap = new HashMap<>();
    Map<String,Item> mappingPopulationMap = new HashMap<>();
    Map<String,Item> genotypingStudyMap = new HashMap<>();

    /**
     * Create a new QTLFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public QTLFileConverter(ItemWriter writer, Model model) {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     * We process the QTL-marker relationships by reading in from a tab-delimited file.
     */
    @Override
    public void process(Reader reader) throws Exception {

        // don't process README files
        if (getCurrentFile().getName().contains("README")) return;

        LOG.info("Processing file "+getCurrentFile().getName()+"...");

        // these are loaded to populate a publication by title rather than DOI or PMID.
        String title = null;
        String journal = null;
        int year = 0;
        String volume = null;
        String pages = null;

        // the Items resulting from the header values
        Item organism = null;
        Item publication = null;
        List<Item> mappingPopulations = new ArrayList<Item>();
	Item genotypingStudy = null;
        
        BufferedReader br = new BufferedReader(reader);
	String line = null;
        while ((line=br.readLine())!=null) {

	    String[] parts = line.split("\t");
            if (line.trim().length()==0 || line.startsWith("#") || parts.length!=2) {
		continue;
	    }
		
            String key = parts[0].trim();
            String value = parts[1].trim();
		
            // header data
            if (key.toLowerCase().equals("taxonid")) {
                String taxonId = value;
                if (organismMap.containsKey(taxonId)) {
                    organism = organismMap.get(taxonId);
                } else {
                    organism = createItem("Organism");
                    organism.setAttribute("taxonId", taxonId);
                    store(organism);
                    organismMap.put(taxonId, organism);
                    LOG.info("Stored organism: "+taxonId);
                }
                
            } else if (key.toLowerCase().equals("pmid")) {
                int pmid = Integer.parseInt(value);
                publication = createItem("Publication");
                publication.setAttribute("pubMedId", String.valueOf(pmid));
                store(publication);
                LOG.info("Stored publication PMID="+pmid);
                
            } else if (key.toLowerCase().equals("doi")) {
                String doi = value;
                publication = createItem("Publication");
                publication.setAttribute("doi", doi);
                store(publication);
                LOG.info("Stored publication DOI="+doi);
                
            } else if (key.toLowerCase().equals("title")) {
                title = value;
                
            } else if (key.toLowerCase().equals("journal")) {
                journal = value;
                
            } else if (key.toLowerCase().equals("year")) {
                year = Integer.parseInt(value);
                
            } else if (key.toLowerCase().equals("volume")) {
                volume = value;
                
            } else if (key.toLowerCase().equals("pages")) {
                pages = value;
                
            } else if (key.toLowerCase().equals("mappingpopulation")) {
                
                String mappingPopulationName = value;
                if (mappingPopulationMap.containsKey(mappingPopulationName)) {
                    Item mappingPopulation = mappingPopulationMap.get(mappingPopulationName);
                    mappingPopulations.add(mappingPopulation);
                } else {
                    Item mappingPopulation = createItem("MappingPopulation");
                    mappingPopulation.setAttribute("primaryIdentifier", mappingPopulationName);
                    mappingPopulation.setReference("organism", organism);
                    mappingPopulationMap.put(mappingPopulationName, mappingPopulation);
                    mappingPopulations.add(mappingPopulation);
                    LOG.info("Storing mapping population:"+mappingPopulationName);
                    if (mappingPopulationName.contains("_x_")) {
                        // add parents
                        String[] strainNames = mappingPopulationName.split("_x_");
                        for (String strainName : strainNames) {
                            if (strainMap.containsKey(strainName)) {
                                Item parent = strainMap.get(strainName);
                                mappingPopulation.addToCollection("parents", parent);
                            } else {
                                Item parent = createItem("Strain");
                                parent.setAttribute("primaryIdentifier", strainName);
                                parent.setReference("organism", organism);
                                store(parent);
                                strainMap.put(strainName, parent);
                                LOG.info("Stored parent: "+strainName);
                                mappingPopulation.addToCollection("parents", parent);
                            }
                        }
                    }
                }
                
            } else if (key.toLowerCase().equals("genotypingstudy")) {

                String genotypingStudyName = value;
                if (genotypingStudyMap.containsKey(genotypingStudyName)) {
                    genotypingStudy = genotypingStudyMap.get(genotypingStudyName);
                } else {
                    genotypingStudy = createItem("GenotypingStudy");
                    genotypingStudy.setAttribute("primaryIdentifier", genotypingStudyName);
                    genotypingStudy.setReference("organism", organism);
                    genotypingStudyMap.put(genotypingStudyName, genotypingStudy);
                    LOG.info("Storing genotyping study:"+genotypingStudyName);
                }

            } else if (key.toLowerCase().equals("description")) {
                
                String description = value;
                for (Item mappingPopulation : mappingPopulations) {
                    mappingPopulation.setAttribute("description", description);
                }
                LOG.info("Set description on mapping populations:"+description);
                
            } else {
                
                // create and store the publication if we've collected at least a title
                if (publication==null && title!=null && title.length()>0) {
                    publication = createItem("Publication");
                    publication.setAttribute("title", title);
                    if (journal!=null && journal.length()>0) publication.setAttribute("journal", journal);
                    if (year!=0) publication.setAttribute("year", String.valueOf(year));
                    if (volume!=null && volume.length()>0) publication.setAttribute("volume", volume);
                    if (pages!=null && pages.length()>0) publication.setAttribute("pages", pages);
                    store(publication);
                    LOG.info("Storing publication: "+title);
                }

                // associate this publication with the mapping populations or genotyping study
                if (publication!=null) {
                    for (Item mappingPopulation : mappingPopulations) {
                        mappingPopulation.addToCollection("publications", publication);
                    }
                    if (genotypingStudy!=null) {
                        genotypingStudy.addToCollection("publications", publication);
                    }
                }

                // load a QTL and phenotype
                String qtlName = key;
                String phenotypeName = value;

                Item qtl = null;
                if (qtlMap.containsKey(qtlName)) {
                    qtl = qtlMap.get(qtlName);
                } else {
                    qtl = createItem("QTL");
                    qtl.setAttribute("primaryIdentifier", qtlName);
                    qtl.setReference("organism", organism);
                    qtlMap.put(qtlName, qtl);
                }

                Item phenotype = null;
                if (phenotypeName.length()>0) {
                    if (phenotypeMap.containsKey(phenotypeName)) {
                        phenotype = phenotypeMap.get(phenotypeName);
                    } else {
                        phenotype = createItem("Phenotype");
                        phenotype.setAttribute("primaryIdentifier", phenotypeName);
                        phenotypeMap.put(phenotypeName, phenotype);
                    }
                    qtl.setReference("phenotype", phenotype);
                }

                // add pub and mapping pops to this QTL's collections
                if (publication!=null) qtl.addToCollection("publications", publication);
                for (Item mappingPopulation : mappingPopulations) {
                    qtl.addToCollection("mappingPopulations", mappingPopulation);
                }
                if (genotypingStudy!=null) {
                    qtl.addToCollection("genotypingStudies", genotypingStudy);
                }
	    }
	}
    }

    /**
     * Store the items we've collected in maps across all the files
     */
    @Override
    public void close() throws ObjectStoreException {
        store(mappingPopulationMap.values());
	store(genotypingStudyMap.values());
        store(qtlMap.values());
        store(phenotypeMap.values());
    }
}
