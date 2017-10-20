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

import org.ncgr.intermine.PublicationTools;
import org.ncgr.pubmed.PubMedSummary;

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.xml.full.Item;

/**
 * Store QTL data and relationships from tab-delimited files.
 *
 * QTLName Parent_1 Parent_2 TraitName Journal Year Volume Page Title PMID
 *
 * Parents are munged into a mapping population name, e.g. Sanzi_x_Vita7 for parents Sanzi and Vita7.
 *
 * @author Sam Hokin, NCGR
 */
public class QTLFileConverter extends BioFileConverter {
	
    private static final Logger LOG = Logger.getLogger(QTLFileConverter.class);

    // Store QTLs in a map since there may be multiple lines (different pubs)
    Map<String,Item> qtlMap = new HashMap<String,Item>();

    // Store mapping populations in a map, presumably far fewer than QTLs
    Map<String,Item> mappingPopulationMap = new HashMap<String,Item>();

    // Store publications in a map
    Map<String,Item> publicationMap = new HashMap<String,Item>();
    Map<String,Item> authorMap = new HashMap<String,Item>();
    
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
        if (getCurrentFile().getName().equals("README")) return;

        LOG.info("Processing file "+getCurrentFile().getName()+"...");

        BufferedReader qtlReader = new BufferedReader(reader);
	String line;
        while ((line=qtlReader.readLine())!=null) {

            QTLRecord rec = new QTLRecord(line);
            if (rec.qtlName!=null) {
                    
                // find the QTL in the map, or add it with qtlName=primaryIdentifier
                Item qtl = null;
                if (qtlMap.containsKey(rec.qtlName)) {
                    qtl = qtlMap.get(rec.qtlName);
                } else {
                    qtl = createItem("QTL");
                    qtl.setAttribute("primaryIdentifier", rec.qtlName);
                    qtl.setAttribute("traitName", rec.traitName);
                    qtlMap.put(rec.qtlName, qtl);
                }

                // mapping population, if parents present
                if (rec.parent1.length()>0 && rec.parent2.length()>0) {
                    String mappingPopulationID = rec.parent1+"_x_"+rec.parent2;
                    Item mappingPopulation = null;
                    if (mappingPopulationMap.containsKey(mappingPopulationID)) {
                        mappingPopulation = mappingPopulationMap.get(mappingPopulationID);
                    } else {
                        mappingPopulation = createItem("MappingPopulation");
                        mappingPopulation.setAttribute("primaryIdentifier", mappingPopulationID);
                        mappingPopulationMap.put(mappingPopulationID, mappingPopulation);
                    }
                    qtl.addToCollection("mappingPopulations", mappingPopulation);
                }

                // publication, if present
                if (rec.pubTitle.length()>0) {
                    Item publication;
                    if (publicationMap.containsKey(rec.pubTitle)) {
                        publication = publicationMap.get(rec.pubTitle);
                    } else {
                        PubMedSummary pms;
                        if (rec.pubPMID==0) {
                            // search on title
                            pms = new PubMedSummary(rec.pubTitle);
                        } else {
                            // search on PMID
                            pms = new PubMedSummary(rec.pubPMID);
                        }
                        if (pms.id>0) {
                            publication = createItem("Publication");
                            publication.setAttribute("title", pms.title);
                            publication.setAttribute("pubMedId", String.valueOf(pms.id));
                            if (pms.doi!=null && pms.doi.length()>0) publication.setAttribute("doi", pms.doi);
                            if (pms.issue!=null && pms.issue.length()>0) publication.setAttribute("issue", pms.issue);
                            if (pms.pages!=null && pms.pages.length()>0) publication.setAttribute("pages", pms.pages);
                            // parse year, month from PubDate
                            if (pms.pubDate!=null && pms.pubDate.length()>0) {
                                String[] dateBits = pms.pubDate.split(" ");
                                publication.setAttribute("year",dateBits[0]);
                                if (dateBits.length>1) publication.setAttribute("month",dateBits[1]);
                            }
                            if (pms.volume.length()>0) publication.setAttribute("volume", pms.volume);
                            if (pms.fullJournalName.length()>0) publication.setAttribute("journal", pms.fullJournalName);
                            // authors collection
                            if (pms.authorList!=null && pms.authorList.size()>0) {
                                boolean firstAuthor = true;
                                for (String author : pms.authorList) {
                                    if (firstAuthor) {
                                        publication.setAttribute("firstAuthor", author);
                                        firstAuthor = false;
                                    }
                                    Item authorItem;
                                    if (authorMap.containsKey(author)) {
                                        authorItem = authorMap.get(author);
                                    } else {
                                        authorItem = createItem("Author");
                                        authorItem.setAttribute("name", author);
                                        authorMap.put(author, authorItem);
                                    }
                                    publication.addToCollection("authors", authorItem);
                                }
                            }
                            // add to map
                            publicationMap.put(rec.pubTitle, publication);
                            // add to QTL collection
                            qtl.addToCollection("publications", publication);
                        }
                    }
                }
                
            }
            
        }
        
        qtlReader.close();

    }

    /**
     * Store the items we've collected in maps and sets
     */
    @Override
    public void close() throws Exception {

        LOG.info("Storing "+qtlMap.size()+" QTL items...");
        store(qtlMap.values());

        LOG.info("Storing "+mappingPopulationMap.size()+" MappingPopulation items...");
        store(mappingPopulationMap.values());

        LOG.info("Storing "+publicationMap.size()+" Publication items...");
        store(publicationMap.values());

        LOG.info("Storing "+authorMap.size()+" Author items...");
        store(authorMap.values());

    }
    
}
