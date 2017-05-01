package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2002-2016 FlyMine, Legume Federation
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */

import java.io.File;
import java.io.IOException;
import java.io.Reader;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.HashSet;
import java.util.Set;

import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.apache.tools.ant.BuildException;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.util.FormattedTextParser;
import org.intermine.xml.full.Item;

import org.ncgr.intermine.PublicationTools;
import org.ncgr.pubmed.PubMedSummary;

/**
 * DataConverter to create items from tab-delimited expression files, like the one downloaded from PvGEA.
 *
 * @author Sam Hokin
 */
public class ExpressionFileConverter extends BioFileConverter {
    
    private static final Logger LOG = Logger.getLogger(ExpressionFileConverter.class);

    Map<String,Item> transcriptMap = new HashMap<String,Item>();

    Set<Item> authorSet = new HashSet<Item>();
    Set<Item> publicationSet = new HashSet<Item>();
    Set<Item> sourceSet = new HashSet<Item>();
    Set<Item> sampleSet = new HashSet<Item>();
    Set<Item> valueSet = new HashSet<Item>();
    
    /**
     * Constructor.
     *
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     * @throws ObjectStoreException os
     */
    public ExpressionFileConverter(ItemWriter writer, Model model) throws ObjectStoreException {
        super(writer, model);
    }

    /**
     * Called for each file found by ant.
     *
     * {@inheritDoc}
     */
    public void process(Reader reader) throws Exception {

        if (getCurrentFile().getName().equals("README")) return;
        LOG.info("Processing expression file:"+getCurrentFile().getName()+"...");
        
        // parse the tab-delimited file
        Iterator<?> tsvIter;
        try {
            tsvIter = FormattedTextParser.parseTabDelimitedReader(reader);
        } catch (Exception ex) {
            throw new RuntimeException("Cannot parse file:"+getCurrentFile().getName(), ex);
        }

        // 1. source name, PMID, geo, bioProject, sra
        String[] sourceLine = (String[]) tsvIter.next();
        String sourceName = sourceLine[0];
        String pmid = sourceLine[1];
        String geo = sourceLine[2];
        String bioProject = sourceLine[3];
        String sra = sourceLine[4];

        LOG.info("Loading expression source:"+sourceName);
        Item source = createItem("ExpressionSource");
        source.setAttribute("primaryIdentifier", sourceName);
        if (geo.length()>0) source.setAttribute("geo", geo);
        if (bioProject.length()>0) source.setAttribute("bioProject", bioProject);
        if (sra.length()>0) source.setAttribute("sra", sra);
        // create and store the publication if it exists; this requires Internet access
        if (pmid!=null && pmid.trim().length()>0) {
            try {
                int id = Integer.parseInt(pmid);
                PubMedSummary pms = new PubMedSummary(id);
                Item publication = createItem("Publication");
                publication.setAttribute("title", pms.title);
                if (pms.doi.length()>0) publication.setAttribute("doi", pms.doi);
                if (pms.issue.length()>0) publication.setAttribute("issue", pms.issue);
                publication.setAttribute("pubMedId", String.valueOf(pms.id));
                publication.setAttribute("pages", pms.pages);
                // parse year, month from PubDate
                String[] dateBits = pms.pubDate.split(" ");
                publication.setAttribute("year",dateBits[0]);
                publication.setAttribute("month",dateBits[1]);
                publication.setAttribute("volume", pms.volume);
                publication.setAttribute("journal", pms.fullJournalName);
                // authors collection
                boolean firstAuthor = true;
                for (String author : pms.authorList) {
                    if (firstAuthor) {
                        publication.setAttribute("firstAuthor", author);
                        firstAuthor = false;
                    }
                    Item authorItem = createItem("Author");
                    authorItem.setAttribute("name", author);
                    authorSet.add(authorItem);
                    publication.addToCollection("authors", authorItem);
                }
                // store it and add reference to ExpressionSource
                publicationSet.add(publication);
                source.setReference("publication", publication);
            } catch (Exception ex) {
                throw new RuntimeException("Cannot create publication with PMID="+pmid, ex);
            }
        }
        sourceSet.add(source);

        // 2. URL
        String[] urlLine = (String[]) tsvIter.next();
        String url = urlLine[0];
        source.setAttribute("url", url);

        // 3. Description
        String[] descriptionLine = (String[]) tsvIter.next();
        String description = descriptionLine[0];
        source.setAttribute("description", description);

        // 4. expression unit
        String[] unitLine = (String[]) tsvIter.next();
        String unit = unitLine[0];
        source.setAttribute("unit", unit);

        // 5. number of samples
        String[] numSamplesLine = (String[]) tsvIter.next();
        int numSamples = Integer.parseInt(numSamplesLine[0]);
        
        // 6. Sample number, name and description in same order as the expression columns
        Item[] samples = new Item[numSamples];
        for (int i=0; i<numSamples; i++) {
            String[] sampleLine = (String[]) tsvIter.next();
            samples[i] = createItem("ExpressionSample");
            samples[i].setAttribute("num", sampleLine[0]);
            samples[i].setAttribute("primaryIdentifier", sampleLine[1]);
            samples[i].setAttribute("description", sampleLine[2]);
            samples[i].setReference("source", source);
            sampleSet.add(samples[i]);
        }

        // 7. expression data
        // If transcriptId does not end in .#, append ".n" to convert from gene to transcript accession, assuming .1, .2, .3 in order of rows in expression file
        while (tsvIter.hasNext()) {
            String[] line = (String[]) tsvIter.next();
            String transcriptId = line[0];
            if (transcriptId.charAt(transcriptId.length()-2)!='.') {
                // it's a gene, not transcript ID, so increment .suffix until we've got a new transcriptId
                int n = 1;
                while (transcriptMap.containsKey(transcriptId+"."+n)) {
                    n++;
                }
                transcriptId += "."+n;
            }
            Item transcript = null;
            if (transcriptMap.containsKey(transcriptId)) {
                // transcript exists from previous expression file
                transcript = transcriptMap.get(transcriptId);
            } else {
                // new transcript
                transcript  = createItem("Transcript");
            }
            transcript.setAttribute("primaryIdentifier", transcriptId);
            transcriptMap.put(transcriptId, transcript);
            // load the expression values
            for (int i=0; i<numSamples; i++) {
                Item value = createItem("ExpressionValue");
                value.setAttribute("value", line[i+1]);
                value.setReference("sample", samples[i]);
                value.setReference("transcript", transcript);
                valueSet.add(value);
            }
        }

    }

    /**
     * Store the items we've held in sets or maps.
     */
    @Override
    public void close() throws Exception {

        LOG.info("Storing "+transcriptMap.size()+" transcript items...");
        for (Item transcript : transcriptMap.values()) store(transcript);

        LOG.info("Storing "+authorSet.size()+" author items...");
        for (Item author : authorSet) store(author);

        LOG.info("Storing "+publicationSet.size()+" publication items...");
        for (Item publication : publicationSet) store(publication);
        
        LOG.info("Storing "+sourceSet.size()+" ExpressionSource items...");
        for (Item source : sourceSet) store(source);
        
        LOG.info("Storing "+sampleSet.size()+" ExpressionSample items...");
        for (Item sample : sampleSet) store(sample);

        LOG.info("Storing "+valueSet.size()+" ExpressionValue items...");
        for (Item value : valueSet) store(value);

    }

}
