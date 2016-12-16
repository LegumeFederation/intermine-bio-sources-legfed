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

    String taxonId;
    
    private static final Logger LOG = Logger.getLogger(ExpressionFileConverter.class);

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
            Item publication = PublicationTools.storePublicationFromPMID(this, Integer.parseInt(pubMedId));
            if (publication!=null) source.setReference("publication", publication);
        }
        store(source);

        // 2. expression unit
        String[] unitLine = (String[]) tsvIter.next();
        String unit = unitLine[0];

        // 3. number of samples
        String[] numSamplesLine = (String[]) tsvIter.next();
        int numSamples = Integer.parseInt(numSamplesLine[0]);
        
        // 4. Sample number, name and description in same order as the expression columns
        Item[] samples = new Item[numSamples];
        for (int i=0; i<numSamples; i++) {
            String[] sampleLine = (String[]) tsvIter.next();
            samples[i] = createItem("ExpressionSample");
            samples[i].setAttribute("num", sampleLine[0]);
            samples[i].setAttribute("primaryIdentifier", sampleLine[1]);
            samples[i].setAttribute("description", sampleLine[2]);
            samples[i].setReference("source", source);
            store(samples[i]);
        }

        // 5. expression data - rely on mRNA merge on primaryIdentifier
        // If mRNAId does not end in .#, append ".1" to convert from gene to mRNA accession
        while (tsvIter.hasNext()) {
            String[] line = (String[]) tsvIter.next();
            String mRNAId = line[0];
            if (mRNAId.charAt(mRNAId.length()-2)!='.') mRNAId += ".1";
            Item mRNA = createItem("MRNA");
            mRNA.setAttribute("primaryIdentifier", mRNAId);
            store(mRNA);
            // load the expression values
            int offset = 1;
            for (int i=0; i<numSamples; i++) {
                Item expressionValue = createItem("ExpressionValue");
                expressionValue.setAttribute("value", line[i+offset]);
                expressionValue.setAttribute("unit", unit);
                expressionValue.setReference("expressionSample", samples[i]);
                expressionValue.setReference("mRNA", mRNA);
                store(expressionValue);
            }
        }

    }

}
