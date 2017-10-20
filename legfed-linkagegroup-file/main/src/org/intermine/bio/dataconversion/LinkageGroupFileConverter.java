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
 * Store Linkage Groups and their Genetic Maps from tab-delimited files. An optional length column may be added.
 *
 * LinkageGroup        GeneticMap      [Length(cM)]
 * GmComposite2003_D1a GmComposite2003 120.89
 *
 * @author Sam Hokin, NCGR
 */
public class LinkageGroupFileConverter extends BioFileConverter {
	
    private static final Logger LOG = Logger.getLogger(LinkageGroupFileConverter.class);

    // Store genetic maps in a map
    Map<String,Item> geneticMapMap = new HashMap<String,Item>();

    /**
     * Create a new LinkageGroupFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public LinkageGroupFileConverter(ItemWriter writer, Model model) {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     * We process the LinkageGroup-marker relationships by reading in from a tab-delimited file.
     */
    @Override
    public void process(Reader reader) throws Exception {

        // don't process README files
        if (getCurrentFile().getName().equals("README")) return;

        LOG.info("Processing file "+getCurrentFile().getName()+"...");

        BufferedReader linkageGroupReader = new BufferedReader(reader);
	String line;
        while ((line=linkageGroupReader.readLine())!=null) {

            if (!line.startsWith("#")) {

                String[] parts = line.split("\t");
                // req'd
                String lgID = parts[0];
                String gmID = parts[1];
                // optional
                double length = 0.0;
                if (parts.length>2) length = Double.parseDouble(parts[2]);

                Item lg = createItem("LinkageGroup");
                lg.setAttribute("primaryIdentifier", lgID);
                if (length>0.0) lg.setAttribute("length", String.valueOf(length));

                Item gm = null;
                if (geneticMapMap.containsKey(gmID)) {
                    gm = geneticMapMap.get(gmID);
                } else {
                    // create and store a new genetic map
                    gm = createItem("GeneticMap");
                    gm.setAttribute("primaryIdentifier", gmID);
                    store(gm);
                    geneticMapMap.put(gmID, gm);
                }

                lg.setReference("geneticMap", gm);
                store(lg);
                
                gm.addToCollection("linkageGroups", lg);

            }

        }
        
        linkageGroupReader.close();

    }

    /**
     * Do nothing, we're storing as we build above.
     */
    @Override
    public void close() throws Exception {

    }
    
}
