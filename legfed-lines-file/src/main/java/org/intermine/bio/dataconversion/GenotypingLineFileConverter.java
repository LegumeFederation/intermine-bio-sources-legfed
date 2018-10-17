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

import java.util.List;
import java.util.LinkedList;

import org.apache.log4j.Logger;

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.xml.full.Item;

/**
 * Loads genotyping lines from a tab-delimited file containing:
 *
 * primaryIdentifier  origin  comment
 *
 * @author Sam Hokin
 */
public class GenotypingLineFileConverter extends BioFileConverter {
	
    private static final Logger LOG = Logger.getLogger(GenotypingLineFileConverter.class);

    // store lines so we only load them once
    List<String> lines = new LinkedList<>();
    
    /**
     * Create a new GenotypingLineFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public GenotypingLineFileConverter(ItemWriter writer, Model model) {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     * Process the lines by reading in from a tab-delimited file.
     */
    @Override
    public void process(Reader reader) throws Exception {

        // don't process README files
        if (getCurrentFile().getName().contains("README")) return;

        LOG.info("Processing file "+getCurrentFile().getName()+"...");

        BufferedReader br = new BufferedReader(reader);
	String line;
        while ((line=br.readLine())!=null) {

            if (line.startsWith("#") || line.trim().length()==0) {

                // do nothing, comment
                
            } else {

                String[] parts = line.split("\t");
                String primaryIdentifier = parts[0];
                String origin = null;
                if (parts.length>1) origin = parts[1];
                String comment = null;
                if (parts.length>2) comment = parts[2];

                // only load fresh records
                if (!lines.contains(primaryIdentifier)) {
                    Item genotypingLine = createItem("GenotypingLine");
                    genotypingLine.setAttribute("primaryIdentifier", primaryIdentifier);
                    if (origin!=null) genotypingLine.setAttribute("origin", origin);
                    if (comment!=null) genotypingLine.setAttribute("comment", comment);
                    store(genotypingLine);
                    lines.add(primaryIdentifier);
                }

            }
            
        }
        
        br.close();

    }
    
}
