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
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.IOException;
import java.io.UnsupportedEncodingException;

import java.net.URL;
import java.net.URLConnection;
import java.net.URLEncoder;

import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.Document;
import org.w3c.dom.NodeList;
import org.w3c.dom.Node;
import org.w3c.dom.Element;

import org.xml.sax.SAXException;


/**
 * Performs a search of NCBI PubMed on journal, year and authors.
 *
 * Prototypical search URL:
 * http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=molecular+breeding[journal]+AND+2012[pdat]+AND+Blair[author]+AND+Galeano[author]
 *
 * @author Sam Hokin
 */
public class PubMedSearch {

    public static String CHARSET = "UTF-8";

    // append query to this URL
    public static String PUBMED_URL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed";
    
    /**
     * Main class for convenient command-line usage
     */
    public static void main(String[] args) {

        // validate args
        if (args.length<3) {
            System.out.println("Usage: PubMedSearch <journal> <year> <author1> [author2] ...");
            System.exit(1);
        }

        // get args
        String journal = args[0];
        int year = Integer.parseInt(args[1]);
        String authors[] = new String[args.length-2];
        for (int i=0; i<authors.length; i++) {
            authors[i] = args[2+i];
        }

        try {
            int pubMedId = getPubMedId(journal, year, authors);
            System.out.println("PubMedID="+pubMedId);
        } catch (Exception ex) {
            System.err.println(ex.getMessage());
        }

    }

    /**
     * Query PubMed for the journal in question, return 0 if none or more than one found, PubMedId if only one found.
     */
    public static int getPubMedId(String journal, int year, String[] authors) throws IOException, UnsupportedEncodingException, ParserConfigurationException, SAXException {

        // form URL
        String url = PUBMED_URL;
        url += "&term=";
        url += URLEncoder.encode(journal,CHARSET)+"[journal]";
        url += "+AND+";
        url += year+"[pdat]";
        for (int i=0; i<authors.length; i++) {
            url += "+AND+";
            url += URLEncoder.encode(authors[i],CHARSET)+"[author]";
        }

        DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
        DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
        Document doc = dBuilder.parse(url);

        //optional, but recommended
        //read this - http://stackoverflow.com/questions/13786607/normalization-in-dom-parsing-with-java-how-does-it-work
        doc.getDocumentElement().normalize();

        int count = Integer.parseInt(doc.getElementsByTagName("Count").item(0).getTextContent());

        if (count==1) {
            int pubMedId = Integer.parseInt(doc.getElementsByTagName("Id").item(0).getTextContent());
            return pubMedId;
        } else {
            return 0;
        }

    }

}
