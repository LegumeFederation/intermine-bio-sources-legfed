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
import java.io.File;
import java.io.FileReader;
import java.io.Reader;
import java.io.IOException;

import java.util.Arrays;
import java.util.Map;
import java.util.HashMap;
import java.util.Set;
import java.util.HashSet;

import org.apache.log4j.Logger;

import org.intermine.dataconversion.FileConverter;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.metadata.StringUtil;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

/**
 * Loads data from LIS datastore files, using the file names to detect what sort of files they are, and,
 * sometimes, what the organism is.
 *
 * @author Sam Hokin
 */
public class DatastoreFileConverter extends FileConverter {
	
    private static final Logger LOG = Logger.getLogger(DatastoreFileConverter.class);

    // store everything in maps
    Map<String,Item> organisms = new HashMap<>();
    Map<String,Item> strains = new HashMap<>();
    Map<String,Item> ontologyTerms = new HashMap<>();
    Map<String,Item> dataSets = new HashMap<>();
    Map<String,Item> geneFamilies = new HashMap<>();
    Map<String,Item> proteinDomains = new HashMap<>();
    
    Set<String> ontologyAnnotations = new HashSet<>(); // concatenation of subject and term identifiers to avoid dupes

    // these are in maps of maps since we can have the same secondaryIdentifier from multiple assembly/annotation versions
    // these are merged on secondaryIdentifier,assemblyVersion,annotationVersion
    Map<String,Map<String,Item>> genes = new HashMap<String,Map<String,Item>>();
    Map<String,Map<String,Item>> proteins = new HashMap<String,Map<String,Item>>();
    Map<String,Map<String,Item>> mRNAs = new HashMap<String,Map<String,Item>>();

    // Taxon IDs of organisms of interest, set in project.xml
    Set<String> taxonIds;

    // ontologies created in constructor for random use
    Item geneOntology;
    Item pfamOntology;
    Item pantherOntology;
    Item kogOntology;
    Item ecOntology;
    Item koOntology;

    // there is only one dataSource, set in project.xml
    Item dataSource;

    // optional datasource and dataset attributes set in project.xml
    String dataSourceName;
    String dataSourceUrl;
    String dataSourceDescription;
    String dataSetUrl;
    String dataSetDescription;

    // map gensp to taxon Ids for associating dir/file names with organisms
    private static Map<String,String> genspTaxonIds = new HashMap<String,String>() {{
            put("aradu","130453");
            put("arahy","3818");
            put("araip","130454");
            put("cajca","3821");
            put("cicar","3827");
            put("glyma","3847");
            put("lotja","34305");
            put("lupan","3871");
            put("medtr","3880");
            put("phavu","3885");
            put("tripr","57577 ");
            put("vigan","3914");
            put("vigra","157791");
            put("vigun","3920");
        }};

    // map Genus_species to taxon Ids for associating dir/file names with organisms
    private static Map<String,String> genusSpeciesTaxonIds = new HashMap<String,String>() {{
            put("Arachis_duranensis","130453");
            put("Arachis_hypogaea","3818");
            put("Arachis_ipaensis","130454");
            put("Cajanus_cajan","3821");
            put("Cicer_arietinum","3827");
            put("Glycine_max","3847");
            put("Lotus_japonicus","34305");
            put("Lupinus angustifolius","3871");
            put("Medicago_truncatula","3880");
            put("Phaseolus_vulgaris","3885");
            put("Trifolium_pratense","57577");
            put("Vigna_angularis","3914");
            put("Vigna_radiata","157791");
            put("Vigna_unguiculata","3920");
        }};

    /**
     * Set the list of taxonIds that should be imported in cases where multiple organism data are present (like gene families).
     *
     * @param taxonIds a space-separated list of taxonIds
     */
    public void setOrganisms(String taxonIdString) {
        this.taxonIds = new HashSet<String>(Arrays.asList(StringUtil.split(taxonIdString, " ")));
        System.out.println("Will store for organisms:"+taxonIds);
    }

    // Set DataSource fields in project.xml
    public void setDataSourceName(String name) {
        this.dataSourceName = name;
    }
    public void setDataSourceUrl(String url) {
        this.dataSourceUrl = url;
    }
    public void setDataSourceDescription(String description) {
        this.dataSourceDescription = description;
    }

    // set DataSet fields in project.xml
    public void setDataSetUrl(String url) {
        this.dataSetUrl = url;
    }
    public void setDataSetDescription(String description) {
        this.dataSetDescription = description;
    }
        
    /**
     * Create a new DatastoreFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public DatastoreFileConverter(ItemWriter writer, Model model) throws ObjectStoreException {
        super(writer, model);
        // GO
        geneOntology = createItem("Ontology");
        geneOntology.setAttribute("name", "GO");
        geneOntology.setAttribute("url", "http://www.geneontology.org");
        store(geneOntology);
        // Pfam
        pfamOntology = createItem("Ontology");
        pfamOntology.setAttribute("name", "Pfam");
        pfamOntology.setAttribute("url", "https://pfam.xfam.org/");
        store(pfamOntology);
        // PANTHER
        pantherOntology = createItem("Ontology");
        pantherOntology.setAttribute("name", "PANTHER");
        pantherOntology.setAttribute("url", "http://www.pantherdb.org/");
        store(pantherOntology);
        // KOG
        kogOntology = createItem("Ontology");
        kogOntology.setAttribute("name", "KOG");
        kogOntology.setAttribute("url", "https://genome.jgi.doe.gov/Tutorial/tutorial/kog.html");
        store(kogOntology);
        // ENZYME
        ecOntology = createItem("Ontology");
        ecOntology.setAttribute("name", "ENZYME");
        ecOntology.setAttribute("url", "https://enzyme.expasy.org/");
        store(ecOntology);
        // KEGG Orthology
        koOntology = createItem("Ontology");
        koOntology.setAttribute("name", "KEGG Orthology");
        koOntology.setAttribute("url", "https://www.genome.jp/kegg/ko.html");
        store(koOntology);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void process(Reader reader) throws Exception {
        // create and store the DataSource if appropriate
        if (dataSource==null && dataSourceName!=null) {
            dataSource = createItem("DataSource");
            dataSource.setAttribute("name", dataSourceName);
            if (dataSourceUrl!=null) dataSource.setAttribute("url", dataSourceUrl);
            if (dataSourceDescription!=null) dataSource.setAttribute("description", dataSourceDescription);
            store(dataSource);
        }

        // process the file
        if (getCurrentFile().getName().contains("description_")) {
            // description_Phaseolus_vulgaris.yml
            printInfoBlurb(getCurrentFile().getName());
            processDescriptionFile(reader);
        } else if (getCurrentFile().getName().contains("strains_")) {
            // strains_Phaseolus_vulgaris.yml
            printInfoBlurb(getCurrentFile().getName());
            processStrainsFile(reader);
        } else if (getCurrentFile().getName().endsWith(".info_annot.txt")) {
            // phavu.G19833.gnm2.ann1.PB8d.info_annot.txt
            printInfoBlurb(getCurrentFile().getName());
            processInfoAnnotFile(reader);
        } else if (getCurrentFile().getName().endsWith(".info_annot_ahrd.tsv")) {
            // legume.genefam.fam1.M65K.info_annot_ahrd.tsv
            String fastaDirname = getCurrentFile().getParent()+"/"+getCurrentFile().getName().replace("info_annot_ahrd.tsv", "family_fasta");
            printInfoBlurb(getCurrentFile().getName());
            printInfoBlurb(fastaDirname);
            processInfoAnnotAhrdFile(reader);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void close() throws Exception {
        System.out.println("Storing all maps...");
        // straight maps
        store(organisms.values());
        store(strains.values());
        store(ontologyTerms.values());
        store(dataSets.values());
        store(geneFamilies.values());
        store(proteinDomains.values());
        // maps of maps
        for (Map<String,Item> geneMap : genes.values()) {
            store(geneMap.values());
        }
        for (Map<String,Item> proteinMap : proteins.values()) {
            store(proteinMap.values());
        }
        for (Map<String,Item> mRNAMap : mRNAs.values()) {
            store(mRNAMap.values());
        }
    }

    /**
     * Print out info about the current file or directory being processed.
     */
    static void printInfoBlurb(String blurb) {
        LOG.info("Processing file/dir "+blurb);
        System.out.println("####################################################################################################################");
        System.out.println("Processing file/dir "+blurb);
        System.out.println("####################################################################################################################");
    }

    /**
     * Get/add the DataSet Item, formed from the current filename.
     */
    Item getDataSet() {
        if (dataSets.containsKey(getCurrentFile().getName())) {
            return dataSets.get(getCurrentFile().getName());
        } else {
            Item dataSet = createItem("DataSet");
            dataSet.setAttribute("name", getCurrentFile().getName());
            if (dataSetUrl!=null) dataSet.setAttribute("url", dataSetUrl);
            if (dataSetDescription!=null) dataSet.setAttribute("description", dataSetDescription);
            if (dataSource!=null) dataSet.setReference("dataSource", dataSource);
            dataSets.put(getCurrentFile().getName(), dataSet);
            return dataSet;
        }
    }

    /**
     * Get/add the DataSet Item, with a given name (when the current filename is not appropriate).
     */
    Item getDataSet(String name) {
        if (dataSets.containsKey(name)) {
            return dataSets.get(name);
        } else {
            Item dataSet = createItem("DataSet");
            dataSet.setAttribute("name", name);
            if (dataSetUrl!=null) dataSet.setAttribute("url", dataSetUrl);
            if (dataSetDescription!=null) dataSet.setAttribute("description", dataSetDescription);
            if (dataSource!=null) dataSet.setReference("dataSource", dataSource);
            dataSets.put(name, dataSet);
            return dataSet;
        }
    }

    /**
     * Get the taxon ID for a string array like ["something","Phaseolus","vulgaris"]
     */
    String getTaxonId(String[] dashparts) {
        String genus = dashparts[1];
        String species = dashparts[2];
        String taxonId = genusSpeciesTaxonIds.get(genus+"_"+species);
        if (taxonId==null) {
            throw new RuntimeException("Taxon ID not available for "+genus+"_"+species);
        }
        return taxonId;
    }

    /**
     * Get the taxon ID for a gensp string like "phavu".
     */
    String getTaxonId(String gensp) {
        String taxonId = genspTaxonIds.get(gensp);
        if (taxonId==null) {
            throw new RuntimeException("Taxon ID not available for "+gensp);
        }
        return taxonId;
    }

    /**
     * Get/add the organism Item associated with the given taxon ID.
     */
    Item getOrganism(String taxonId) {
        Item organism;
        if (organisms.containsKey(taxonId)) {
            organism = organisms.get(taxonId);
        } else {
            organism = createItem("Organism");
            organism.setAttribute("taxonId", taxonId);
            organisms.put(taxonId, organism);
        }
        return organism;
    }

    /**
     * Get/add the strain Item associated with the given strain name.
     * Sets the organism reference if created.
     */
    Item getStrain(String strainName, Item organism) {
        Item strain;
        if (strains.containsKey(strainName)) {
            strain = strains.get(strainName);
        } else {
            strain = createItem("Strain");
            strain.setAttribute("primaryIdentifier", strainName);
            strain.setReference("organism", organism);
            strains.put(strainName, strain);
        }
        return strain;
    }

    /**
     * Process an organism description file, which is in YAML format:
     *
     * description_Phaseolus_vulgaris.yml
     *
     * %YAML 1.2
     * ##### Phaseolus vulgaris	
     * organism.taxid:	3885
     * organism.genus:	Phaseolus
     * organism.species:	vulgaris
     * organism.abbrev:	phavu
     * organism.commonName:	common bean
     * organism.description:	Common bean was likely domesticated independently both in Central America and in the Andes....
     */
    void processDescriptionFile(Reader reader) throws IOException {
        // get the organism
        String[] dotparts = getCurrentFile().getName().split("\\.");
        String[] dashparts = dotparts[0].split("_");
        Item organism = getOrganism(getTaxonId(dashparts));
        // now load the attributes
        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("%")) continue;
            if (line.startsWith("#")) continue;
            String[] parts = line.split("\t");
            if (parts.length>1) {
                String attributeName = parts[0].replace("organism.","").replace(":","");
                // munge
                if (attributeName.equals("taxid")) attributeName = "taxonId";
                if (attributeName.equals("abbrev")) attributeName = "abbreviation";
                String attributeValue = parts[1].trim();
                if (attributeValue.length()>0) {
                    organism.setAttribute(attributeName, attributeValue);
                }
            }
        }
        br.close();
    }

    /**
     * Process a strains file, which is in YAML format:
     *
     * strains_Phaseolus_vulgaris.yml
     *
     * %YAML 1.2
     * #####
     * strain.identifier:	G19833
     * strain.accession:	
     * strain.name:	G19833
     * strain.origin:	Peru
     * strain.description:	Andean landrace G19833 was selected for genome sequencing partly due to its resistance to numerous diseases...
     * #####
     * strain.identifier:	BAT93
     * strain.accession:	PI 633451
     * strain.name:	BAT93
     * strain.origin:	CIAT
     * strain.description:	Accession BAT93 is a Mesomarican line that has been used in numerous breeding projects and trait-mapping studies.
     */
    void processStrainsFile(Reader reader) throws IOException {
        // get the organism
        String[] dotparts = getCurrentFile().getName().split("\\.");
        String[] dashparts = dotparts[0].split("_");
        Item organism = getOrganism(getTaxonId(dashparts));
        // spin through the strain sections
        BufferedReader br = new BufferedReader(reader);
        String strainName = null; // keep track of current strain name
        String line = null;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.startsWith("%") || line.startsWith("#####")) {
                continue;
            } else {
                // strain attributes
                String[] parts = line.split("\t");
                if (parts.length>1) {
                    String attributeName = parts[0].replace("strain.","").replace(":","");
                    String attributeValue = parts[1].trim();
                    if (attributeName.equals("identifier")) {
                        strainName = attributeValue;
                    }
                    Item strain = getStrain(strainName, organism);
                    // special case identifier --> primaryIdentifier (Annotatable)
                    if (attributeName.equals("identifier")) attributeName = "primaryIdentifier";
                    if (attributeValue.length()>0) {
                        strain.setAttribute(attributeName, attributeValue);
                    }
                }
            }
        }
        br.close();
    }

    /**
     * Process an info_annot.txt file which contains relationships between genes, transcripts, proteins and ontology terms.
     * This will also link genes to proteins. The file name starts with gensp.strainName.assemblyVersion.annotationVersion.
     * 0     1      2    3    4    5          6
     * phavu.G19833.gnm2.ann1.PB8d.info_annot.txt
     *
     * pacId locusName transcriptName peptideName Pfam Panther KOG ec KO GO Best-hit-arabi-name arabi-symbol arabi-defline
     * 37170591 Phvul.001G000400 Phvul.001G000400.1 Phvul.001G000400.1.p PF00504 PTHR21649,PTHR21649:SF24 1.10.3.9 K14172 GO:0016020,GO:0009765 AT1G76570.1 Chlorophyll family protein
     */
    void processInfoAnnotFile(Reader reader) throws IOException, ObjectStoreException {
        // get the dot-separated file information
        String[] dotparts = getCurrentFile().getName().split("\\.");
        String gensp = dotparts[0];
        String strainName = dotparts[1];
        String assemblyVersion = dotparts[2];
        String annotationVersion = dotparts[3];
        String versionKey = assemblyVersion+"."+annotationVersion;
        // get the dataSet, organism and strain
        Item dataSet = getDataSet();
        dataSet.setAttribute("version", versionKey);
        Item organism = getOrganism(getTaxonId(gensp));
        Item strain = getStrain(strainName, organism);
        // spin through the file
        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            // comment line
            if (line.startsWith("#")) continue;
            // String bestHitAtName;
            // String bestHitAtSymbol;
            // String bestHitAtDefline;
            InfoAnnotRecord record = new InfoAnnotRecord(line);
            if (record.pacId!=null) {
                // the gene
                Item gene = null;
                String geneIdentifier = record.locusName;
                if (genes.containsKey(geneIdentifier)) {
                    Map<String,Item> geneMap = genes.get(geneIdentifier);
                    if (geneMap.containsKey(versionKey)) {
                        gene = geneMap.get(versionKey);
                    } else {
                        // Phvul.002G040500
                        gene = createItem("Gene");
                        gene.setAttribute("secondaryIdentifier", geneIdentifier);
                        gene.setReference("organism", organism);
                        gene.setReference("strain", strain);
                        gene.addToCollection("dataSets", dataSet);
                        gene.setAttribute("assemblyVersion", assemblyVersion);
                        gene.setAttribute("annotationVersion", annotationVersion);
                        geneMap.put(versionKey, gene);
                    }
                } else {
                    // Phvul.002G040500
                    gene = createItem("Gene");
                    gene.setAttribute("secondaryIdentifier", geneIdentifier);
                    gene.setReference("organism", organism);
                    gene.setReference("strain", strain);
                    gene.addToCollection("dataSets", dataSet);
                    gene.setAttribute("assemblyVersion", assemblyVersion);
                    gene.setAttribute("annotationVersion", annotationVersion);
                    Map<String,Item> geneMap = new HashMap<>();
                    geneMap.put(versionKey, gene);
                    genes.put(geneIdentifier, geneMap);
                }
                // the protein
                Item protein = null;
                String proteinIdentifier = record.peptideName;
                if (proteins.containsKey(proteinIdentifier)) {
                    Map<String,Item> proteinMap = proteins.get(proteinIdentifier);
                    if (proteinMap.containsKey(versionKey)) {
                        protein = proteinMap.get(versionKey);
                    } else {
                        // Phvul.002G040500.1.p --> Phvul.002G040500.1
                        protein = createItem("Protein");
                        protein.setAttribute("secondaryIdentifier", proteinIdentifier);
                        protein.setReference("organism", organism);
                        protein.setReference("strain", strain);
                        protein.addToCollection("dataSets", dataSet);
                        protein.setAttribute("assemblyVersion", assemblyVersion);
                        protein.setAttribute("annotationVersion", annotationVersion);
                        protein.addToCollection("genes", gene);
                        proteinMap.put(versionKey, protein);
                    }
                } else {
                    // Phvul.002G040500.1.p --> Phvul.002G040500.1
                    protein = createItem("Protein");
                    protein.setAttribute("secondaryIdentifier", proteinIdentifier);
                    protein.setReference("organism", organism);
                    protein.setReference("strain", strain);
                    protein.addToCollection("dataSets", dataSet);
                    protein.setAttribute("assemblyVersion", assemblyVersion);
                    protein.setAttribute("annotationVersion", annotationVersion);
                    protein.addToCollection("genes", gene);
                    Map<String,Item> proteinMap = new HashMap<>();
                    proteinMap.put(versionKey, protein);
                    proteins.put(proteinIdentifier, proteinMap);
                }
                // the transcript = mRNA
                Item mRNA = null;
                String mRNAIdentifier = record.transcriptName;
                if (mRNAs.containsKey(mRNAIdentifier)) {
                    Map<String,Item> mRNAMap = mRNAs.get(mRNAIdentifier);
                    if (mRNAMap.containsKey(versionKey)) {
                        mRNA = mRNAMap.get(versionKey);
                    } else {
                        // Phvul.002G040500.1
                        mRNA = createItem("MRNA");
                        mRNA.setAttribute("secondaryIdentifier", mRNAIdentifier);
                        mRNA.setReference("organism", organism);
                        mRNA.setReference("strain", strain);
                        mRNA.addToCollection("dataSets", dataSet);
                        mRNA.setAttribute("assemblyVersion", assemblyVersion);
                        mRNA.setAttribute("annotationVersion", annotationVersion);
                        mRNA.setReference("gene", gene);
                        mRNA.setReference("protein", protein);
                        mRNAMap.put(versionKey, mRNA);
                    }
                } else {
                    // Phvul.002G040500.1
                    mRNA = createItem("MRNA");
                    mRNA.setAttribute("secondaryIdentifier", mRNAIdentifier);
                    mRNA.setReference("organism", organism);
                    mRNA.setReference("strain", strain);
                    mRNA.addToCollection("dataSets", dataSet);
                    mRNA.setAttribute("assemblyVersion", assemblyVersion);
                    mRNA.setAttribute("annotationVersion", annotationVersion);
                    mRNA.setReference("gene", gene);
                    mRNA.setReference("protein", protein);
                    Map<String,Item> mRNAMap = new HashMap<>();
                    mRNAMap.put(versionKey, mRNA);
                    mRNAs.put(mRNAIdentifier, mRNAMap);
                }
                // GO terms
                for (String identifier : record.GO) {
                    Item goTerm = null;
                    if (ontologyTerms.containsKey(identifier)) {
                        goTerm = ontologyTerms.get(identifier);
                    } else {
                        goTerm = createItem("OntologyTerm");
                        goTerm.setAttribute("identifier", identifier);
                        goTerm.setReference("ontology", geneOntology);
                        ontologyTerms.put(identifier, goTerm);
                    }
                    String annotKey = identifier+"_"+versionKey+"_"+geneIdentifier;
                    if (!ontologyAnnotations.contains(annotKey)) {
                        Item goAnnotation = createItem("OntologyAnnotation");
                        goAnnotation.setReference("subject", gene);
                        goAnnotation.setReference("ontologyTerm", goTerm);
                        goAnnotation.addToCollection("dataSets", dataSet);
                        store(goAnnotation);
                        ontologyAnnotations.add(annotKey);
                    }
                }
                // Pfam terms
                for (String identifier : record.Pfam) {
                    Item pfamTerm = null;
                    if (ontologyTerms.containsKey(identifier)) {
                        pfamTerm = ontologyTerms.get(identifier);
                    } else {
                        pfamTerm = createItem("OntologyTerm");
                        pfamTerm.setAttribute("identifier", identifier);
                        pfamTerm.setReference("ontology", pfamOntology);
                        ontologyTerms.put(identifier, pfamTerm);
                    }
                    String annotKey = identifier+"_"+versionKey+"_"+proteinIdentifier;
                    if (!ontologyAnnotations.contains(annotKey)) {
                        Item pfamAnnotation = createItem("OntologyAnnotation");
                        pfamAnnotation.setReference("subject", protein);
                        pfamAnnotation.setReference("ontologyTerm", pfamTerm);
                        pfamAnnotation.addToCollection("dataSets", dataSet);
                        store(pfamAnnotation);
                        ontologyAnnotations.add(annotKey);
                    }
                }
                // Panther terms
                for (String identifier : record.Panther) {
                    Item pantherTerm = null;
                    if (ontologyTerms.containsKey(identifier)) {
                        pantherTerm = ontologyTerms.get(identifier);
                    } else {
                        pantherTerm = createItem("OntologyTerm");
                        pantherTerm.setAttribute("identifier", identifier);
                        pantherTerm.setReference("ontology", pantherOntology);
                        ontologyTerms.put(identifier, pantherTerm);
                    }
                    String annotKey = identifier+"_"+versionKey+"_"+proteinIdentifier;
                    if (!ontologyAnnotations.contains(annotKey)) {
                        Item pantherAnnotation = createItem("OntologyAnnotation");
                        pantherAnnotation.setReference("subject", protein);
                        pantherAnnotation.setReference("ontologyTerm", pantherTerm);
                        pantherAnnotation.addToCollection("dataSets", dataSet);
                        store(pantherAnnotation);
                        ontologyAnnotations.add(annotKey);
                    }                }
                // KOG terms
                for (String identifier : record.KOG) {
                    Item kogTerm = null;
                    if (ontologyTerms.containsKey(identifier)) {
                        kogTerm = ontologyTerms.get(identifier);
                    } else {
                        kogTerm = createItem("OntologyTerm");
                        kogTerm.setAttribute("identifier", identifier);
                        kogTerm.setReference("ontology", kogOntology);
                        ontologyTerms.put(identifier, kogTerm);
                    }
                    String annotKey = identifier+"_"+versionKey+"_"+proteinIdentifier;
                    if (!ontologyAnnotations.contains(annotKey)) {
                        Item kogAnnotation = createItem("OntologyAnnotation");
                        kogAnnotation.setReference("subject", protein);
                        kogAnnotation.setReference("ontologyTerm", kogTerm);
                        kogAnnotation.addToCollection("dataSets", dataSet);
                        store(kogAnnotation);
                        ontologyAnnotations.add(annotKey);
                    }
                }
                // ec terms
                for (String identifier : record.ec) {
                    Item ecTerm = null;
                    if (ontologyTerms.containsKey(identifier)) {
                        ecTerm = ontologyTerms.get(identifier);
                    } else {
                        ecTerm = createItem("OntologyTerm");
                        ecTerm.setAttribute("identifier", identifier);
                        ecTerm.setReference("ontology", ecOntology);
                        ontologyTerms.put(identifier, ecTerm);
                    }
                    String annotKey = identifier+"_"+versionKey+"_"+proteinIdentifier;
                    if (!ontologyAnnotations.contains(annotKey)) {
                        Item ecAnnotation = createItem("OntologyAnnotation");
                        ecAnnotation.setReference("subject", protein);
                        ecAnnotation.setReference("ontologyTerm", ecTerm); 
                        ecAnnotation.addToCollection("dataSets", dataSet);
                        store(ecAnnotation);
                        ontologyAnnotations.add(annotKey);
                    }
                }
                // KO terms
                for (String identifier : record.KO) {
                    Item koTerm = null;
                    if (ontologyTerms.containsKey(identifier)) {
                        koTerm = ontologyTerms.get(identifier);
                    } else {
                        koTerm = createItem("OntologyTerm");
                        koTerm.setAttribute("identifier", identifier);
                        koTerm.setReference("ontology", koOntology);
                        ontologyTerms.put(identifier, koTerm);
                    }
                    String annotKey = identifier+"_"+versionKey+"_"+geneIdentifier;
                    if (!ontologyAnnotations.contains(annotKey)) {
                        Item koAnnotation = createItem("OntologyAnnotation");
                        koAnnotation.setReference("subject", gene);
                        koAnnotation.setReference("ontologyTerm", koTerm);
                        koAnnotation.addToCollection("dataSets", dataSet);
                        store(koAnnotation);
                        ontologyAnnotations.add(annotKey);
                    }
                }
            }
        }
        br.close();
    }

    /**
     * Process an info_annot_ahrd.tsv file which contains gene families and semi-colon separated groups of ontology terms.
     * 0     1       2    3    4               5
     * zzzzz.genefam.fam1.M65K.info_annot_ahrd.tsv
     * zzzzz.genefam.fam1.M65K.family_fasta/          contains individual gene family FASTAs
     *
     * legfed_v1_0.L_H6Q7T0-consensus type I inositol-1,4,5-trisphosphate 5-phosphatase;
     *             ^^^^^^^^           IPR005135 (Endonuclease/exonuclease/phosphatase), IPR015943 (WD40/YVTN repeat-like-containing domain);
     *                                GO:0005515 (protein binding), GO:0046856 (phosphatidylinositol dephosphorylation)
     *
     * zzzzz.genefam.fam1.M65K.family_fasta/L_H6Q7T0
     *                                      ^^^^^^^^
     * >lupan.Lup015831.1
     * ----CASFAKLT--TLSPHWIGNNSFSSRRGGSSPLTATRRVSLPIRASSYSDELVQTAK
     * TIASPGRGILAIDESNATCGKRLASIGLDNTEVNRQAYRQLLLTTPGLGEYISGAILFEE
     * ...
     * >phavu.Phvul.007G033800.1
     * -----------------------------------TFSPRRVSLPIRASSYQQELVQTAK
     * SIASPGRGILAIDESNATCGKRLASIGLDNTEVNRQAYRQLLLTTPGLGEYISGAILFEE
     * ...
     */
    void processInfoAnnotAhrdFile(Reader reader) throws IOException, ObjectStoreException {
        // get the dataSet version from the filename
        String[] dotparts = getCurrentFile().getName().split("\\.");
        String dataSetVersion = dotparts[2];
        Item dataSet = getDataSet();
        dataSet.setAttribute("version", dataSetVersion);
        // get the FASTA directory
        String fastaDirname = getCurrentFile().getParent()+"/"+getCurrentFile().getName().replace("info_annot_ahrd.tsv", "family_fasta");
        // spin through the lines
        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            // comment line
            if (line.startsWith("#")) continue;
            // parse record and create items
            InfoAnnotAhrdRecord record = new InfoAnnotAhrdRecord(line);
            Item geneFamily = createItem("GeneFamily");
            geneFamily.setAttribute("identifier", record.identifier);
            geneFamily.setAttribute("version", record.version);
            geneFamily.setAttribute("description", record.description);
            geneFamily.setReference("dataSet", dataSet);
            geneFamilies.put(record.identifier, geneFamily);
            // GO terms
            for (String identifier : record.go.keySet()) {
                String description = record.go.get(identifier);
                Item goTerm = null;
                if (ontologyTerms.containsKey(identifier)) {
                    goTerm = ontologyTerms.get(identifier);
                } else {
                    goTerm = createItem("OntologyTerm");
                    goTerm.setAttribute("identifier", identifier);
                    goTerm.setAttribute("description", description);
                    goTerm.setReference("ontology", geneOntology);
                    ontologyTerms.put(identifier, goTerm);
                }
                Item goAnnotation = createItem("OntologyAnnotation");
                goAnnotation.setReference("subject", geneFamily);
                goAnnotation.setReference("ontologyTerm", goTerm);
                goAnnotation.addToCollection("dataSets", dataSet);
                store(goAnnotation);
            }
            // interpro domains
            for (String identifier : record.interpro.keySet()) {
                String description = record.interpro.get(identifier);
                Item proteinDomain = null;
                if (proteinDomains.containsKey(identifier)) {
                    proteinDomain = proteinDomains.get(identifier);
                } else {
                    proteinDomain = createItem("ProteinDomain");
                    proteinDomain.setAttribute("primaryIdentifier", identifier);
                    proteinDomain.setAttribute("description", description);
                    proteinDomains.put(identifier, proteinDomain);
                }
                proteinDomain.addToCollection("geneFamilies", geneFamily);
            }
            // load the gene family FASTA if present to link proteins of the desired organism
            String fastaFilename = fastaDirname+"/"+record.identifier;
            File fastaFile = new File(fastaFilename);
            if (fastaFile.exists()) {
                BufferedReader fbr = new BufferedReader(new FileReader(fastaFile));
                String fline = null;
                while ((fline=fbr.readLine())!=null) {
                    if (fline.startsWith(">")) {
                        String genspName = fline.substring(1);
                        String[] parts = genspName.split("\\.");
                        String gensp = parts[0];
                        String name = parts[1];
                        for (int i=2; i<parts.length; i++) name += "."+parts[i];
                        String taxonId = genspTaxonIds.get(gensp);
                        if (taxonId!=null && taxonIds.contains(taxonId)) {
                            if (proteins.containsKey(name)) {
                                // update proteins we already have (hopefully most of them)
                                for (Item protein : proteins.get(name).values()) {
                                    protein.setReference("geneFamily", geneFamily);
                                    protein.addToCollection("dataSets", dataSet);
                                }
                            } else {
                                // create protein without assemblyVersion, annotationVersion if not already loaded (boo)
                                Item protein = createItem("Protein");
                                protein.setAttribute("secondaryIdentifier", name);
                                protein.setReference("geneFamily", geneFamily);
                                protein.addToCollection("dataSets", dataSet);
                                Map<String,Item> proteinMap = new HashMap<>();
                                proteinMap.put(dataSetVersion, protein);
                                proteins.put(name, proteinMap);
                            }
                        }
                    }
                }
                fbr.close();
            }
        }
        br.close();
    }
}
