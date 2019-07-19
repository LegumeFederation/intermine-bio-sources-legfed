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
import java.io.InputStream;

import java.util.Arrays;
import java.util.Map;
import java.util.HashMap;
import java.util.Set;
import java.util.HashSet;
import java.util.Properties;

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
    Map<String,Item> linkageGroups = new HashMap<>();
    Map<String,Item> geneticMarkers = new HashMap<>();
    Map<String,Item> chromosomes = new HashMap<>(); // contains both Chromosome and Supercontig Items
    Map<String,Item> genes = new HashMap<>();
    Map<String,Item> proteins = new HashMap<>();
    Map<String,Item> mRNAs = new HashMap<>();
    
    Set<String> ontologyAnnotations = new HashSet<>(); // concatenation of subject and term identifiers to avoid dupes


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

    // map gensp and Genus_species to taxonId, etc.
    Map<String,String> genspTaxonId = new HashMap<>();
    Map<String,String> taxonIdGenus = new HashMap<>();
    Map<String,String> taxonIdSpecies = new HashMap<>();

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

        // get the gensp to taxonId map and other organism deets
        DatastoreUtils dsu = new DatastoreUtils();
        genspTaxonId = dsu.getGenspTaxonId();
        taxonIdGenus = dsu.getTaxonIdGenus();
        taxonIdSpecies = dsu.getTaxonIdSpecies();
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
        } else if (getCurrentFile().getName().endsWith(".cmap.txt")) {
            // phavu.mixed.map1.7PMp.cmap.txt
            printInfoBlurb(getCurrentFile().getName());
            processCmapFile(reader);
        } else if (getCurrentFile().getName().endsWith(".map.gff3")) {
            // phavu.mixed.map1.7PMp.map.gff3
            printInfoBlurb(getCurrentFile().getName());
            processMapGFFFile(reader);
        } else if (getCurrentFile().getName().endsWith(".flanking_seq.fna")) {
            // phavu.mixed.map1.7PMp.flanking_seq.fna
            printInfoBlurb(getCurrentFile().getName());
            processFlankingSeqFastaFile(reader);
        } else if (getCurrentFile().getName().endsWith(".info_descriptors.txt")) {
            // arahy.Tifrunner.gnm1.ann1.CCJH.info_descriptors.txt
            printInfoBlurb(getCurrentFile().getName());
            processInfoDescriptorsFile(reader);
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
        store(linkageGroups.values());
        store(geneticMarkers.values());
        store(chromosomes.values());
        store(genes.values());
        store(proteins.values());
        store(mRNAs.values());
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
     * Get the taxon ID for a gensp string like "phavu".
     */
    String getTaxonId(String gensp) {
        String taxonId = genspTaxonId.get(gensp);
        if (taxonId==null) {
            throw new RuntimeException("Taxon ID not available for "+gensp);
        }
        return taxonId;
    }

    /**
     * Get the Genus for a gensp string like "phavu".
     */
    String getGenus(String gensp) {
        String taxonId = getTaxonId(gensp);
        if (taxonIdGenus.containsKey(taxonId)) {
            return taxonIdGenus.get(taxonId);
        } else {
            throw new RuntimeException("Genus not available for taxon ID "+taxonId);
        }
    }
    
    /**
     * Get the species for a gensp string like "phavu".
     */
    String getSpecies(String gensp) {
        String taxonId = genspTaxonId.get(gensp);
        if (taxonIdSpecies.containsKey(taxonId)) {
            return taxonIdSpecies.get(taxonId);
        } else {
            throw new RuntimeException("Species not available for taxon ID "+taxonId);
        }
    }

    /**
     * Get/add the organism Item associated with the given gensp value (e.g. "phavu").
     * Returns null if the gensp isn't resolvable.
     */
    Item getOrganism(String gensp) {
        Item organism = null;
        if (organisms.containsKey(gensp)) {
            organism = organisms.get(gensp);
        } else {
            String taxonId = getTaxonId(gensp);
            String genus = getGenus(gensp);
            String species = getSpecies(gensp);
            organism = createItem("Organism");
            organism.setAttribute("abbreviation", gensp);
            organism.setAttribute("taxonId", taxonId);
            organism.setAttribute("genus", genus);
            organism.setAttribute("species", species);
            organism.setAttribute("name", genus+" "+species);
            organisms.put(gensp, organism);
        }
        return organism;
    }

    /**
     * Get the organism Item for a string array like ["something","Genus","species"]
     */
    Item getOrganism(String[] threeparts) {
        String genus = threeparts[1];
        String species = threeparts[2];
        String gensp = genus.toLowerCase().substring(0,3) + species.toLowerCase().substring(0,2);
        return getOrganism(gensp);
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
     * Get/add a Gene Item, keyed by secondaryIdentifier (!)
     */
    public Item getGene(String secondaryIdentifier) {
        if (secondaryIdentifier==null || secondaryIdentifier.trim().length()==0) {
            throw new RuntimeException("secondaryIdentifier is null or empty in getGene.");
        }
        Item gene;
        if (genes.containsKey(secondaryIdentifier)) {
            gene = genes.get(secondaryIdentifier);
        } else {
            // phavu.Phvul.002G040500
            gene = createItem("Gene");
            gene.setAttribute("secondaryIdentifier", secondaryIdentifier);
            genes.put(secondaryIdentifier, gene);
        }
        return gene;
    }

    /**
     * Get/add a Protein Item, keyed by secondaryIdentifier (!)
     */
    public Item getProtein(String secondaryIdentifier) {
        Item protein;
        if (proteins.containsKey(secondaryIdentifier)) {
            protein = proteins.get(secondaryIdentifier);
        } else {
            // Phvul.002G040500.1.p --> Phvul.002G040500.1
            protein = createItem("Protein");
            protein.setAttribute("secondaryIdentifier", secondaryIdentifier);
            proteins.put(secondaryIdentifier, protein);
        }
        return protein;
    }

    /**
     * Get/add a ProteinDomain Item.
     */
    public Item getProteinDomain(String identifier) {
        Item proteinDomain;
        if (proteinDomains.containsKey(identifier)) {
            proteinDomain = proteinDomains.get(identifier);
        } else {
            proteinDomain = createItem("ProteinDomain");
            proteinDomain.setAttribute("primaryIdentifier", identifier);
        }
        return proteinDomain;
    }

    /**
     * Get/add an OntologyTerm Item, keyed by identifier
     */
    public Item getOntologyTerm(String identifier) {
        Item ontologyTerm;
        if (ontologyTerms.containsKey(identifier)) {
            ontologyTerm = ontologyTerms.get(identifier);
        } else {
            ontologyTerm = createItem("OntologyTerm");
            ontologyTerm.setAttribute("identifier", identifier);
            ontologyTerms.put(identifier, ontologyTerm);
        }
        return ontologyTerm;
    }

    /**
     * Get/add an MRNA Item, keyed by secondaryIdentifier (!)
     */
    public Item getMRNA(String secondaryIdentifier) {
        Item mRNA;
        if (mRNAs.containsKey(secondaryIdentifier)) {
            mRNA = mRNAs.get(secondaryIdentifier);
        } else {
            // Phvul.002G040500.1
            mRNA = createItem("MRNA");
            mRNA.setAttribute("secondaryIdentifier", secondaryIdentifier);
            mRNAs.put(secondaryIdentifier, mRNA);
        }
        return mRNA;
    }

    /**
     * Get/add a GeneFamily Item.
     */
    public Item getGeneFamily(String identifier) {
        Item geneFamily;
        if (geneFamilies.containsKey(identifier)) {
            geneFamily = geneFamilies.get(identifier);
        } else {
            geneFamily = createItem("GeneFamily");
            geneFamily.setAttribute("identifier", identifier);
            geneFamilies.put(identifier, geneFamily);
        }
        return geneFamily;
    }
    
    /**
     * Get/add a LinkageGroup Item
     */
    public Item getLinkageGroup(String primaryIdentifier) {
        if (linkageGroups.containsKey(primaryIdentifier)) {
            return linkageGroups.get(primaryIdentifier);
        } else {
            Item lg = createItem("LinkageGroup");
            lg.setAttribute("primaryIdentifier", primaryIdentifier);
            linkageGroups.put(primaryIdentifier, lg);
            return lg;
        }
    }

    /**
     * Get/add a GeneticMarker Item
     */
    public Item getGeneticMarker(String primaryIdentifier) {
        if (geneticMarkers.containsKey(primaryIdentifier)) {
            return geneticMarkers.get(primaryIdentifier);
        } else {
            Item gm = createItem("GeneticMarker");
            gm.setAttribute("primaryIdentifier", primaryIdentifier);
            geneticMarkers.put(primaryIdentifier, gm);
            return gm;
        }
    }

    /**
     * Get/add a Chromosome/Supercontig Item
     */
    public Item getChromosomeOrSupercontig(String primaryIdentifier, Item organism) {
        if (chromosomes.containsKey(primaryIdentifier)) {
            return chromosomes.get(primaryIdentifier);
        } else {
            Item chr;
            // HACK: change the className to "Supercontig" if identifier contains "scaffold" or ends in "sc" etc.
            String[] dotparts = primaryIdentifier.split("\\.");
            String lastPart = dotparts[dotparts.length-1];
            if (primaryIdentifier.toLowerCase().contains("scaffold") || lastPart.contains("sc")) {
                chr = createItem("Supercontig");
                // DEBUG
                System.out.println("Adding Supercontig "+primaryIdentifier);
            } else {
                chr = createItem("Chromosome");
                // DEBUG
                System.out.println("Adding Chromosome "+primaryIdentifier);
            }
            chr.setAttribute("primaryIdentifier", primaryIdentifier);
            chr.setReference("organism", organism);
            chromosomes.put(primaryIdentifier, chr);
            return chr;
        }
    }

    ////////////////////////////////////////////////////////
    //////////////////// FILE PROCESSORS ///////////////////
    ////////////////////////////////////////////////////////

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
        Item organism = getOrganism(dashparts);
        String genus = null;   // for forming Organism.name
        String species = null; // for forming Organism.name
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
                    if (attributeName.equals("genus")) genus = attributeValue;
                    if (attributeName.equals("species")) species = attributeValue;
                }
            }
        }
        br.close();
        if (genus!=null && species!=null) organism.setAttribute("name", genus+" "+species);
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
        Item organism = getOrganism(dashparts);
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
        Item organism = getOrganism(gensp);
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
                String geneIdentifier = gensp+"."+record.locusName;
                Item gene = getGene(geneIdentifier);
                gene.setAttribute("assemblyVersion", assemblyVersion);
                gene.setAttribute("annotationVersion", annotationVersion);
                gene.setReference("organism", organism);
                gene.setReference("strain", strain);
                gene.addToCollection("dataSets", dataSet);
                // the protein
                String proteinIdentifier = gensp+"."+record.peptideName;
                Item protein = getProtein(proteinIdentifier);
                protein.setAttribute("assemblyVersion", assemblyVersion);
                protein.setAttribute("annotationVersion", annotationVersion);
                protein.setReference("organism", organism);
                protein.setReference("strain", strain);
                protein.addToCollection("genes", gene);
                protein.addToCollection("dataSets", dataSet);
                // the transcript = mRNA
                String mRNAIdentifier = gensp+"."+record.transcriptName;
                Item mRNA = getMRNA(mRNAIdentifier);
                mRNA.setAttribute("assemblyVersion", assemblyVersion);
                mRNA.setAttribute("annotationVersion", annotationVersion);
                mRNA.setReference("gene", gene); 
                mRNA.setReference("organism", organism);
                mRNA.setReference("strain", strain);
                mRNA.addToCollection("dataSets", dataSet);
                mRNA.setReference("protein", protein);
                // GO terms
                for (String identifier : record.GO) {
                    Item goTerm = getOntologyTerm(identifier);
                    goTerm.setReference("ontology", geneOntology);
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
                    Item pfamTerm = getOntologyTerm(identifier);
                    pfamTerm.setReference("ontology", pfamOntology);
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
                    Item pantherTerm = getOntologyTerm(identifier);
                    pantherTerm.setReference("ontology", pantherOntology);
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
                    Item kogTerm = getOntologyTerm(identifier);
                    kogTerm.setReference("ontology", kogOntology);
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
                    Item ecTerm = getOntologyTerm(identifier);
                    ecTerm.setReference("ontology", ecOntology);
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
                    Item koTerm = getOntologyTerm(identifier);
                    koTerm.setReference("ontology", koOntology);
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
     * 0      1       2    3    4               5
     * legfed.genefam.fam1.M65K.info_annot_ahrd.tsv
     *
     * legfed_v1_0.L_H6Q7T0-consensus type I inositol-1,4,5-trisphosphate 5-phosphatase;
     *             ^^^^^^^^           IPR005135 (Endonuclease/exonuclease/phosphatase), IPR015943 (WD40/YVTN repeat-like-containing domain);
     *                                GO:0005515 (protein binding), GO:0046856 (phosphatidylinositol dephosphorylation)
     *
     * legfed.genefam.fam1.M65K.family_fasta/L_H6Q7T0
     *                                       ^^^^^^^^
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
        String dataSetName = getCurrentFile().getName();
        String dataSetVersion = dotparts[2];
        Item dataSet = getDataSet(dataSetName);
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
            Item geneFamily = getGeneFamily(record.identifier);
            geneFamily.setAttribute("version", record.version);
            geneFamily.setAttribute("description", record.description);
            geneFamily.setReference("dataSet", dataSet);
            // GO terms
            for (String identifier : record.go.keySet()) {
                String description = record.go.get(identifier);
                Item goTerm = getOntologyTerm(identifier);
                goTerm.setAttribute("description", description);
                goTerm.setReference("ontology", geneOntology);
                Item goAnnotation = createItem("OntologyAnnotation");
                goAnnotation.setReference("subject", geneFamily);
                goAnnotation.setReference("ontologyTerm", goTerm);
                goAnnotation.addToCollection("dataSets", dataSet);
                store(goAnnotation);
            }
            // interpro domains
            for (String identifier : record.interpro.keySet()) {
                Item proteinDomain = getProteinDomain(identifier);
                String description = record.interpro.get(identifier);
                proteinDomain.setAttribute("description", description);
                proteinDomain.addToCollection("geneFamilies", geneFamily);
            }
            // load the gene family FASTA if present to link proteins
            String fastaFilename = fastaDirname+"/"+record.identifier;
            File fastaFile = new File(fastaFilename);
            if (fastaFile.exists()) {
                BufferedReader fbr = new BufferedReader(new FileReader(fastaFile));
                String fline = null;
                while ((fline=fbr.readLine())!=null) {
                    if (fline.startsWith(">")) {
                        String name = fline.substring(1);
                        String[] parts = name.split("\\.");
                        String gensp = parts[0];
                        Item organism = getOrganism(gensp);
                        Item protein = getProtein(name);
                        protein.setReference("geneFamily", geneFamily);
                        protein.addToCollection("dataSets", dataSet);
                    }
                }
                fbr.close();
            }
        }
        br.close();
    }

    /**
     * Process a genetic Cmap file.
     *
     * 0     1     2    3    4    5
     * phavu.mixed.map1.7PMp.cmap.txt
     *
     * map_acc                map_name map_start map_stop feature_acc                feature_name feature_aliases feature_start feature_stop feature_type_acc is_landmark
     * PvCookUCDavis2009_Pv09 Pv09     0         66.1     PvCook...Pv09_Pv_TOG913042 Pv_TOG913042 TOG913042       66.1          66.1         SNP              0
     * PvCookUCDavis2009_Pv09 Pv09     0         66.1     PvCook...Pv09_Pv_TOG899751 Pv_TOG899751 TOG899751       44.4          44.4         SNP              0
     */
    void processCmapFile(Reader reader) throws IOException, ObjectStoreException {
        // get the identifiers from the filename
        String[] dotparts = getCurrentFile().getName().split("\\.");
        String dataSetName = getCurrentFile().getName();
        String gensp = dotparts[0];
        String parents = dotparts[1];
        String dataSetVersion = dotparts[2];
        // create Items
        Item dataSet = getDataSet(dataSetName);
        dataSet.setAttribute("version", dataSetVersion);
        Item organism = getOrganism(gensp);
        // spin through the lines
        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            CMapRecord record = new CMapRecord(line);
            if (!record.hasData()) continue;
            Item lg = getLinkageGroup(record.map_acc);
            lg.setAttribute("secondaryIdentifier", record.map_name);
            lg.setAttribute("length", String.valueOf(record.map_stop));
            lg.setReference("organism", organism);
            lg.setReference("dataSet", dataSet);
            String fullname = record.feature_name;
            String name = fullname;
            if (record.feature_aliases!=null) name = record.feature_aliases;
            Item gm = getGeneticMarker(name);
            gm.setAttribute("type", String.valueOf(record.feature_type_acc));
            gm.setReference("organism", organism);
            gm.addToCollection("dataSets", dataSet);
            Item lgp = createItem("LinkageGroupPosition");
            lgp.setAttribute("position", String.valueOf(record.feature_start));
            lgp.setReference("linkageGroup", lg);
            store(lgp);
            gm.addToCollection("linkageGroupPositions", lgp);
            lg.addToCollection("markers", gm);
        }
        br.close();
    }

    /**
     * Process a flanking sequence FASTA file -- store the flanking sequence as the genetic marker's "sequence".
     *
     * 0     1     2    3    4            5
     * phavu.mixed.map1.7PMp.flanking_seq.fna
     *
     * >TOG905303_749
     * GACACGTAACTGAAATTTCACYATCCCTTACTGTTTCTCAAATCTTAGGATGAAAATAGATGCAATTGCTGGAAGTTCCCAATATTTTCTTTTRCACCTATGTACTCAGGTAATTTTATAT
     */
    void processFlankingSeqFastaFile(Reader reader) throws IOException, ObjectStoreException {
        // get the identifiers from the filename
        String[] dotparts = getCurrentFile().getName().split("\\.");
        String dataSetName = getCurrentFile().getName();
        String gensp = dotparts[0];
        String parents = dotparts[1];
        String dataSetVersion = dotparts[2];
        // create Items
        Item dataSet = getDataSet(dataSetName);
        dataSet.setAttribute("version", dataSetVersion);
        Item organism = getOrganism(gensp);
        // spin through the lines
        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            if (!line.startsWith(">")) continue;
            String fullname = line.substring(1);
            String name = fullname;
            String[] underscoreparts = fullname.split("_");
            if (underscoreparts.length==3) {
                // strip front and back of Pv_TOG905303_749
                name = underscoreparts[1];
            } else if (underscoreparts.length==2) {
                // strip back of TOG905303_749
                name = underscoreparts[0];
            }
            Item gm = getGeneticMarker(name);
            if (!name.equals(fullname)) gm.setAttribute("secondaryIdentifier", fullname);
            String residues = br.readLine();
            Item sequence = createItem("Sequence");
            sequence.setAttribute("residues", residues);
            sequence.setAttribute("length", String.valueOf(residues.length()));
            store(sequence);
            gm.setReference("sequence", sequence);
        }
        br.close();
    }

    /**
     * Process a genetic map GFF file.
     *
     * 0     1     2    3    4   5
     * phavu.mixed.map1.7PMp.map.gff3
     *
     * ##gff-version 3
     * ##date Sat Jul  8 11:12:08 2017
     * ##source gbrowse gbgff gff3 dumper, from https://legumeinfo.org/genomes/gbrowse/Pv1.0 with post-processing
     * phavu.G19833.gnm1.Chr01 blastn  SNP     242423  242423  .       +       .       Name=Pv_TOG905303_749;ID=1;Note=LG01 cM 0.0;alleles=T/G
     */
    void processMapGFFFile(Reader reader) throws IOException, ObjectStoreException {
        // get the identifiers from the filename
        String[] dotparts = getCurrentFile().getName().split("\\.");
        String dataSetName = getCurrentFile().getName();
        String gensp = dotparts[0];
        String parents = dotparts[1];
        String dataSetVersion = dotparts[2];
        // create Items
        Item dataSet = getDataSet(dataSetName);
        dataSet.setAttribute("version", dataSetVersion);
        Item organism = getOrganism(gensp);
        // spin through the lines
        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            MapGFFRecord record = new MapGFFRecord(line);
            if (!record.hasData()) continue;
            Item chr = getChromosomeOrSupercontig(record.chr, organism);
            chr.addToCollection("dataSets", dataSet);
            Item gm = getGeneticMarker(record.name);
            gm.setReference("organism", organism);
            gm.addToCollection("dataSets", dataSet);
            if (!record.name.equals(record.fullname)) gm.setAttribute("secondaryIdentifier", record.fullname);
            gm.setAttribute("type", record.type);
            gm.setAttribute("length", String.valueOf(record.end-record.start+1));
            if (record.alleles!=null) gm.setAttribute("alleles", record.alleles);
            Item location = createItem("Location");
            location.setAttribute("strand", record.strand);
            location.setAttribute("start", String.valueOf(record.start));
            location.setAttribute("end", String.valueOf(record.end));
            location.setReference("feature", gm);
            location.addToCollection("dataSets", dataSet);
            location.setReference("locatedOn", chr);
            store(location);
            if (chr.getClassName().equals("Supercontig")) {
                gm.setReference("supercontig", chr);
                gm.setReference("supercontigLocation", location);
            } else {
                gm.setReference("chromosome", chr);
                gm.setReference("chromosomeLocation", location);
            }
        }
        br.close();
    }

    /**
     * Get the gene descriptions and ontology annotations from an info_descriptors.txt file.
     *
     * 0     1         2    3    4    5                6
     * arahy.Tifrunner.gnm1.ann1.CCJH.info_descriptors.txt
     * arahy.Tifrunner.gnm1.ann1.GU6A2U RING finger protein 5-like [Glycine max]; IPR013083 (Zinc finger, RING...); GO:0005515 (protein binding), GO:0008270 (zinc ion binding)
     */
    void processInfoDescriptorsFile(Reader reader) throws IOException, ObjectStoreException {
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
        Item organism = getOrganism(gensp);
        Item strain = getStrain(strainName, organism);
        // spin through the file
        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            // comment line
            if (line.startsWith("#")) continue;
            // parse record and update items
            InfoDescriptorsRecord record = new InfoDescriptorsRecord(line);
            Item gene = getGene(record.identifier);
            gene.setAttribute("description", record.description);
            gene.addToCollection("dataSets", dataSet);
            // GO terms
            for (String identifier : record.go.keySet()) {
                String description = record.go.get(identifier);
                Item goTerm = getOntologyTerm(identifier);
                goTerm.setAttribute("description", description);
                goTerm.setReference("ontology", geneOntology);
                Item goAnnotation = createItem("OntologyAnnotation");
                goAnnotation.setReference("subject", gene);
                goAnnotation.setReference("ontologyTerm", goTerm);
                goAnnotation.addToCollection("dataSets", dataSet);
                store(goAnnotation);
            }
            // ON HOLD: genes don't have proteinDomains collection and we don't have gene family IDs here!
            // // interpro domains
            // for (String identifier : record.interpro.keySet()) {
            //     String description = record.interpro.get(identifier);
            //     Item proteinDomain = getProteinDomain(identifier);
            //     proteinDomain.setAttribute("description", description);
            // }
        }
    }
}
