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
import java.util.List;
import java.util.ArrayList;
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

    // store Items in maps if they may be read more than once
    Map<String,Item> organisms = new HashMap<>();
    Map<String,Item> strains = new HashMap<>();
    Map<String,Item> ontologyTerms = new HashMap<>();
    Map<String,Item> dataSets = new HashMap<>();
    Map<String,Item> geneFamilies = new HashMap<>();
    Map<String,Item> proteinDomains = new HashMap<>();
    Map<String,Item> chromosomes = new HashMap<>(); // contains both Chromosome and Supercontig Items
    Map<String,Item> genes = new HashMap<>();
    Map<String,Item> proteins = new HashMap<>();
    Map<String,Item> mRNAs = new HashMap<>();
    Map<String,Item> linkageGroups = new HashMap<>();
    Map<String,Item> geneticMarkers = new HashMap<>();
    Map<String,Item> phenotypes = new HashMap<>();
    Map<String,Item> publicationsByDoi = new HashMap<>();   // this will store the same pub twice if it's DOI-referenced in one case
    Map<Integer,Item> publicationsByPmid = new HashMap<>(); // and PMID-referenced in another case
    Map<String,Item> ontologyAnnotations = new HashMap<>(); // keyed by identifier_version_subject

    // ontologies created in constructor for random use
    Item geneOntology;
    Item pfamOntology;
    Item pantherOntology;
    Item kogOntology;
    Item ecOntology;
    Item koOntology;

    // there is only one dataSource, set in project.xml
    Item dataSource;

    // DataSource is set once in project.xml (DataSet instances are formed from file names.)
    String dataSourceName;
    String dataSourceUrl;
    String dataSourceDescription;

    // instantiate for non-static utility methods
    DatastoreUtils dsu;

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
        // KEGG Ontology
        koOntology = createItem("Ontology");
        koOntology.setAttribute("name", "KEGG Orthology");
        koOntology.setAttribute("url", "https://www.genome.jp/kegg/ko.html");
        store(koOntology);
        
        // for non-static utility methods
        dsu = new DatastoreUtils();
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

    /**
     * {@inheritDoc}
     */
    @Override
    public void process(Reader reader) throws Exception {
        if (dataSource==null) {
            // set defaults for LIS if not given
            if (dataSourceName==null) dataSourceName = DatastoreUtils.DEFAULT_DATASOURCE_NAME;
            if (dataSourceName.equals(DatastoreUtils.DEFAULT_DATASOURCE_NAME)) {
                dataSourceUrl = DatastoreUtils.DEFAULT_DATASOURCE_URL;
                dataSourceDescription = DatastoreUtils.DEFAULT_DATASOURCE_DESCRIPTION;
            }
            // create the DataSource once
            dataSource = createItem("DataSource");
            dataSource.setAttribute("name", dataSourceName);
            if (dataSourceUrl!=null) dataSource.setAttribute("url", dataSourceUrl);
            if (dataSourceDescription!=null) dataSource.setAttribute("description", dataSourceDescription);
            store(dataSource);
        }
        // process the file
        if (getCurrentFile().getName().contains("description_")) {
            // description_Phaseolus_vulgaris.yml
            processDescriptionFile(reader);
        } else if (getCurrentFile().getName().contains("strains_")) {
            // strains_Phaseolus_vulgaris.yml
            processStrainsFile(reader);
        } else if (getCurrentFile().getName().endsWith(".info_annot.txt")) {
            // phavu.G19833.gnm2.ann1.PB8d.info_annot.txt
            processInfoAnnotFile(reader);
        } else if (getCurrentFile().getName().endsWith(".info_annot_ahrd.tsv")) {
            // legume.genefam.fam1.M65K.info_annot_ahrd.tsv
            String fastaDirname = getCurrentFile().getParent()+"/"+getCurrentFile().getName().replace("info_annot_ahrd.tsv", "family_fasta");
            printInfoBlurb(fastaDirname);
            processInfoAnnotAhrdFile(reader);
        } else if (getCurrentFile().getName().endsWith(".cmap.txt")) {
            // phavu.mixed.map1.7PMp.cmap.txt
            processCmapFile(reader);
        } else if (getCurrentFile().getName().endsWith(".map.gff3")) {
            // phavu.mixed.map1.7PMp.map.gff3
            processMapGFFFile(reader);
        } else if (getCurrentFile().getName().endsWith(".flanking_seq.fna")) {
            // phavu.mixed.map1.7PMp.flanking_seq.fna
            processFlankingSeqFastaFile(reader);
        } else if (getCurrentFile().getName().endsWith(".info_descriptors.txt")) {
            // arahy.Tifrunner.gnm1.ann1.CCJH.info_descriptors.txt
            processInfoDescriptorsFile(reader);
        } else if (getCurrentFile().getName().endsWith(".gwas.txt")) {
            // COMMENTED OUT UNTIL GWAS FILE FORMAT GETS RESOLVED
            // // glyma.mixed.LBC20180525.4.gwas.txt
            // printInfoBlurb(getCurrentFile().getName());
            // processGwasFile(reader);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void close() throws Exception {
        System.out.println("Storing all Item maps...");
        store(organisms.values());
        store(strains.values());
        store(ontologyTerms.values());
        store(dataSets.values());
        store(geneFamilies.values());
        store(proteinDomains.values());
        store(chromosomes.values());
        store(genes.values());
        store(proteins.values());
        store(mRNAs.values());
        store(linkageGroups.values());
        store(geneticMarkers.values());
        store(phenotypes.values());
        store(publicationsByDoi.values());
        store(publicationsByPmid.values());
        store(ontologyAnnotations.values());
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
     * Get/add the DataSet Item, formed from input parameters and the instance dataSource.
     */
    Item getDataSet(String name, String url, String description, String version) {
        if (dataSets.containsKey(name)) {
            return dataSets.get(getCurrentFile().getName());
        } else {
            Item dataSet = createItem("DataSet");
            dataSet.setAttribute("name", getCurrentFile().getName());
            if (url!=null) dataSet.setAttribute("url", url);
            if (description!=null) dataSet.setAttribute("description", description);
            if (version!=null) dataSet.setAttribute("version", version);
            // dataSource must already be set
            dataSet.setReference("dataSource", dataSource);
            dataSets.put(getCurrentFile().getName(), dataSet);
            return dataSet;
        }
    }
    
    /**
     * Get/add the DataSet Item, with attributes formed from the current filename.
     * getName()   glyma.Wm82.gnm2.ann1.RVB6.info_annot.txt
     * getParent() /data/public/Glycine_max/Wm82.gnm2.ann1.RVB6/
     * https://legumeinfo.org/data/public/Glycine_max/Wm82.gnm2.ann1.RVB6/glyma.Wm82.gnm2.ann1.RVB6.info_annot.txt.gz
     */
    Item getDataSet() {
        String name = getCurrentFile().getName();     
        String parent = getCurrentFile().getParent();
        String url = "https://legumeinfo.org"+parent+name;
        String description = "LIS datastore file.";
        return getDataSet(name, url, description, null);
    }

    /**
     * Get/add the DataSet Item, with version attribute supplied.
     */
    Item getDataSet(String version) {
        String name = getCurrentFile().getName();
        String parent = getCurrentFile().getParent();
        String url = "https://legumeinfo.org"+parent+"/"+name+".gz";
        String description = "LIS datastore file.";
        return getDataSet(name, url, description, version);
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
            String taxonId = dsu.getTaxonId(gensp);
            String genus = dsu.getGenus(gensp);
            String species = dsu.getSpecies(gensp);
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
    Item getStrain(String strainId, Item organism) {
        Item strain;
        if (strains.containsKey(strainId)) {
            strain = strains.get(strainId);
        } else {
            strain = createItem("Strain");
            strain.setAttribute("identifier", strainId);
            strain.setReference("organism", organism);
            strains.put(strainId, strain);
        }
        return strain;
    }

    /**
     * Get/add a Gene Item, keyed by primaryIdentifier
     */
    public Item getGene(String primaryIdentifier) {
        Item gene;
        if (genes.containsKey(primaryIdentifier)) {
            gene = genes.get(primaryIdentifier);
        } else {
            // phavu.Phvul.002G040500
            gene = createItem("Gene");
            gene.setAttribute("primaryIdentifier", primaryIdentifier);
            genes.put(primaryIdentifier, gene);
        }
        return gene;
    }

    /**
     * Get/add a Protein Item, keyed by primaryIdentifier
     */
    public Item getProtein(String primaryIdentifier) {
        Item protein;
        if (proteins.containsKey(primaryIdentifier)) {
            protein = proteins.get(primaryIdentifier);
        } else {
            protein = createItem("Protein");
            protein.setAttribute("primaryIdentifier", primaryIdentifier);
            proteins.put(primaryIdentifier, protein);
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
     * Get/add an MRNA Item, keyed by primaryIdentifier (!)
     */
    public Item getMRNA(String primaryIdentifier) {
        Item mRNA;
        if (mRNAs.containsKey(primaryIdentifier)) {
            mRNA = mRNAs.get(primaryIdentifier);
        } else {
            // Phvul.002G040500.1
            mRNA = createItem("MRNA");
            mRNA.setAttribute("primaryIdentifier", primaryIdentifier);
            mRNAs.put(primaryIdentifier, mRNA);
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
     * Get/add a GeneticMarker Item based on secondaryIdentifier match (since it's the same marker regardless of assembly/annotation in this loader).
     * Plus a lot of the files that this loader loads do not have "full yuck" identifiers (e.g. SoyBase GWAS).
     */
    public Item getGeneticMarker(String secondaryIdentifier) {
        if (geneticMarkers.containsKey(secondaryIdentifier)) {
            return geneticMarkers.get(secondaryIdentifier);
        } else {
            Item gm = createItem("GeneticMarker");
            gm.setAttribute("secondaryIdentifier", secondaryIdentifier);
            geneticMarkers.put(secondaryIdentifier, gm);
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
            } else {
                chr = createItem("Chromosome");
            }
            chr.setAttribute("primaryIdentifier", primaryIdentifier);
            chr.setReference("organism", organism);
            chromosomes.put(primaryIdentifier, chr);
            return chr;
        }
    }

    /**
     * Get/add a Publication Item based on DOI.
     */
    public Item getPublication(String doi) {
        if (publicationsByDoi.containsKey(doi)) {
            return publicationsByDoi.get(doi);
        } else {
            Item publication = createItem("Publication");
            publication.setAttribute("doi", doi);
            publicationsByDoi.put(doi, publication);
            return publication;
        }
    }
        
    /**
     * Get/add a Publication Item based on DOI.
     */
    public Item getPublication(int pmid) {
        if (publicationsByPmid.containsKey(pmid)) {
            return publicationsByPmid.get(pmid);
        } else {
            Item publication = createItem("Publication");
            publication.setAttribute("pubMedId", String.valueOf(pmid));
            publicationsByPmid.put(pmid, publication);
            return publication;
        }
    }

    /**
     * Get/add a Phenotype base on primary identifier.
     */
    public Item getPhenotype(String name) {
        if (phenotypes.containsKey(name)) {
            return phenotypes.get(name);
        } else {
            Item phenotype = createItem("Phenotype");
            phenotype.setAttribute("name", name);
            phenotypes.put(name, phenotype);
            return phenotype;
        }
    }
    
    /**
     * Form a key for an OntologyAnnotation for dupe avoidance.
     */
    String formAnnotKey(String termIdentifier, String versionKey, String subjectId) { 
        return termIdentifier+"_"+versionKey+"_"+subjectId;
    }

    /**
     * Form a primaryIdentifier from the core identifier, gensp, strain, assembly and annotation
     */
    String formPrimaryIdentifier(String gensp, String strain, String assembly, String annotation, String id) {
        return gensp+"."+strain+"."+assembly+"."+annotation+"."+id;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////// FILE PROCESSORS ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
        Map<String,String> attributes = new HashMap<>();
        String line = null;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#####")) {
                // new strain section, store previous strain
                if (attributes.size()>0) {
                    String strainId = attributes.get("identifier");
                    Item strain = getStrain(strainId, organism);
                    for (String name : attributes.keySet()) {
                        String value = attributes.get(name);
                        strain.setAttribute(name, value);
                    }
                }
                // clear attributes map
                attributes = new HashMap<>();
            } else if (line.startsWith("#") || line.startsWith("%")) {
                // comment
                continue;
            } else {
                // put strain attribute into map
                String[] parts = line.split("\t");
                if (parts.length>1) {
                    String attributeName = parts[0].replace("strain.","").replace(":","");
                    String attributeValue = parts[1].trim();
                    attributes.put(attributeName, attributeValue);
                }
            }
        }
        br.close();
        // last one
        if (attributes.size()>0) {
            String strainId = attributes.get("identifier");
            Item strain = getStrain(strainId, organism);
            for (String name : attributes.keySet()) {
                String value = attributes.get(name);
                strain.setAttribute(name, value);
            }
        }
    }

    /**
     * Process an info_annot.txt file which contains relationships between genes, transcripts, proteins and ontology terms.
     * This will also link genes to proteins. The file name starts with gensp.strain.assembly.annotation.
     * 0     1      2    3    4    5          6
     * phavu.G19833.gnm2.ann1.PB8d.info_annot.txt
     *
     * Note: gensp.strain.assembly.annotation must be prepended to names in this file, which only contain the core IDs.
     *
     * pacId locusName transcriptName peptideName Pfam Panther KOG ec KO GO Best-hit-arabi-name arabi-symbol arabi-defline
     * 37170591 Phvul.001G000400 Phvul.001G000400.1 Phvul.001G000400.1.p PF00504 PTHR21649,PTHR21649:SF24 1.10.3.9 K14172 GO:0016020,GO:0009765 AT1G76570.1 Chlorophyll family protein
     */
    void processInfoAnnotFile(Reader reader) throws IOException, ObjectStoreException {
        // get the dot-separated file information
        String[] dotparts = getCurrentFile().getName().split("\\.");
        String gensp = dotparts[0];
        String strainId = dotparts[1];
        String assemblyVersion = dotparts[2];
        String annotationVersion = dotparts[3];
        String versionKey = assemblyVersion+"."+annotationVersion;
        // get the dataSet, organism and strain
        Item dataSet = getDataSet(versionKey);
        Item organism = getOrganism(gensp);
        Item strain = getStrain(strainId, organism);
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
                String geneIdentifier = formPrimaryIdentifier(gensp, strainId, assemblyVersion, annotationVersion, record.locusName);
                Item gene = getGene(geneIdentifier);
                gene.setAttribute("assemblyVersion", assemblyVersion);
                gene.setAttribute("annotationVersion", annotationVersion);
                gene.setReference("organism", organism);
                gene.setReference("strain", strain);
                gene.addToCollection("dataSets", dataSet);
                // the protein
                String proteinIdentifier = formPrimaryIdentifier(gensp, strainId, assemblyVersion, annotationVersion, record.peptideName);
                Item protein = getProtein(proteinIdentifier);
                protein.setAttribute("assemblyVersion", assemblyVersion);
                protein.setAttribute("annotationVersion", annotationVersion);
                protein.setReference("organism", organism);
                protein.setReference("strain", strain);
                protein.addToCollection("genes", gene);
                protein.addToCollection("dataSets", dataSet);
                // the transcript = mRNA
                String mRNAIdentifier = formPrimaryIdentifier(gensp, strainId, assemblyVersion, annotationVersion, record.transcriptName);
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
                    String annotKey = formAnnotKey(identifier, versionKey, geneIdentifier);
                    if (!ontologyAnnotations.containsKey(annotKey)) {
                        Item goAnnotation = createItem("OntologyAnnotation");
                        goAnnotation.setReference("subject", gene);
                        goAnnotation.setReference("ontologyTerm", goTerm);
                        goAnnotation.addToCollection("dataSets", dataSet);
                        ontologyAnnotations.put(annotKey, goAnnotation);
                    }
                }
                // Pfam terms
                for (String identifier : record.Pfam) {
                    Item pfamTerm = getOntologyTerm(identifier);
                    pfamTerm.setReference("ontology", pfamOntology);
                    String annotKey = formAnnotKey(identifier, versionKey, proteinIdentifier);
                    if (!ontologyAnnotations.containsKey(annotKey)) {
                        Item pfamAnnotation = createItem("OntologyAnnotation");
                        pfamAnnotation.setReference("subject", protein);
                        pfamAnnotation.setReference("ontologyTerm", pfamTerm);
                        pfamAnnotation.addToCollection("dataSets", dataSet);
                        ontologyAnnotations.put(annotKey, pfamAnnotation);
                    }
                }
                // Panther terms
                for (String identifier : record.Panther) {
                    Item pantherTerm = getOntologyTerm(identifier);
                    pantherTerm.setReference("ontology", pantherOntology);
                    String annotKey = formAnnotKey(identifier, versionKey, proteinIdentifier);
                    if (!ontologyAnnotations.containsKey(annotKey)) {
                        Item pantherAnnotation = createItem("OntologyAnnotation");
                        pantherAnnotation.setReference("subject", protein);
                        pantherAnnotation.setReference("ontologyTerm", pantherTerm);
                        pantherAnnotation.addToCollection("dataSets", dataSet);
                        ontologyAnnotations.put(annotKey, pantherAnnotation);
                    }                }
                // KOG terms
                for (String identifier : record.KOG) {
                    Item kogTerm = getOntologyTerm(identifier);
                    kogTerm.setReference("ontology", kogOntology);
                    String annotKey = formAnnotKey(identifier, versionKey, proteinIdentifier);
                    if (!ontologyAnnotations.containsKey(annotKey)) {
                        Item kogAnnotation = createItem("OntologyAnnotation");
                        kogAnnotation.setReference("subject", protein);
                        kogAnnotation.setReference("ontologyTerm", kogTerm);
                        kogAnnotation.addToCollection("dataSets", dataSet);
                        ontologyAnnotations.put(annotKey, kogAnnotation);
                    }
                }
                // ec terms
                for (String identifier : record.ec) {
                    Item ecTerm = getOntologyTerm(identifier);
                    ecTerm.setReference("ontology", ecOntology);
                    String annotKey = formAnnotKey(identifier, versionKey, proteinIdentifier);
                    if (!ontologyAnnotations.containsKey(annotKey)) {
                        Item ecAnnotation = createItem("OntologyAnnotation");
                        ecAnnotation.setReference("subject", protein);
                        ecAnnotation.setReference("ontologyTerm", ecTerm); 
                        ecAnnotation.addToCollection("dataSets", dataSet);
                        ontologyAnnotations.put(annotKey, ecAnnotation);
                    }
                }
                // KO terms
                for (String identifier : record.KO) {
                    Item koTerm = getOntologyTerm(identifier);
                    koTerm.setReference("ontology", koOntology);
                    String annotKey = formAnnotKey(identifier, versionKey, geneIdentifier);
                    if (!ontologyAnnotations.containsKey(annotKey)) {
                        Item koAnnotation = createItem("OntologyAnnotation");
                        koAnnotation.setReference("subject", gene);
                        koAnnotation.setReference("ontologyTerm", koTerm);
                        koAnnotation.addToCollection("dataSets", dataSet);
                        ontologyAnnotations.put(annotKey, koAnnotation);
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
        String dataSetVersion = dotparts[2];
        Item dataSet = getDataSet(dataSetVersion);
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
            }
            // interpro protein domains
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
     * Process a genetic map markers file.
     *
     * phavu.mixed.map1.7PMp.markers.txt
     * LG	                ID	                                NAME	        ALIAS	        START	STOP	TYPE	LANDMARK
     * PvCookUCDavis2009_Pv01	PvCookUCDavis2009_Pv01_Pv_TOG894002	Pv_TOG894002	TOG894002	55.6	55.6	SNP	0
     */
    void processMarkersFile(Reader reader) throws IOException, ObjectStoreException {
        // get the identifiers from the filename
        String[] dotparts = getCurrentFile().getName().split("\\.");
        String gensp = dotparts[0];
        String parents = dotparts[1];
        String dataSetVersion = dotparts[2];
        // create Items
        Item dataSet = getDataSet(dataSetVersion);
        Item organism = getOrganism(gensp);
        // spin through the lines
        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            // CODE HERE!
        }
        br.close();
    }

    /**
     * Process a genetic map linkage groups file.
     *
     * phavu.mixed.map1.7PMp.lg.txt 
     * ID	                NAME	LENGTH
     * PvCookUCDavis2009_Pv01	Pv01	83.5
     */
    void processLgFile(Reader reader) throws IOException, ObjectStoreException {
        // get the identifiers from the filename
        String[] dotparts = getCurrentFile().getName().split("\\.");
        String gensp = dotparts[0];
        String parents = dotparts[1];
        String dataSetVersion = dotparts[2];
        // create Items
        Item dataSet = getDataSet(dataSetVersion);
        Item organism = getOrganism(gensp);
        // spin through the lines
        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            // CODE HERE!
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
        String gensp = dotparts[0];
        String parents = dotparts[1];
        String dataSetVersion = dotparts[2];
        // create Items
        Item dataSet = getDataSet(dataSetVersion);
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
            Item gm = getGeneticMarker(fullname);
            gm.setAttribute("primaryIdentifier", name);
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
        String gensp = dotparts[0];
        String parents = dotparts[1];
        String dataSetVersion = dotparts[2];
        // create Items
        Item dataSet = getDataSet(dataSetVersion);
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
            gm.setAttribute("primaryIdentifier", record.fullname);
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
        String strainId = dotparts[1];
        String assemblyVersion = dotparts[2];
        String annotationVersion = dotparts[3];
        String versionKey = assemblyVersion+"."+annotationVersion;
        // get the dataSet, organism and strain
        Item dataSet = getDataSet(versionKey);
        Item organism = getOrganism(gensp);
        Item strain = getStrain(strainId, organism);
        // spin through the file
        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            // comment line
            if (line.startsWith("#")) continue;
            // parse record and update items
            InfoDescriptorsRecord record = new InfoDescriptorsRecord(line);
            String geneId = record.identifier;
            Item gene = getGene(geneId);
            gene.setAttribute("description", record.description);
            gene.addToCollection("dataSets", dataSet);
            // GO terms
            for (String identifier : record.go.keySet()) {
                String annotKey = formAnnotKey(identifier, "0", geneId);
                if (!ontologyAnnotations.containsKey(annotKey)) {
                    String description = record.go.get(identifier);
                    Item goTerm = getOntologyTerm(identifier);
                    goTerm.setAttribute("description", description);
                    goTerm.setReference("ontology", geneOntology);
                    Item goAnnotation = createItem("OntologyAnnotation");
                    goAnnotation.setReference("subject", gene);
                    goAnnotation.setReference("ontologyTerm", goTerm);
                    goAnnotation.addToCollection("dataSets", dataSet);
                    ontologyAnnotations.put(annotKey, goAnnotation);
                }
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

    /**
     * Get the GWAS experiment and associated marker-phenotype relations from a gwas.txt file.
     * 0     1     2           3 4    5
     * glyma.mixed.KGK20170707.1.gwas.txt
     *
     * [ ] GWAS 
     *   [ ] Number Germplasm Tested Integer 
     *   [ ] Number Loci Tested Integer 
     *   [ ] Platform Details
     *   [ ] Platform Name
     *   [ ] Primary Identifier
     *   [+] Ontology Annotations Ontology Annotation 
     *   [+] Publications Publication 
     *   [+] Results GWAS Result  
     *
     * #Name=KGK20170707.1
     * #PlatformName=Illumina GoldenGate
     * #NumberLociTested=1142
     * #NumberGermplasmTested=219
     * #Assembly=Wm82.a2
     * #DOI=10.1111/pbr.12305
     * CHR                  BP      MARKER            PVAL    BPEND   PHENOTYPE           ONTOLOGY_IDENTIFIER
     * glyma.Wm82.gnm2.Gm02 2380973 glyma.ss107912620 3.27E-7 2380973 Stem diameter, main SOY:0001624
     */
    void processGwasFile(Reader reader) throws IOException, ObjectStoreException {
        // get the dot-separated file information
        String[] dotparts = getCurrentFile().getName().split("\\.");
        String gensp = dotparts[0];
        String populations = dotparts[1];
        String experiment = dotparts[2]+"."+dotparts[3];
        // get the dataSet
        Item dataSet = getDataSet();
        // create this GWAS
        Item gwas = createItem("GWAS");
        gwas.setReference("dataSet", dataSet);
        // spin through the file
        boolean pubFound = false;
        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            // metadata
            if (line.startsWith("#Name")) {
                String[] parts = line.split("=");
                gwas.setAttribute("primaryIdentifier", parts[1]);
                continue;
            } else if (line.startsWith("#PlatformName")) {
                String[] parts = line.split("=");
                gwas.setAttribute("platformName", parts[1]);
                continue;
            } else if (line.startsWith("#PlatformDetails")) {
                String[] parts = line.split("=");
                gwas.setAttribute("platformDetails", parts[1]);
                continue;
            } else if (line.startsWith("#NumberLociTested")) {
                String[] parts = line.split("=");
                gwas.setAttribute("numberLociTested", parts[1]);
                continue;
            } else if (line.startsWith("#NumberGermplasmTested")) {
                String[] parts = line.split("=");
                gwas.setAttribute("numberGermplasmTested", parts[1]);
                continue;
            } else if (line.startsWith("#Assembly")) {
                // ignore, we don't pull genomic info from GWAS files
                continue;
            } else if (line.startsWith("#DOI") && !pubFound) {
                String[] parts = line.split("=");
                String doi = parts[1];
                Item publication = getPublication(doi);
                gwas.addToCollection("publications", publication);
                pubFound = true;
                continue;
            } else if (line.startsWith("#PMID") && !pubFound) {
                String[] parts = line.split("=");
                int pmid = Integer.parseInt(parts[1]);
                Item publication = getPublication(pmid);
                gwas.addToCollection("publications", publication);
                pubFound = true;
                continue;
            } else if (line.startsWith("#")) {
                // an actual comment
                continue;
            } else if (line.startsWith("CHR")) {
                // the column header
                continue;
            }
            // parse record and update items; we don't use the location info here, since we load marker locations elsewhere
            // marker, pval, phenotype, ontologyIdentifier
            GwasRecord record = new GwasRecord(line);
            Item gwasResult = createItem("GWASResult");
            // assume GWAS file has marker = secondaryIdentifier e.g. glyma.ss107919771
            Item marker = getGeneticMarker(record.marker);
            gwasResult.setReference("marker", marker);
            if (record.phenotype!=null) {
                Item phenotype = getPhenotype(record.phenotype);
                gwasResult.setReference("phenotype", phenotype);
                if (record.ontologyIdentifier!=null) {
                    String annotKey = formAnnotKey(record.ontologyIdentifier, "0", record.phenotype);
                    if (!ontologyAnnotations.containsKey(annotKey)) {
                        Item ontologyAnnotation = createItem("OntologyAnnotation");
                        Item ontologyTerm = getOntologyTerm(record.ontologyIdentifier);
                        ontologyAnnotation.setReference("ontologyTerm", ontologyTerm);
                        ontologyAnnotation.setReference("subject", phenotype);
                        ontologyAnnotations.put(annotKey, ontologyAnnotation);
                    }
                }
            }
            // not all records have a pval entry, which comes across as 0 if missing
            if (record.pval>0) {
                gwasResult.setAttribute("pValue", String.valueOf(record.pval));
            }
            gwasResult.setReference("study", gwas);
            // store the GWASResult here since it's a unique load per GWAS file
            store(gwasResult);
        }
        // store the GWAS here since it's a unique load per GWAS file
        store(gwas);
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
}
