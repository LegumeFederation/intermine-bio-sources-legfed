package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2002-2019 FlyMine
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.NoSuchElementException;

import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.apache.tools.ant.BuildException;
import org.biojava.nbio.core.exceptions.ParserException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.io.DNASequenceCreator;
import org.biojava.nbio.core.sequence.io.FastaReader;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.biojava.nbio.core.sequence.io.PlainFastaHeaderParser;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.intermine.bio.util.OrganismData;
import org.intermine.bio.util.OrganismRepository;
import org.intermine.metadata.Util;
import org.intermine.model.InterMineObject;
import org.intermine.model.bio.BioEntity;
import org.intermine.model.bio.DataSet;
import org.intermine.model.bio.DataSource;
import org.intermine.model.bio.Organism;
import org.intermine.model.bio.Strain;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.objectstore.query.PendingClob;
import org.intermine.task.FileDirectDataLoaderTask;

/**
 * A task that can read a set of FASTA files and create the corresponding Sequence objects in an
 * ObjectStore.
 *
 * This minor tweak changes the class to Supercontig if the identifier contains "scaffold" or other LIS indicators.
 *
 * @author Kim Rutherford
 * @author Peter Mclaren
 * @author Sam Hokin
 */

public class LegfedFastaLoaderTask extends FileDirectDataLoaderTask {
    private static final Logger LOG = Logger.getLogger(LegfedFastaLoaderTask.class);

    private String sequenceType = "dna";
    private String classAttribute = "primaryIdentifier";
    private Organism org;
    private Strain strain;
    private String strainName;
    private String assemblyVersion;
    private String annotationVersion; // for proteins
    private String className;
    private int storeCount = 0;
    private String dataSourceName = null;
    private String dataSourceUrl = null;
    private DataSource dataSource = null;
    private String dataSetTitle;
    private String dataSetUrl;
    private String dataSetVersion; // from file name
    private Map<String, DataSet> dataSets = new HashMap<String, DataSet>();
    private String fastaTaxonId = null;
    private String fastaGensp = null;
   
    // map gensp to taxonId
    Map<String,String> genspTaxonId = new HashMap<>();

    /**
     * Append this suffix to the identifier of the BioEnitys that are stored.
     */
    private String idSuffix = "";

    //Set this if we want to do some testing...
    private File[] files = null;

    /**
     * Set the sequence type to be passed to the FASTA parser.  The default is "dna".
     * @param sequenceType the sequence type
     */
    public void setSequenceType(String sequenceType) {
        if ("${fasta.sequenceType}".equals(sequenceType)) {
            this.sequenceType = "dna";
        } else {
            this.sequenceType = sequenceType;
        }
    }

    /**
     * Set the suffix to add to identifiers from the FASTA file when creating BioEntities.
     * @param idSuffix the suffix
     */
    public void setIdSuffix(String idSuffix) {
        this.idSuffix = idSuffix;
    }

    /**
     * The class name to use for objects created during load.  Generally this is
     * "org.intermine.model.bio.LocatedSequenceFeature" or "org.intermine.model.bio.Protein"
     * @param className the class name
     */
    public void setClassName(String className) {
        this.className = className;
    }

    /**
     * Return the class name set with setClassName().
     * @return the class name
     */
    public String getClassName() {
        return className;
    }

    /**
     * The attribute of the class created to set with the identifying field.  If not set will
     * be 'primaryIdentifier'.
     * @param classAttribute the class name
     */
    public void setClassAttribute(String classAttribute) {
        this.classAttribute = classAttribute;
    }

    /**
     * DataSource.name for any bioentities created
     * @param dataSourceName name of datasource for items created
     */
    public void setDataSourceName(String dataSourceName) {
        this.dataSourceName = dataSourceName;
    }

    /**
     * DataSource.url for any bioentities created
     * @param dataSourceUrl url of datasource for items created
     */
    public void setDataSourceUrl(String dataSourceUrl) {
        this.dataSourceUrl = dataSourceUrl;
    }

    /**
     * If a value is specified this title will used when a DataSet is created.
     * @param dataSetTitle the title of the DataSets of any new features
     */
    public void setDataSetTitle(String dataSetTitle) {
        this.dataSetTitle = dataSetTitle;
    }

    /**
     * If a value is specified this url will used when a DataSet is created.
     * @param dataSetUrl the title of the DataSets of any new features
     */
    public void setDataSetUrl(String dataSetUrl) {
        this.dataSetUrl = dataSetUrl;
    }

    /**
     * Directly set the array of files to read from.  Use this for testing with junit.
     * @param files the File objects
     */
    protected void setFileArray(File[] files) {
        this.files = files;
    }

    /**
     * Process and load all of the fasta files.
     */
    @Override
    public void process() {
        long start = System.currentTimeMillis();
        // get the organism data
        DatastoreUtils dsu = new DatastoreUtils();
        genspTaxonId = dsu.getGenspTaxonId();
        try {
            storeCount++;
            super.process();
            getIntegrationWriter().commitTransaction();
            getIntegrationWriter().beginTransaction();
            getDirectDataLoader().close();
        } catch (ObjectStoreException e) {
            throw new BuildException("failed to store object", e);
        }
        long now = System.currentTimeMillis();
        LOG.info("Finished dataloading " + storeCount + " objects at " + ((60000L * storeCount)
                                                                          / (now - start)) + " objects per minute (" + (now - start)
                 + " ms total) for source " + sourceName);
    }

    /**
     * Be sure to close the data loader so the last batch gets stored. only needed for tests
     * since the data loading task usually does that for hte live builds.
     * @throws ObjectStoreException if we can't store to db
     */
    public void close() throws ObjectStoreException {
        // store any data left over
        getDirectDataLoader().close();
    }

    /**
     * @throws BuildException if an ObjectStore method fails
     */
    @Override
    public void execute() {
        // don't configure dynamic attributes if this is a unit test!
        if (getProject() != null) {
            configureDynamicAttributes(this);
        }
        if (className == null) {
            throw new RuntimeException("className needs to be set");
        }
        if (files != null) {
            // setFiles() is used only for testing
            for (int i = 0; i < files.length; i++) {
                processFile(files[i]);
            }
        } else {
            // this will call processFile() for each file
            super.execute();
        }
    }

    /**
     * Handles each fasta file. Factored out so we can supply files for testing.
     *
     * @param file the File to process.
     * @throws BuildException if the is a problem
     */
    @Override
    public void processFile(File file) {
        try {
            System.out.println("##############################################################################################################################");
            System.out.println("Reading "+sequenceType+" sequences from: "+file);
            System.out.println("##############################################################################################################################");
            LOG.info("LegfedFastaLoaderTask loading file "+file.getName());
            // pull the organism, strain, assembly and annotation version (if present) from the file name
            // 0     1      2    3    4    5       6
            // phavu.G19833.gnm2.ann1.PB8d.protein.faa
            // 0     1      2    3    4           5
            // phavu.G19833.gnm2.fC0g.genome_main.fna
            String[] parts = file.getName().split("\\.");
            fastaGensp = parts[0];
            fastaTaxonId = genspTaxonId.get(fastaGensp);
            if (fastaTaxonId==null) {
                throw new RuntimeException("Organism "+fastaGensp+" is not in the genspTaxonId map.");
            }
            strainName = parts[1];
            assemblyVersion = parts[2];
            dataSetVersion = assemblyVersion;
            if (parts.length==7) {
                annotationVersion = parts[3];
                dataSetVersion += "."+annotationVersion;
            }
            if (sequenceType.equalsIgnoreCase("dna")) {
                FastaReader<DNASequence, NucleotideCompound> aFastaReader
                    = new FastaReader<DNASequence, NucleotideCompound>(file,
                                                                       new PlainFastaHeaderParser<DNASequence, NucleotideCompound>(),
                                                                       new DNASequenceCreator(AmbiguityDNACompoundSet.getDNACompoundSet()));
                LinkedHashMap<String, DNASequence> b = aFastaReader.process();
                for (Entry<String, DNASequence> entry : b.entrySet()) {
                    Sequence bioJavaSequence = entry.getValue();
                    processSequence(getOrganism(bioJavaSequence), getStrain(bioJavaSequence), bioJavaSequence);
                }
            } else {
                LinkedHashMap<String, ProteinSequence> b =
                    FastaReaderHelper.readFastaProteinSequence(file);
                for (Entry<String, ProteinSequence> entry : b.entrySet()) {
                    Sequence bioJavaSequence = entry.getValue();
                    processSequence(getOrganism((ProteinSequence) bioJavaSequence), getStrain(bioJavaSequence), bioJavaSequence);
                }
            }
        } catch (ParserException e) {
            throw new BuildException("Sequence not in FASTA format or wrong alphabet for: "+file, e);
        } catch (NoSuchElementException e) {
            throw new BuildException("No FASTA sequences in: "+file, e);
        } catch (FileNotFoundException e) {
            throw new BuildException("Problem reading file - file not found: "+file, e);
        } catch (ObjectStoreException e) {
            throw new BuildException("ObjectStore problem while processing: "+file, e);
        } catch (IOException e) {
            throw new BuildException("Error while closing FileReader for: "+file, e);
        }
    }

    /**
     * Get and store() the Organism object to reference when creating new objects.
     * @param bioJavaSequence the biojava sequence to be parsed
     * @throws ObjectStoreException if there is a problem
     * @return the new Organism
     */
    protected Organism getOrganism(Sequence bioJavaSequence) throws ObjectStoreException {
        if (org == null) {
            org = getDirectDataLoader().createObject(Organism.class);
            org.setTaxonId(fastaTaxonId);
            org.setGensp(fastaGensp);
            getDirectDataLoader().store(org);
        }
        return org;
    }

    /**
     * Get and store() the Strain object to reference when creating new objects.
     * @param bioJavaSequence the biojava sequence to be parsed (not used)
     * @throws ObjectStoreException if there is a problem
     * @return the new Strain
     */
    protected Strain getStrain(Sequence bioJavaSequence) throws ObjectStoreException {
        if (strain == null) {
            strain = getDirectDataLoader().createObject(Strain.class);
            strain.setPrimaryIdentifier(strainName);
            if (org!=null) strain.setOrganism(org);
            getDirectDataLoader().store(strain);
        }
        return strain;
    }

    /**
     * Create a Sequence and an object of type className for the given BioJava Sequence.
     * @param organism the Organism to reference from new objects
     * @param strain the Strain to reference from new objects (null to avoid reference)
     * @param bioJavaSequence the Sequence object
     * @throws ObjectStoreException if store() fails
     */
    private void processSequence(Organism organism, Strain strain, Sequence bioJavaSequence) throws ObjectStoreException {
        // some fasta files are not filtered - they contain sequences from organisms not
        // specified in project.xml
        if (organism==null) {
            return;
        }

        org.intermine.model.bio.Sequence bioSequence = getDirectDataLoader().createObject(org.intermine.model.bio.Sequence.class);

        String sequence = bioJavaSequence.getSequenceAsString();
        String md5checksum = Util.getMd5checksum(sequence);
        
        bioSequence.setResidues(new PendingClob(sequence));
        bioSequence.setLength(bioJavaSequence.getLength());
        bioSequence.setMd5checksum(md5checksum);

        // the identifier
        String attributeValue = getIdentifier(bioJavaSequence);

        // HACK: don't allow spaces or tabs in sequence primary identifiers; set symbol=extra part
        String symbol = null;
        String[] spaceChunks = attributeValue.split(" ");
        if (spaceChunks.length>1) {
            attributeValue = spaceChunks[0];
            symbol = spaceChunks[1];
        }
        String[] tabChunks = attributeValue.split("\t");
        if (tabChunks.length>1) {
            attributeValue = tabChunks[0];
            symbol = tabChunks[1];
        }

        // HACK: change the className to "Supercontig" if identifier contains "scaffold" or ends in "sc" etc.
        String[] dotparts = attributeValue.split("\\.");
        String lastPart = dotparts[dotparts.length-1];
        if (attributeValue.toLowerCase().contains("scaffold")
            || lastPart.contains("sc")
            || lastPart.contains("pilon")
            ) {
            className = "org.intermine.model.bio.Supercontig";
        }

        Class<? extends InterMineObject> imClass;
        Class<?> c;
        try {
            c = Class.forName(className);
            if (InterMineObject.class.isAssignableFrom(c)) {
                imClass = (Class<? extends InterMineObject>) c;
            } else {
                throw new RuntimeException("Feature className must be a valid class in the model that inherits from InterMineObject, but was: " + className);
            }
        } catch (ClassNotFoundException e1) {
            throw new RuntimeException("unknown class: " + className + " while creating new Sequence object");
        }

        // create the object that has the sequence
        BioEntity imo = (BioEntity) getDirectDataLoader().createObject(imClass);

        try {
            imo.setFieldValue(classAttribute, attributeValue);
        } catch (Exception e) {
            throw new IllegalArgumentException("Error setting: "+className+"."+classAttribute+" to: "+attributeValue+". Does the attribute exist?");
        }
        
        try {
            imo.setFieldValue("sequence", bioSequence);
        } catch (Exception e) {
            throw new IllegalArgumentException("Error setting: "+className+".sequence to: "+attributeValue+". Does the sequence attribute exist?");
        }

        imo.setOrganism(organism);
        if (strain!=null) imo.setStrain(strain);
        if (assemblyVersion!=null) imo.setAssemblyVersion(assemblyVersion);
        if (annotationVersion!=null) imo.setAnnotationVersion(annotationVersion);

        try {
            imo.setFieldValue("length", new Integer(bioSequence.getLength()));
        } catch (Exception e) {
            throw new IllegalArgumentException("Error setting: "+className+".length to: "+bioSequence.getLength()+". Does the attribute exist?");
        }

        try {
            imo.setFieldValue("md5checksum", md5checksum);
        } catch (Exception e) {
            // Ignore - we don't care if the field doesn't exist.
        }

        try {
            if (symbol!=null) imo.setFieldValue("symbol", symbol);
        } catch (Exception e) {
            throw new IllegalArgumentException("Error setting: "+className+".symbol to: "+attributeValue+". Does the symbol attribute exist?");
        }


        if (StringUtils.isEmpty(dataSetTitle)) {
            throw new RuntimeException("DataSet title (legfed-fasta.dataSetTitle) not set.");
        }

        extraProcessing(bioJavaSequence, bioSequence, imo, organism, strain, getDataSet());

        DataSet dataSet = getDataSet();
        imo.addDataSets(dataSet);
        try {
            getDirectDataLoader().store(bioSequence);
            getDirectDataLoader().store(imo);
            storeCount += 2;
        } catch (ObjectStoreException e) {
            throw new BuildException("store failed", e);
        }
    }

    /**
     * Return the DataSet to add to each object.
     * @return the DataSet
     * @throws ObjectStoreException if there is an ObjectStore problem
     */
    public DataSet getDataSet() throws ObjectStoreException {
        if (dataSets.containsKey(dataSetTitle)) {
            return dataSets.get(dataSetTitle);
        }
        DataSet dataSet = getDirectDataLoader().createObject(DataSet.class);
        dataSet.setName(dataSetTitle);
        if (dataSetUrl!=null) dataSet.setUrl(dataSetUrl);
        if (dataSetVersion!=null) dataSet.setVersion(dataSetVersion);
        if (dataSourceName!=null) dataSet.setDataSource(getDataSource());
        getDirectDataLoader().store(dataSet);
        dataSets.put(dataSetTitle, dataSet);
        return dataSet;
    }


    /**
     * Do any extra processing needed for this record (extra attributes, objects, references, etc.).
     * This method is called before the new objects are stored.
     * @param bioJavaSequence the BioJava Sequence
     * @param imSequence the IntermMine Sequence
     * @param bioEntity the object that references the sequence
     * @param organism the Organism object for the new InterMineObject
     * @param dataSet the DataSet object
     * @throws ObjectStoreException if a store() fails during processing
     */
    protected void extraProcessing(Sequence bioJavaSequence, org.intermine.model.bio.Sequence imSequence,
                                   BioEntity bioEntity, Organism organism, Strain strain, DataSet dataSet)
        throws ObjectStoreException {
        // create a secondaryIdentifier omitting the strain, assembly version and annotation version (1, 2, 3 below).
        // 0     1      2    3    4     5          6
        // phavu.G19833.gnm2.ann1.Phvul.010G034400.1
        String identifier = getIdentifier(bioJavaSequence);
        String[] parts = identifier.split("\\.");
        String secondaryIdentifier = parts[0];                   // phavu
        if (parts.length>4) secondaryIdentifier += "."+parts[4]; // phavu.Phvul
        if (parts.length>5) secondaryIdentifier += "."+parts[5]; // phavu.Phvul.010G034400
        if (parts.length>6) secondaryIdentifier += "."+parts[6]; // phavu.Phvul.010G034400.1
        if (parts.length>7) secondaryIdentifier += "."+parts[7]; // phavu.Phvul.010G034400.1.other
        if (parts.length>8) secondaryIdentifier += "."+parts[8]; // phavu.Phvul.010G034400.1.other.stuff
        if (parts.length>9) secondaryIdentifier += "."+parts[9]; // phavu.Phvul.010G034400.1.other.stuff.here
        bioEntity.setSecondaryIdentifier(secondaryIdentifier);
    }

    /**
     * For the given BioJava Sequence object, return an identifier to be used when creating
     * the corresponding BioEntity.
     * if | is present the middle bit is returned, eg sp|Q9V8R9-2|41_DROME
     * @param bioJavaSequence the Sequenece
     * @return an identifier
     */
    protected String getIdentifier(Sequence bioJavaSequence) {
        String name = bioJavaSequence.getAccession().getID() + idSuffix;
        // getID does not seem to work properly
        // quick fix to get only the primaryidentifier
        if (name.contains(" ")) {
            String[] bits = name.split(" ");
            name = bits[0];
        }
        // description_line=sp|Q9V8R9-2|41_DROME
        if (name.contains("|")) {
            String[] bits = name.split("\\|");
            if (bits.length < 2) {
                return null;
            }
            name = bits[1];
        }
        return name;
    }

    private DataSource getDataSource() throws ObjectStoreException {
        if (StringUtils.isEmpty(dataSourceName)) {
            throw new RuntimeException("dataSourceName not set");
        }
        if (dataSource == null) {
            dataSource = getDirectDataLoader().createObject(DataSource.class);
            dataSource.setName(dataSourceName);
            if (dataSourceUrl!=null) dataSource.setUrl(dataSourceUrl);
            getDirectDataLoader().store(dataSource);
            storeCount += 1;
        }
        return dataSource;
    }

    /**
     * Get and store() the Organism object to reference when creating new objects.
     * @param bioJavaSequence the biojava sequence to be parsed
     * @throws ObjectStoreException if there is a problem
     * @return the new Organism
     */
    protected Organism getOrganism(ProteinSequence bioJavaSequence)
        throws ObjectStoreException {
        if (org == null) {
            org = getDirectDataLoader().createObject(Organism.class);
            org.setTaxonId(fastaTaxonId);
            org.setGensp(fastaGensp);
            getDirectDataLoader().store(org);
        }
        return org;
    }
}

