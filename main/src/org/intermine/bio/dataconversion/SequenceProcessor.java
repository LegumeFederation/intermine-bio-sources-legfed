package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2002-2015 FlyMine
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */

import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Map.Entry;

import org.apache.commons.collections.keyvalue.MultiKey;
import org.apache.commons.collections.map.MultiKeyMap;
import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.intermine.bio.chado.config.ConfigAction;
import org.intermine.bio.chado.config.CreateSynonymAction;
import org.intermine.bio.chado.config.DoNothingAction;
import org.intermine.bio.chado.config.SetFieldConfigAction;
import org.intermine.bio.util.OrganismData;
import org.intermine.metadata.ClassDescriptor;
import org.intermine.metadata.FieldDescriptor;
import org.intermine.metadata.MetaDataException;
import org.intermine.metadata.ReferenceDescriptor;
import org.intermine.model.bio.SequenceFeature;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.metadata.StringUtil;
import org.intermine.metadata.TypeUtil;
import org.intermine.xml.full.Item;
import org.intermine.xml.full.Reference;
import org.intermine.xml.full.ReferenceList;

/**
 * A spin of Kim Rutherford's SequenceProcessor to support the fields present in the Legfed chado data.
 * Note that the Legfed chado does NOT have proteins; rather, it has polypeptides with the following relations in feature_relationship:
 * <pre>
 * polypeptide derives_from mRNA
 * mRNA part_of gene
 * </pre>
 * As a result, we also tie polypeptides to genes ("producedby") via the intermediate mRNA as did Kim for proteins and genes.
 *
 * This processor is only for genomic quantities. Genetic quantites such as QTLs and linkage groups are dealt with in GeneticProcessor, as well as 
 * genetic relationships like genetic marker positions on QTLs (Nearest Marker, Flanking Marker Low, Flanking Marker High).
 *
 * @author Sam Hokin
 */
public class SequenceProcessor extends ChadoProcessor {
	
    private static final Logger LOG = Logger.getLogger(SequenceProcessor.class);

    // set true for extra verbose debug lines
    private static final boolean DEBUG = true;
    
    // a map from chado feature id to FeatureData objects, populated by processFeatureTable()
    // and used to get object types, Item IDs etc. (see FeatureData)
    protected Map<Integer, FeatureData> featureMap = new HashMap<Integer, FeatureData>();

    // we don't configure anything by default, so the process methods do their default actions
    private static final MultiKeyMap DEFAULT_CONFIG = new MultiKeyMap();

    // A map from sourceFeature uniqueName to chado feature_ids, populated by processFeatureTable()
    private Map<Integer, Map<String, Integer>> sourceFeatureMaps = new HashMap<Integer, Map<String, Integer>>();

    // a map from chado pubmed id to item identifier for the publication
    private Map<Integer, String> publications = new HashMap<Integer, String>();

    // the name of the temporary table we create from the feature table to speed up processing
    private String tempFeatureTableName = null;

    // a list of all the cvterm relation names we wish to process, including those generated locally
    // NOTE: "producedby" is generated locally to tie polypeptides to genes
    private static final List<String> DESIRED_RELATIONS = Arrays.asList("producedby", "part_of", "derives_from");
	
    // a list of the one-way relation cvterm names; these will only be processed in the subject-object direction
    private static final List<String> PARTOF_RELATIONS = Arrays.asList("producedby", "part_of", "derives_from");

    // desired feature cvterm types to query from the chado.feature table INCLUDING those that appear as sources in chado.featureloc
    private static final List<String> DESIRED_FEATURES = Arrays.asList(
								       "gene", "mRNA", "exon", 
								       "protein_hmm_match", "protein_match", 
								       "chromosome", "polypeptide", "polypeptide_domain",
								       "supercontig", "consensus_region", "genetic_marker"
								       );

    // source feature cvterm types, such as chromosome
    private static final List<String> SRC_FEATURE_TYPES = Arrays.asList( "chromosome", "polypeptide", "polypeptide_domain", "supercontig", "consensus_region" );
    
    // the InterMine info for mapping features to source feature locations; all must be IN SAME ORDER (I know, I should use maps)
    private static final String[] SRC_FEATURE_IM_CLASSES =    { "Chromosome",         "Polypeptide",         "PolypeptideDomain" ,        "Supercontig",         "ConsensusRegion" };
    private static final String[] SRC_FEATURE_IM_REFERENCES = { "chromosome",         "polypeptide",         "polypeptideDomain" ,        "supercontig",         "consensusRegion" };
    private static final String[] SRC_FEATURE_IM_LOCATIONS =  { "chromosomeLocation", "polypeptideLocation", "polypeptideDomainLocation", "supercontigLocation", "consensusRegionLocation" }; 
	
    // Avoid explosion of log messages by only logging missing collections once
    private Set<String> loggedMissingCols = new HashSet<String>();

    // Count specific debug lines
    private HashMap<String, Integer> debugLineMap = new HashMap<String, Integer>();

    /**
     * An action that makes a synonym.
     */
    protected static final ConfigAction CREATE_SYNONYM_ACTION = new CreateSynonymAction();

    /**
     * An action that does nothing - used to ignore a synonym/dbxref/whatever instead of doing the
     * default.
     */
    protected static final ConfigAction DO_NOTHING_ACTION = new DoNothingAction();

    // the prefix to use when making a temporary table, unix epoch is appended to make unique
    private static final String TEMP_FEATURE_TABLE_NAME_PREFIX = "intermine_chado_features_temp";

    static final String PRIMARY_IDENTIFIER_STRING = "primaryIdentifier";
    static final String SECONDARY_IDENTIFIER_STRING = "secondaryIdentifier";
    static final String SYMBOL_STRING = "symbol";
    static final String NAME_STRING = "name";
    static final String SEQUENCE_STRING = "sequence";
    static final String LENGTH_STRING = "length";
    static final String SOURCE_STRING = "source";

    /**
     * Create a new SequenceProcessor
     * @param chadoDBConverter the ChadoDBConverter that is controlling this processor
     */
    public SequenceProcessor(ChadoDBConverter chadoDBConverter) {
        super(chadoDBConverter);
        synchronized (this) {
            tempFeatureTableName  = TEMP_FEATURE_TABLE_NAME_PREFIX+"_"+System.currentTimeMillis()/1000L;
        }
    }

    /**
     * Initialise FeatureMap with features that have already been processed.  This is optional, it
     * permits subclasses to avoid processing the same features in multiple runs.
     * @param initialMap map of chado feature id to FeatureData objects
     */
    protected void initialiseFeatureMap(Map<Integer, FeatureData> initialMap) {
        featureMap.putAll(initialMap);
    }

    /**
     * Return the config Map.
     * @param taxonId return the configuration for this organism
     * @return the Map from configuration key to a list of actions
     */
    @SuppressWarnings("unchecked")
    protected Map<MultiKey, List<ConfigAction>> getConfig(int taxonId) {
        return DEFAULT_CONFIG;
    }

    /**
     * {@inheritDoc}
     * We process the chado database by reading each table in turn (feature, pub, featureloc, etc.)
     * Each row of each table is read and stored if appropriate.
     */
    @Override
    public void process(Connection connection) throws Exception {
        // overridden by subclasses if necessary
        earlyExtraProcessing(connection); if (DEBUG) logDebugLines();
		
        createFeatureTempTable(connection); if (DEBUG) logDebugLines();

        processFeatureTable(connection); if (DEBUG) logDebugLines();
        processFeatureCVTermTable(connection); if (DEBUG) logDebugLines();
        processPubTable(connection); if (DEBUG) logDebugLines();
        processLocationTable(connection); if (DEBUG) logDebugLines();
        processRelationTable(connection, true); if (DEBUG) logDebugLines();
        processRelationTable(connection, false); if (DEBUG) logDebugLines();
        processDbxrefTable(connection); if (DEBUG) logDebugLines();
        processSynonymTable(connection); if (DEBUG) logDebugLines();
        processFeaturePropTable(connection); if (DEBUG) logDebugLines();

        // overridden by subclasses if necessary
        extraProcessing(connection, featureMap); if (DEBUG) logDebugLines();
        // overridden by subclasses if necessary
        finishedProcessing(connection, featureMap); if (DEBUG) logDebugLines();
    }

    /**
     * Query the feature table and store features as object of the appropriate type in the
     * object store.
     * @param connection
     * @throws SQLException
     * @throws ObjectStoreException
     */
    private void processFeatureTable(Connection connection) throws SQLException, ObjectStoreException {
	Set<String> sourceFeatureTypesSet = new HashSet<String>(SRC_FEATURE_TYPES);
        ResultSet res = getFeatureTableResultSet(connection);
        int count = 0;
        while (res.next()) {
            // load vars from the ResultSet
            Integer featureId = new Integer(res.getInt("feature_id"));
            String name = res.getString("name");
            String uniqueName = res.getString("uniquename");
            String type = res.getString("type");
            String residues = res.getString("residues");
            String checksum = res.getString("md5checksum");
            Integer organismId = new Integer(res.getInt("organism_id"));
            if (sourceFeatureTypesSet.contains(type)) {
                addToSourceFeatureMaps(organismId, uniqueName, featureId);
            }
            int seqlen = res.getInt("seqlen"); // zero if missing
            // process and store this feature
            if (processAndStoreFeature(featureId, uniqueName, name, seqlen, residues, checksum, type, organismId)) {
		count++;
	    }
        }
        LOG.info("created " + count + " features");
        res.close();
    }


    /**
     * Add the given sourceFeature feature_id, uniqueName and organismId to sourceFeatureMaps.
     */
    private void addToSourceFeatureMaps(Integer organismId, String chrUniqueName, Integer chrId) {
        Map<String, Integer> sourceFeatureMap;
        if (sourceFeatureMaps.containsKey(organismId)) {
            sourceFeatureMap = sourceFeatureMaps.get(organismId);
        } else {
            sourceFeatureMap = new HashMap<String, Integer>();
            sourceFeatureMaps.put(organismId, sourceFeatureMap);
        }
        sourceFeatureMap.put(chrUniqueName, chrId);
    }

    /**
     * Create and store a new InterMineObject given data from a row of the feature table in a
     * Chado database.
     * @param featureId the chado id from the feature table
     * @param uniqueName the uniquename from Chado
     * @param name the name from Chado
     * @param seqlen the sequence length from Chado
     * @param residues the residues from Chado
     * @param md5checksum the MD5 checksum of the residues
     * @param chadoType the type of the feature from the feature + cvterm tables
     * @param organismId the chado organism id
     * @throws ObjectStoreException if there is a problem while storing
     */
    private boolean processAndStoreFeature(Integer featureId, String uniqueName, String name, int seqlen, String residues, String md5checksum, String chadoType, Integer organismId) throws ObjectStoreException {
		
        if (featureMap.containsKey(featureId)) {
            return false;
        }

        FeatureData fdat = makeFeatureData(featureId.intValue(), chadoType, uniqueName, name, md5checksum, seqlen, organismId.intValue());

        if (fdat == null) {
            return false;
        }

        String fixedUniqueName = fixIdentifier(fdat, uniqueName);

        if (seqlen > 0) {
            setAttributeIfNotSet(fdat, "length", String.valueOf(seqlen));
        }
		
        ChadoDBConverter chadoDBConverter = getChadoDBConverter();
        String dataSourceName = chadoDBConverter.getDataSourceName();
        MultiKey nameKey = new MultiKey("feature", fdat.getInterMineType(), dataSourceName, "name");
        OrganismData orgData = fdat.getOrganismData();
        List<ConfigAction> nameActionList = getConfig(orgData.getTaxonId()).get(nameKey);

        Set<String> fieldValuesSet = new HashSet<String>();
        String fixedName = fixIdentifier(fdat, name);

        // using the configuration, set a field to be the feature name
        if (!StringUtils.isBlank(fixedName)) {
            if (nameActionList == null || nameActionList.size() == 0) {
                fieldValuesSet.add(fixedName);
                setAttributeIfNotSet(fdat, SECONDARY_IDENTIFIER_STRING, fixedName);
            } else {
                for (ConfigAction action: nameActionList) {
                    if (action instanceof SetFieldConfigAction) {
                        SetFieldConfigAction attrAction =
			    (SetFieldConfigAction) action;
                        if (attrAction.isValidValue(fixedName)) {
                            String newFieldValue = attrAction.processValue(fixedName);
                            setAttributeIfNotSet(fdat, attrAction.getFieldName(), newFieldValue);
                            fieldValuesSet.add(newFieldValue);
                        }
                    }
                }
            }
        }

        MultiKey uniqueNameKey = new MultiKey("feature", fdat.getInterMineType(), dataSourceName, "uniquename");
        List<ConfigAction> uniqueNameActionList = getConfig(fdat.getOrganismData().getTaxonId()).get(uniqueNameKey);
        if (uniqueNameActionList == null || uniqueNameActionList.size() == 0) {
            // default: set primaryIdentifier to be the uniquename
            setAttributeIfNotSet(fdat, "primaryIdentifier", fixedUniqueName);
            fieldValuesSet.add(fixedUniqueName);
        } else {
            // using the configuration, set a field to be the feature name
            for (ConfigAction action: uniqueNameActionList) {
                if (action instanceof SetFieldConfigAction) {
                    SetFieldConfigAction attrAction = (SetFieldConfigAction) action;
                    if (attrAction.isValidValue(fixedUniqueName)) {
                        String newFieldValue =
			    attrAction.processValue(fixedUniqueName);
                        setAttributeIfNotSet(fdat, attrAction.getFieldName(), newFieldValue);
                        fieldValuesSet.add(newFieldValue);
                    }
                }
            }
        }

        // set the sequence if there is one
        if (fdat.checkField(SEQUENCE_STRING) && residues!=null && residues.length()>0) {
            if (!fdat.getFlag(SEQUENCE_STRING)) {
                Item sequence = getChadoDBConverter().createItem("Sequence");
                sequence.setAttribute("residues", residues);
                sequence.setAttribute("length", String.valueOf(seqlen));
                Reference chrReference = new Reference();
                chrReference.setName(SEQUENCE_STRING);
                chrReference.setRefId(sequence.getIdentifier());
                getChadoDBConverter().store(chrReference, fdat.getIntermineObjectId());
                getChadoDBConverter().store(sequence);
                fdat.setFlag(SEQUENCE_STRING, true);
            }
        }

        // create a synonym for name, if configured
        if (!StringUtils.isBlank(name)) {
            if (nameActionList != null) {
                for (ConfigAction action : nameActionList) {
                    if (action instanceof CreateSynonymAction) {
                        CreateSynonymAction createSynonymAction = (CreateSynonymAction) action;
                        if (createSynonymAction.isValidValue(fixedName)) {
                            String processedName = createSynonymAction.processValue(fixedName);
                            if (!fdat.getExistingSynonyms().contains(processedName)) {
                                Item nameSynonym = createSynonym(fdat, processedName);
                                if (nameSynonym != null) {
                                    getChadoDBConverter().store(nameSynonym);
                                }
                            }
                        }
                    }
                }
            }
        }

        addToFeatureMap(featureId, fdat);

        return true;
    }

    /**
     * Add feature data to FeatureMap, can be overidden by subclasses that need to store some
     * features in additional maps.
     * @param featureId the chado feature id
     * @param fdat feature information
     */
    protected void addToFeatureMap(Integer featureId, FeatureData fdat) {
        featureMap.put(featureId, fdat);
    }

    /**
     * Set the given attribute if the FeatureData says it's not set, then set the flag in
     * FeatureData to say it's set.
     */
    private void setAttributeIfNotSet(FeatureData fdat, final String attributeName,
				      final String value) throws ObjectStoreException {
        if (!fdat.getFlag(attributeName)) {
            setAttribute(fdat.getIntermineObjectId(), attributeName, value);
            fdat.setFlag(attributeName, true);
        }
    }

    /**
     * Create and store a new Item, returning a FeatureData object for the feature.
     *
     * @param featureId the chado id from the feature table
     * @param chadoType the type of the feature from the feature + cvterm tables
     * @param uniqueName the uniquename from chado
     * @param name the name from chado
     * @param md5checksum the checksum from the chado feature able
     * @param seqlen the length from the feature table
     * @param organismId the organism id of the feature from chado
     * @return a FeatureData object
     * @throws ObjectStoreException if there is a problem while storing
     */
    protected FeatureData makeFeatureData(int featureId, String chadoType, String uniqueName, String name, String md5checksum, int seqlen, int organismId) throws ObjectStoreException {
        String interMineType = TypeUtil.javaiseClassName(fixFeatureType(chadoType));
        OrganismData organismData = getChadoDBConverter().getChadoIdToOrgDataMap().get(new Integer(organismId));

        Item feature = makeFeature(new Integer(featureId), chadoType, interMineType, name, uniqueName, seqlen, organismData.getTaxonId());
        if (feature == null) {
            return null;
        }
        int taxonId = organismData.getTaxonId();
        FeatureData fdat = new FeatureData();
        Item organismItem = getChadoDBConverter().getOrganismItem(taxonId);
        feature.setReference("organism", organismItem);
        if (feature.checkAttribute("md5checksum")) {
            feature.setAttribute("md5checksum", md5checksum);
        }
        BioStoreHook.setSOTerm(getChadoDBConverter(), feature, chadoType, getChadoDBConverter().getSequenceOntologyRefId());
        fdat.setFieldExistenceFlags(feature);

        fdat.setIntermineObjectId(store(feature, taxonId));
        fdat.setItemIdentifier(feature.getIdentifier());
        fdat.setUniqueName(uniqueName);
        fdat.setChadoFeatureName(name);
        fdat.setInterMineType(feature.getClassName());
        fdat.organismData = organismData;
        fdat.setMd5checksum(md5checksum);
        return fdat;
    }

    /**
     * Store the feature Item.
     * @param feature the Item
     * @param taxonId the taxon id of this feature
     * @return the database id of the new Item
     * @throws ObjectStoreException if an error occurs while storing
     */
    protected Integer store(Item feature, int taxonId) throws ObjectStoreException {
        return getChadoDBConverter().store(feature);
    }

    /**
     * Make a new feature
     * @param featureId the chado feature id
     * @param chadoFeatureType the chado feature type (a SO term)
     * @param interMineType the InterMine type of the feature
     * @param name the name
     * @param uniqueName the uniquename
     * @param seqlen the sequence length (if known)
     * @param taxonId the NCBI taxon id of the current feature
     * @return the new feature
     */
    protected Item makeFeature(Integer featureId, String chadoFeatureType, String interMineType,
			       String name, String uniqueName,
			       int seqlen, int taxonId) {
        return getChadoDBConverter().createItem(interMineType);
    }

    /**
     * Get a list of the chado/so types of the desired features we wish to load.
     * @return the list of features
     */
    protected List<String> getFeatures() {
        return DESIRED_FEATURES;
    }

    /**
     * Return a list of types where one logical feature is represented as multiple rows in the
     * feature table.  An example is UTR features in flybase - if a UTR spans multiple exons, each
     * part is a separate row in the feature table.  In InterMine we represent the UTR as a single
     * object with a start and end so we need special handling for these types.
     * @return a list of segmented feature type
     */
    protected List<String> getSegmentedFeatures() {
        return new ArrayList<String>();
    }

    /**
     * Fix types from the feature table, perhaps by changing non-SO type into their SO equivalent.
     * Types that don't need fixing will be returned unchanged.
     * @param type the input type
     * @return the fixed type
     */
    protected String fixFeatureType(String type) {
        if ("five_prime_untranslated_region".equals(type)) {
            return "five_prime_UTR";
        }
        if ("three_prime_untranslated_region".equals(type)) {
            return "three_prime_UTR";
        }
        if ("full_transcript".equals(type)) {
            return "mature_transcript";
        }
        return type;
    }

    /**
     * Do any extra processing that is needed before the converter starts querying features
     * @param connection the Connection
     * @throws ObjectStoreException if there is a object store problem
     * @throws SQLException if there is a database problem
     */
    protected void earlyExtraProcessing(Connection connection) throws ObjectStoreException, SQLException {
        // override in subclasses as necessary
    }

    /**
     * Do any extra processing for this database, after all other processing is done
     * @param connection the Connection
     * @param featureDataMap a map from chado feature_id to data for that feature
     * @throws ObjectStoreException if there is a problem while storing
     * @throws SQLException if there is a problem
     */
    protected void extraProcessing(Connection connection, Map<Integer, FeatureData> featureDataMap)
        throws ObjectStoreException, SQLException {
        // override in subclasses as necessary
    }

    /**
     * Perform any actions needed after all processing is finished.
     * @param connection the Connection
     * @param featureDataMap a map from chado feature_id to data for that feature
     * @throws SQLException if there is a problem
     */
    protected void finishedProcessing(Connection connection, Map<Integer, FeatureData> featureDataMap) throws SQLException {
        // connection will be null for tests
        if (connection != null) {
            String query = "DROP TABLE " + tempFeatureTableName;
            Statement stmt = connection.createStatement();
            LOG.info("executing: " + query);
            stmt.execute(query);
        }
    }

    /**
     * Process a featureloc table and create Location objects.
     * @param connection the Connection
     * @param res a ResultSet that has the columns: featureloc_id, feature_id, srcfeature_id, fmin, fmax, strand
     * @throws SQLException if there is a problem while querying
     * @throws ObjectStoreException if there is a problem while storing
     */
    protected void processLocationTable(Connection connection) throws SQLException, ObjectStoreException {
	LOG.info("processLocationTable:START");

        ResultSet res = getFeatureLocResultSet(connection);

	int count = 0;
        int featureWarnings = 0;
	// put source feature InterMine reference and location strings into maps
	Map<String,String> srcFeatureIMReferenceMap = new HashMap<String,String>();
	Map<String,String> srcFeatureIMLocationMap = new HashMap<String,String>();
	for (int i=0; i<SRC_FEATURE_IM_CLASSES.length; i++) {
	    srcFeatureIMReferenceMap.put(SRC_FEATURE_IM_CLASSES[i], SRC_FEATURE_IM_REFERENCES[i]);
	    srcFeatureIMLocationMap.put(SRC_FEATURE_IM_CLASSES[i], SRC_FEATURE_IM_LOCATIONS[i]);
	}
        while (res.next()) {
            Integer featureLocId = new Integer(res.getInt("featureloc_id"));
            Integer featureId = new Integer(res.getInt("feature_id"));
            Integer srcFeatureId = new Integer(res.getInt("srcfeature_id"));
            int start = res.getInt("fmin") + 1;
            int end = res.getInt("fmax");
            int strand = res.getInt("strand");
            if (featureMap.containsKey(srcFeatureId)) {
                FeatureData srcFeatureData = featureMap.get(srcFeatureId);
                if (featureMap.containsKey(featureId)) {
                    FeatureData featureData = featureMap.get(featureId);
                    int taxonId = featureData.organismData.getTaxonId();
                    Item location = makeLocation(start, end, strand, srcFeatureData, featureData, taxonId, featureId);
		    if (location!=null) {
			getChadoDBConverter().store(location); // location should never be null
		    } else {
			LOG.warn("processLocationTable: NULL location for featureLocId="+featureLocId);
		    }
                    final String featureClassName = getModel().getPackageName()+"."+featureData.getInterMineType();
                    Class<?> featureClass;
                    try {
                        featureClass = Class.forName(featureClassName);
                    } catch (ClassNotFoundException e) {
			LOG.error("processLocationTable:can't find "+featureClassName);
                        throw new RuntimeException("processLocationTable:unable to find class object for setting the feature reference", e);
                    }
		    Integer featureIntermineObjectId = featureData.getIntermineObjectId();
		    // locate this feature on the srcFeature
		    if (srcFeatureIMReferenceMap.containsKey(srcFeatureData.getInterMineType())) {
			Reference srcReference = new Reference();
			srcReference.setName(srcFeatureIMReferenceMap.get(srcFeatureData.getInterMineType()));
			srcReference.setRefId(srcFeatureData.getItemIdentifier());
			getChadoDBConverter().store(srcReference, featureIntermineObjectId);
			if (location!=null) {
			    Reference locReference = new Reference();
			    locReference.setName(srcFeatureIMLocationMap.get(srcFeatureData.getInterMineType()));
			    locReference.setRefId(location.getIdentifier());
			    getChadoDBConverter().store(locReference, featureIntermineObjectId);
			}
		    }
		    // determine length from start and end if not given in chado record
		    if (!featureData.getFlag(FeatureData.LENGTH_SET)) {
			setAttribute(featureData.getIntermineObjectId(), "length", String.valueOf(end - start + 1));
		    }
		    // count total features placed on locations
                    count++;
                } else {
                    if (featureWarnings<=20) {
                        if (featureWarnings<20) {
                            LOG.warn("processLocationTable:featureId ("+featureId+") from location "+featureLocId+" was not found in the feature table");
                        } else {
                            LOG.warn("processLocationTable:further location warnings ignored.");
                        }
                        featureWarnings++;
                    }
                }
            } else {
		LOG.error("processLocationTable:srcfeature_id ("+srcFeatureId+") from location "+featureLocId+" was not found in the feature table.");
                // don't error out because of a missing source feature
		// throw new RuntimeException("processLocationTable:srcfeature_id ("+srcFeatureId+") from location "+featureLocId+" was not found in the feature table.");
            }
        } // while res.next()
	res.close();
        LOG.info("processLocationTable:DONE, created "+count+" locations.");
    }

    /**
     * Make a Location of a feature on a SequenceFeature.
     * @param start the start position
     * @param end the end position
     * @param strand the strand
     * @param srcFeatureData the FeatureData for the src feature
     * @param featureData the FeatureData for the SequenceFeature
     * @param taxonId of the source feature
     * @param featureId id of feature
     * @return the new Location object
     * @throws ObjectStoreException if there is a problem while storing
     */
    protected Item makeLocation(int start, int end, int strand, FeatureData srcFeatureData, FeatureData featureData, int taxonId, int featureId) throws ObjectStoreException {
        Item location = getChadoDBConverter().makeLocation(srcFeatureData.getItemIdentifier(), featureData.getItemIdentifier(), start, end, strand, taxonId);
        return location;
    }

    /**
     * Use the feature_relationship table to set relations (references and collections) between features.
     */
    private void processRelationTable(Connection connection, boolean subjectFirst) throws SQLException, ObjectStoreException {
	LOG.info("processRelationTable:START; subjectFirst="+subjectFirst);
		
        ResultSet res = getFeatureRelationshipResultSet(connection, subjectFirst);
        Integer lastSubjectId = null;

        // Map from relation type to Map from object type to FeatureData - used to collect up all
        // the collection/reference information for ONE subject feature
        Map<String, Map<String, List<FeatureData>>> relTypeMap = new HashMap<String, Map<String, List<FeatureData>>>();
		
        int featureWarnings = 0;
        int collectionWarnings = 0;
        int count = 0;
        int collectionTotal = 0;

        while (res.next()) {
			
	    // load the data
            Integer featRelationshipId = new Integer(res.getInt("feature_relationship_id"));
            Integer firstFeatureId = new Integer(res.getInt("feature1_id"));
            Integer secondFeatureId = new Integer(res.getInt("feature2_id"));
            String relationTypeName = res.getString("type_name");

            if (lastSubjectId != null && !firstFeatureId.equals(lastSubjectId)) {
                if (!processCollectionData(lastSubjectId, relTypeMap, collectionWarnings, subjectFirst)) collectionWarnings++;
                collectionTotal += relTypeMap.size();
                relTypeMap = new HashMap<String, Map<String, List<FeatureData>>>();
            }

	    // skip out if this is an undesired relation
	    if (!DESIRED_RELATIONS.contains(relationTypeName)) {
		if (DEBUG) debugToLog("processRelationTable: skipped undesired relation "+relationTypeName);
		continue;
	    }
			
	    // special case for part_of type relations - they are directional so only keep subjectFirst
            if (PARTOF_RELATIONS.contains(relationTypeName) && !subjectFirst) {
		if (DEBUG) debugToLog("processRelationTable: skipped relation "+relationTypeName+" because subjectFirst="+subjectFirst);
                continue;
            }

	    if (DEBUG) debugToLog("processRelationTable: processing relation "+relationTypeName+"; subjectFirst="+subjectFirst);

            if (featureMap.containsKey(firstFeatureId) && featureMap.containsKey(secondFeatureId)) {
		// both first and second feature found in featureMap
		FeatureData objectFeatureData = featureMap.get(secondFeatureId);
		Map<String, List<FeatureData>> objectClassFeatureDataMap;
		if (relTypeMap.containsKey(relationTypeName)) {
		    objectClassFeatureDataMap = relTypeMap.get(relationTypeName);
		} else {	
		    // THIS ALWAYS HAPPENS
		    objectClassFeatureDataMap = new HashMap<String, List<FeatureData>>();
		    relTypeMap.put(relationTypeName, objectClassFeatureDataMap);
		}
		List<FeatureData> featureDataList;
		String objectFeatureType = objectFeatureData.getInterMineType();
		if (objectClassFeatureDataMap.containsKey(objectFeatureType)) {
		    featureDataList = objectClassFeatureDataMap.get(objectFeatureType);
		} else {
		    // THIS ALWAYS HAPPENS
		    featureDataList = new ArrayList<FeatureData>();
		    objectClassFeatureDataMap.put(objectFeatureType, featureDataList);
		}
		featureDataList.add(objectFeatureData);
	    } else if (!featureMap.containsKey(secondFeatureId)) {
		// second feature missing from featureMap!
		if (featureWarnings<20) LOG.warn("object_id " + secondFeatureId + " from feature_relationship " + featRelationshipId + " was not found in the feature table");
		if (featureWarnings==20) LOG.warn("further feature_relationship warnings ignored");
		featureWarnings++;
	    } else {
		// first feature missing from featureMap!
                if (featureWarnings<20) LOG.warn("subject_id " + firstFeatureId + " from feature_relationship " + featRelationshipId + " was not found in the feature table");
		if (featureWarnings==20) LOG.warn("further feature_relationship warnings ignored");
		featureWarnings++;
            }
            count++;
            lastSubjectId = firstFeatureId;
        }
	// get last one
        if (lastSubjectId != null) {
            processCollectionData(lastSubjectId, relTypeMap, collectionWarnings, subjectFirst);
            collectionTotal += relTypeMap.size();
        }

	res.close();
        LOG.info("processRelationTable:END; subjectFirst="+subjectFirst+"; processed "+count+" relations; total collection elements created: "+collectionTotal);
    }

    /**
     * Create collections and references for the Item given by chadoSubjectId.
     * @param collectionWarnings
     */
    private boolean processCollectionData(Integer chadoSubjectId, Map<String, Map<String, List<FeatureData>>> relTypeMap, int collectionWarnings, boolean subjectIsFirst) throws ObjectStoreException {
        FeatureData subjectData = featureMap.get(chadoSubjectId);
        if (subjectData == null) {
            if (collectionWarnings<20) LOG.warn("unknown feature " + chadoSubjectId + " passed to processCollectionData - ignoring");
	    if (collectionWarnings==20) LOG.warn("ignoring further unknown feature warnings from processCollectionData()");
	    return false;
        }

        // map from collection name to list of item ids
        Map<String, List<String>> collectionsToStore = new HashMap<String, List<String>>();

        String subjectInterMineType = subjectData.getInterMineType();
        ClassDescriptor cd = getModel().getClassDescriptorByName(subjectInterMineType);
        Integer intermineObjectId = subjectData.getIntermineObjectId();
        for (Map.Entry<String, Map<String, List<FeatureData>>> entry: relTypeMap.entrySet()) {
            String relationType = entry.getKey();
            if (DEBUG) debugToLog("processCollectionData:relationType="+relationType);
            Map<String, List<FeatureData>> objectClassFeatureDataMap = entry.getValue();

            Set<Entry<String, List<FeatureData>>> mapEntries = objectClassFeatureDataMap.entrySet();
            for (Map.Entry<String, List<FeatureData>> featureDataMap: mapEntries) {
                String objectClass = featureDataMap.getKey();
                if (DEBUG) debugToLog("processCollectionData:objectClass="+objectClass);
                List<FeatureData> featureDataCollection = featureDataMap.getValue();
                List<FieldDescriptor> fds = null;

                FeatureData subjectFeatureData = featureMap.get(chadoSubjectId);
                // key example: ("relationship", "Polypeptide", "producedby", "MRNA")
                String relType;
                if (subjectIsFirst) {
                    relType = "relationship";
                } else {
                    relType = "rev_relationship";
                }
                MultiKey key = new MultiKey(relType, subjectFeatureData.getInterMineType(), relationType, objectClass);
                List<ConfigAction> actionList = getConfig(subjectData.organismData.getTaxonId()).get(key);

                if (actionList != null) {
                    if (actionList.size()==0 || actionList.size()==1 && actionList.get(0) instanceof DoNothingAction) {
                        // do nothing
			if (DEBUG) debugToLog("processCollectionData:"+key.toString()+" actionList!=null, do nothing");
                        continue;
                    }
		    if (DEBUG) debugToLog("processCollectionData:"+key.toString()+" actionList!=null, building fds");
                    fds = new ArrayList<FieldDescriptor>();
                    for (ConfigAction action: actionList) {
                        if (action instanceof SetFieldConfigAction) {
                            SetFieldConfigAction setAction = (SetFieldConfigAction) action;
                            String fieldName = setAction.getFieldName();
                            FieldDescriptor fd = cd.getFieldDescriptorByName(fieldName);
                            if (fd == null) {
                                throw new RuntimeException("can't find field " + fieldName + " in class " + cd + " configured for " + key);
                            }
                            fds.add(fd);
                        }
                    }
                    if (fds.size() == 0) {
                        throw new RuntimeException("no actions found for " + key);
                    }
                } else {
		    if (DEBUG) debugToLog("processCollectionData:"+key.toString()+" actionList==null");
                    if (PARTOF_RELATIONS.contains(relationType)) {
                        // special case for part_of relations - try to find a reference or
                        // collection that has a name that looks right for these objects (of class
                        // objectClass).  eg.  If the subject is a Transcript and the objectClass
                        // is Exon then find collections called "exons", "geneParts" (GenePart is
                        // a superclass of Exon)
                        fds = getReferenceForRelationship(objectClass, cd);
                    } else {
                        continue;
                    }
                }

                if (fds.size() == 0) {
                    if (!loggedMissingCols.contains(subjectInterMineType + relationType)) {
                        LOG.error("can't find collection for type "+objectClass+" with relationship "+relationType+" in "+subjectInterMineType+" (was processing feature "+chadoSubjectId+")");
                        loggedMissingCols.add(subjectInterMineType + relationType);
                    }
                    continue;
                }

                for (FieldDescriptor fd: fds) {

                    if (DEBUG) debugToLog("processCollectionData:fd="+fd.getName());
                    
                    if (fd.isReference()) {
                        if (objectClassFeatureDataMap.size() > 1) {
                            throw new RuntimeException("found more than one object for reference "+fd+" in class "+subjectInterMineType+" current subject identifier: "+subjectData.getUniqueName());
                        }
                        if (objectClassFeatureDataMap.size() == 1) {
                            // create a Reference
                            Reference reference = new Reference();
                            reference.setName(fd.getName());
                            FeatureData referencedFeatureData = featureDataCollection.get(0);
                            reference.setRefId(referencedFeatureData.getItemIdentifier());
                            // store the Reference
                            getChadoDBConverter().store(reference, intermineObjectId);
                            // special case for 1-1 relations - we need to set the reverse reference
                            ReferenceDescriptor rd = (ReferenceDescriptor) fd;
                            ReferenceDescriptor reverseRD = rd.getReverseReferenceDescriptor();
                            if (reverseRD != null && !reverseRD.isCollection()) {
                                // create the reverse Reference
                                Reference revReference = new Reference();
                                revReference.setName(reverseRD.getName());
                                revReference.setRefId(subjectData.getItemIdentifier());
                                Integer refObjectId = referencedFeatureData.getIntermineObjectId();
				// store the reverse Reference
                                getChadoDBConverter().store(revReference, refObjectId);
                            }
                        }
                        
                    } else {
                        // store a Collection
                        List<String> itemIds;
                        if (collectionsToStore.containsKey(fd.getName())) {
                            itemIds = collectionsToStore.get(fd.getName());
                        } else {
                            itemIds = new ArrayList<String>();
                            collectionsToStore.put(fd.getName(), itemIds);
                        }
                        for (FeatureData featureData: featureDataCollection) {
                            itemIds.add(featureData.getItemIdentifier());
                        }
                    }
                }
            }
        }

        for (Map.Entry<String, List<String>> entry: collectionsToStore.entrySet()) {
            ReferenceList referenceList = new ReferenceList();
            String collectionName = entry.getKey();
            referenceList.setName(collectionName);
            List<String> idList = entry.getValue();
            referenceList.setRefIds(idList);
            getChadoDBConverter().store(referenceList, intermineObjectId);

            // if there is a field called <classname>Count that matches the name of the collection we just stored, set it
            String countName;
            if (collectionName.endsWith("s")) {
                countName = collectionName.substring(0, collectionName.length() - 1);
            } else {
                countName = collectionName;
            }
            countName += "Count";
            if (cd.getAttributeDescriptorByName(countName, true) != null) {
                setAttribute(intermineObjectId, countName, String.valueOf(idList.size()));
            }
        }

        return true;
    }

    /**
     * Search ClassDescriptor cd class for refs/collections with the right name for the objectType
     * eg. find CDSs collection for objectType = CDS and find gene reference for objectType = Gene.
     */
    private List<FieldDescriptor> getReferenceForRelationship(String objectType, ClassDescriptor cd) {
        List<FieldDescriptor> fds = new ArrayList<FieldDescriptor>();
        LinkedHashSet<String> allClasses = new LinkedHashSet<String>();
        allClasses.add(objectType);
        try {
            Set<String> parentClasses = ClassDescriptor.findSuperClassNames(getModel(), objectType);
            allClasses.addAll(parentClasses);
        } catch (MetaDataException e) {
            throw new RuntimeException("class not found in the model", e);
        }

        for (String clsName: allClasses) {
            List<String> possibleRefNames = new ArrayList<String>();
            String unqualifiedClsName = TypeUtil.unqualifiedName(clsName);
            possibleRefNames.add(unqualifiedClsName);
            possibleRefNames.add(unqualifiedClsName + 's');
            possibleRefNames.add(StringUtil.decapitalise(unqualifiedClsName));
            possibleRefNames.add(StringUtil.decapitalise(unqualifiedClsName) + 's');
            for (String possibleRefName: possibleRefNames) {
		if (DEBUG) debugToLog("getReferenceForRelationship("+objectType+","+cd.getName()+"):possibleRefName="+possibleRefName);
                FieldDescriptor fd = cd.getFieldDescriptorByName(possibleRefName);
                if (fd != null) {
                    fds.add(fd);
                }
            }
        }
        return fds;
    }

    @SuppressWarnings("boxing")
    private void processDbxrefTable(Connection connection)
        throws SQLException, ObjectStoreException {

        ResultSet res = getDbxrefResultSet(connection);
        Set<String> existingAttributes = new HashSet<String>();
        Integer currentFeatureId = null;
        int count = 0;

        while (res.next()) {
            Integer featureId = new Integer(res.getInt("feature_id"));
            String accession = res.getString("accession");
            String dbName = res.getString("db_name");
            Boolean isCurrent = res.getBoolean("is_current");

            if (currentFeatureId != null && currentFeatureId != featureId) {
                existingAttributes = new HashSet<String>();
            }

            if (featureMap.containsKey(featureId)) {
                FeatureData fdat = featureMap.get(featureId);
                MultiKey key = new MultiKey("dbxref", fdat.getInterMineType(), dbName, isCurrent);

                if (accession == null) {
                    throw new RuntimeException("found null accession in dbxref table for database "
					       + dbName + ".");
                }
                accession  = fixIdentifier(fdat, accession);

                int taxonId = fdat.organismData.getTaxonId();
                Map<MultiKey, List<ConfigAction>> orgConfig = getConfig(taxonId);
                List<ConfigAction> actionList = orgConfig.get(key);

                if (actionList == null) {
                    // try ignoring isCurrent
                    MultiKey key2 = new MultiKey("dbxref", fdat.getInterMineType(), dbName, null);
                    actionList = orgConfig.get(key2);
                }

                if (actionList == null) {
                    // no actions configured for this synonym
                    continue;
                }
                Set<String> fieldsSet = new HashSet<String>();
                for (ConfigAction action: actionList) {
                    if (action instanceof SetFieldConfigAction) {
                        SetFieldConfigAction setAction = (SetFieldConfigAction) action;
                        if (!existingAttributes.contains(setAction.getFieldName())) {
                            if (setAction.isValidValue(accession)) {
                                String newFieldValue = setAction.processValue(accession);
                                setAttribute(fdat.getIntermineObjectId(), setAction.getFieldName(), newFieldValue);
                                existingAttributes.add(setAction.getFieldName());
                                fieldsSet.add(newFieldValue);
                                if ("primaryIdentifier".equals(setAction.getFieldName())) {
                                    fdat.setFlag(FeatureData.IDENTIFIER_SET, true);
                                }
                            }
                        }
                    }
                }

                for (ConfigAction action: actionList) {
                    if (action instanceof CreateSynonymAction) {
                        CreateSynonymAction createSynonymAction = (CreateSynonymAction) action;
                        if (!createSynonymAction.isValidValue(accession)) {
                            continue;
                        }
                        String newFieldValue = createSynonymAction.processValue(accession);
                        if (fdat.getExistingSynonyms().contains(newFieldValue)) {
                            continue;
                        }
                        Item synonym = createSynonym(fdat, newFieldValue);
                        if (synonym != null) {
                            getChadoDBConverter().store(synonym);
                            count++;
                        }
                    }
                }
            }

            currentFeatureId = featureId;
        }

        LOG.info("created " + count + " synonyms from the dbxref table");
        res.close();
    }

    private void processFeaturePropTable(Connection connection)
        throws SQLException, ObjectStoreException {
        ResultSet res = getFeaturePropResultSet(connection);
        int count = 0;
        while (res.next()) {
            Integer featureId = new Integer(res.getInt("feature_id"));
            String identifier = res.getString("value");

            if (identifier == null) {
                continue;
            }

            String propTypeName = res.getString("type_name");

            if (featureMap.containsKey(featureId)) {

                FeatureData fdat = featureMap.get(featureId);
                MultiKey key = new MultiKey("prop", fdat.getInterMineType(), propTypeName);
                int taxonId = fdat.organismData.getTaxonId();
                List<ConfigAction> actionList = getConfig(taxonId).get(key);
                if (actionList == null) {
                    // no actions configured for this prop
                    continue;
                }

                Set<String> fieldsSet = new HashSet<String>();

                for (ConfigAction action: actionList) {
                    if (action instanceof SetFieldConfigAction) {
                        SetFieldConfigAction setAction = (SetFieldConfigAction) action;
                        if (setAction.isValidValue(identifier)) {
                            String newFieldValue = setAction.processValue(identifier);
                            setAttribute(fdat.getIntermineObjectId(), setAction.getFieldName(),
					 newFieldValue);
                            fieldsSet.add(newFieldValue);

                            if ("primaryIdentifier".equals(setAction.getFieldName())) {
                                fdat.setFlag(FeatureData.IDENTIFIER_SET, true);
                            }
                        }
                    }
                }

                for (ConfigAction action: actionList) {
                    if (action instanceof CreateSynonymAction) {
                        CreateSynonymAction synonymAction = (CreateSynonymAction) action;
                        if (!synonymAction.isValidValue(identifier)) {
                            continue;
                        }
                        String newFieldValue = synonymAction.processValue(identifier);
                        Set<String> existingSynonyms = fdat.getExistingSynonyms();
                        if (existingSynonyms.contains(newFieldValue)) {
                            continue;
                        }
                        Item synonym = createSynonym(fdat, newFieldValue);
                        if (synonym != null) {
                            getChadoDBConverter().store(synonym);
                            count++;
                        }
                    }
                }
            }
        }
        LOG.info("created " + count + " synonyms from the featureprop table");
        res.close();
    }

    /**
     * @param identifier identifier for term, eg. FB00004958
     * @return an id representing the term object
     * @throws ObjectStoreException if somethign goes wrong
     */
    protected String makeAnatomyTerm(String identifier)
        throws ObjectStoreException {
        // override in subclasses as necessary
        return null;
    }


    /**
     * Read the feature, feature_cvterm and cvterm tables, then set fields, create synonyms or
     * create objects based on the cvterms.
     * @param connection the Connection
     */
    private void processFeatureCVTermTable(Connection connection)
        throws SQLException, ObjectStoreException {
        ResultSet res = getFeatureCVTermResultSet(connection);
        int count = 0;
        Integer previousFeatureId = null;

        // map from reference/collection name to list of Items to store in the reference or collection
        Map<String, List<Item>> dataMap = new HashMap<String, List<Item>>();

        while (res.next()) {
            Integer featureId = new Integer(res.getInt("feature_id"));
            String cvtermName = res.getString("cvterm_name");
            String cvName = res.getString("cv_name");

            FeatureData fdat = featureMap.get(featureId);

            if (fdat == null) {
                continue;
            }

            if (!featureId.equals(previousFeatureId) && previousFeatureId != null) {
                processCVTermRefCols(previousFeatureId, dataMap);
                dataMap = new HashMap<String, List<Item>>();
            }

            MultiKey key = new MultiKey("cvterm", fdat.getInterMineType(), cvName);

            int taxonId = fdat.organismData.getTaxonId();

            List<ConfigAction> actionList = getConfig(taxonId).get(key);
            if (actionList == null) {
                // no actions configured for this prop
                continue;
            }

            Set<String> fieldsSet = new HashSet<String>();

            for (ConfigAction action: actionList) {
                if (action instanceof SetFieldConfigAction) {
                    SetFieldConfigAction setAction = (SetFieldConfigAction) action;
                    if (setAction.isValidValue(cvtermName)) {
                        String newFieldValue = setAction.processValue(cvtermName);
                        setAttribute(fdat.getIntermineObjectId(), setAction.getFieldName(),
				     newFieldValue);

                        fieldsSet.add(newFieldValue);
                        if ("primaryIdentifier".equals(setAction.getFieldName())) {
                            fdat.setFlag(FeatureData.IDENTIFIER_SET, true);
                        }
                    }
                } else {
                    if (action instanceof CreateSynonymAction) {
                        CreateSynonymAction synonymAction = (CreateSynonymAction) action;
                        if (!synonymAction.isValidValue(cvtermName)) {
                            continue;
                        }
                        String newFieldValue = synonymAction.processValue(cvtermName);
                        Set<String> existingSynonyms = fdat.getExistingSynonyms();
                        if (existingSynonyms.contains(newFieldValue)) {
                            continue;
                        }
                        Item synonym = createSynonym(fdat, newFieldValue);
                        if (synonym != null) {
                            getChadoDBConverter().store(synonym);
                            count++;
                        }
                    } else {
                        // TODO fixme
			//                        if (action instanceof CreateCollectionAction) {
			//                            CreateCollectionAction cca = (CreateCollectionAction) action;
			//
			//                            Item item = null;
			//                            String fieldName = cca.getFieldName();
			//                            String className = cca.getClassName();
			//                            if (cca.createSingletons()) {
			//                                MultiKey singletonKey =
			//                                        new MultiKey(className, fieldName, cvtermName);
			//                                item = (Item) singletonMap.get(singletonKey);
			//                            }
			//                            if (item == null) {
			//                                item = getChadoDBConverter().createItem(className);
			//                                item.setAttribute(fieldName, cvtermName);
			//                                getChadoDBConverter().store(item);
			//                                if (cca.createSingletons()) {
			//                                    singletonMap.put(key, item);
			//                                }
			//                            }
			//
			//                            String referenceName = cca.getReferenceName();
			//                            List<Item> itemList;
			//                            // creating collection, already seen this ref
			//                            if (dataMap.containsKey(referenceName)) {
			//                                itemList = dataMap.get(referenceName);
			//                                // new collection
			//                            } else {
			//                                itemList = new ArrayList<Item>();
			//                                dataMap.put(referenceName, itemList);
			//                            }
			//                            itemList.add(item);
			//                        }
                    }
                }
            }

            previousFeatureId = featureId;
        }

        if (previousFeatureId != null) {
            processCVTermRefCols(previousFeatureId, dataMap);
        }

        LOG.info("created " + count + " synonyms from the feature_cvterm table");
        res.close();
    }


    /**
     * Given the object id and a map of reference/collection names to Items, store the Items in the
     * reference or collection of the object.
     */
    private void processCVTermRefCols(Integer chadoObjectId, Map<String, List<Item>> dataMap)
        throws ObjectStoreException {

        FeatureData fdat = featureMap.get(chadoObjectId);
        String interMineType = fdat.getInterMineType();
        ClassDescriptor cd = getModel().getClassDescriptorByName(interMineType);
        for (String referenceName: dataMap.keySet()) {
            FieldDescriptor fd = cd.getFieldDescriptorByName(referenceName);
            if (fd == null) {
                throw new RuntimeException("failed to find " + referenceName + " in "
					   + interMineType);
            }
            List<Item> itemList = dataMap.get(referenceName);
            Integer intermineObjectId = fdat.getIntermineObjectId();
            if (fd.isReference()) {
                if (itemList.size() > 1) {
                    throw new RuntimeException("found more than one object for reference "
					       + fd + " in class "
					       + interMineType + " items: " + itemList);
                }
                Item item = itemList.iterator().next();
                Reference reference = new Reference();
                reference.setName(fd.getName());
                String itemIdentifier = item.getIdentifier();
                reference.setRefId(itemIdentifier);
                getChadoDBConverter().store(reference, intermineObjectId);

                // XXX FIXME TODO: special case for 1-1 relations - we need to set the reverse
                // reference

            } else {
                ReferenceList referenceList = new ReferenceList();
                referenceList.setName(referenceName);
                for (Item item: itemList) {
                    referenceList.addRefId(item.getIdentifier());
                }
                getChadoDBConverter().store(referenceList, intermineObjectId);
            }
        }
    }

    @SuppressWarnings("boxing")
    private void processSynonymTable(Connection connection)
        throws SQLException, ObjectStoreException {
        ResultSet res = getSynonymResultSet(connection);
        Set<String> existingAttributes = new HashSet<String>();
        Integer currentFeatureId = null;
        int count = 0;
        while (res.next()) {
            Integer featureId = new Integer(res.getInt("feature_id"));
            String identifier = res.getString("synonym_name");
            String synonymTypeName = res.getString("type_name");

            // change type so synonyms will merge correctly
            if ("primaryAccession".equals(synonymTypeName)) {
                synonymTypeName = "accession";
            } else if ("primaryIdentifier".equals(synonymTypeName)) {
                synonymTypeName = "identifier";
            }

            Boolean isCurrent = res.getBoolean("is_current");

            // it is a not null in db
            if (identifier == null) {
                throw new RuntimeException("found null synonym name in synonym table.");
            }
            identifier = fixIdentifier(featureMap.get(featureId), identifier);

            if (currentFeatureId != null && currentFeatureId != featureId) {
                existingAttributes = new HashSet<String>();
            }

            if (featureMap.containsKey(featureId)) {
                FeatureData fdat = featureMap.get(featureId);
                identifier = fixIdentifier(fdat, identifier);
                MultiKey key =
		    new MultiKey("synonym", fdat.getInterMineType(),
				 synonymTypeName, isCurrent);
                int taxonId = fdat.organismData.getTaxonId();
                Map<MultiKey, List<ConfigAction>> orgConfig = getConfig(taxonId);
                List<ConfigAction> actionList = orgConfig.get(key);

                if (actionList == null) {
                    // try ignoring isCurrent
                    MultiKey key2 =
			new MultiKey("synonym", fdat.getInterMineType(), synonymTypeName, null);
                    actionList = orgConfig.get(key2);
                }
                if (actionList == null) {
                    // no actions configured for this synonym
                    continue;
                }
                for (ConfigAction action: actionList) {
                    if (action instanceof SetFieldConfigAction) {
                        SetFieldConfigAction setAction = (SetFieldConfigAction) action;
                        if (!existingAttributes.contains(setAction.getFieldName())
			    && setAction.isValidValue(identifier)) {
                            String newFieldValue = setAction.processValue(identifier);
                            setAttribute(fdat.getIntermineObjectId(), setAction.getFieldName(),
					 newFieldValue);
                            existingAttributes.add(setAction.getFieldName());
                            if ("primaryIdentifier".equals(setAction.getFieldName())) {
                                fdat.setFlag(FeatureData.IDENTIFIER_SET, true);
                            }
                        }
                    }
                }

                for (ConfigAction action: actionList) {
                    if (action instanceof CreateSynonymAction) {
                        CreateSynonymAction createSynonymAction = (CreateSynonymAction) action;
                        if (!createSynonymAction.isValidValue(identifier)) {
                            continue;
                        }
                        String newFieldValue = createSynonymAction.processValue(identifier);
                        if (fdat.getExistingSynonyms().contains(newFieldValue)) {
                            continue;
                        }
                        Item synonym =
			    createSynonym(fdat, newFieldValue);
                        if (synonym != null) {
                            getChadoDBConverter().store(synonym);
                            count++;
                        }
                    }
                }
            }
            currentFeatureId = featureId;
        }

        LOG.info("created " + count + " synonyms from the synonym table");
        res.close();
    }

    /**
     * Process the identifier and return a "cleaned" version.  Implement in sub-classes to fix
     * data problem.
     * @param featureData the FeatureData object of the feature that this identifier came from
     * @param identifier the identifier
     * @return a cleaned identifier
     */
    protected String fixIdentifier(FeatureData featureData, String identifier) {
        return identifier;
    }

    /**
     * Process the identifier and return a "cleaned" version.  Implement in sub-classes to fix
     * data problem.
     * @param featureData the FeatureData object of the feature that this identifier came from
     * @param identifier the identifier
     * @param prefix Needed for modencode to distinguish gene names coming from a gene model
     * @return a cleaned identifier
     */
    protected String fixIdentifier(FeatureData featureData, String identifier, String prefix) {
        return identifier;
    }

    private void processPubTable(Connection connection)
        throws SQLException, ObjectStoreException {
        ResultSet res = getPubResultSet(connection);

        List<String> currentPublicationIds = new ArrayList<String>();
        Integer lastPubFeatureId = null;
        int featureWarnings = 0;
        int count = 0;

        while (res.next()) {
            Integer featureId = new Integer(res.getInt("feature_id"));
            if (!featureMap.containsKey(featureId)) {
                if (featureWarnings <= 20) {
                    if (featureWarnings < 20) {
                        LOG.warn("feature " + featureId + " not found in features Map while "
				 + "processing publications");
                    } else {
                        LOG.warn("further feature id warnings ignored in processPubTable()");
                    }
                    featureWarnings++;
                }
                continue;
            }
            Integer pubMedId = fixPubMedId(res.getString("pub_db_identifier"));
            if (lastPubFeatureId != null && !featureId.equals(lastPubFeatureId)) {
                makeFeaturePublications(lastPubFeatureId, currentPublicationIds);
                currentPublicationIds = new ArrayList<String>();
            }
            String publicationId = makePublication(pubMedId);
            currentPublicationIds.add(publicationId);
            lastPubFeatureId = featureId;
            count++;
        }

        if (lastPubFeatureId != null) {
            makeFeaturePublications(lastPubFeatureId, currentPublicationIds);
        }
        LOG.info("Created " + count + " publications");
        res.close();
    }

    /**
     * Parser a pubmed id from a results string, extracted to a method so subclasses can override
     * and fix prefixed pubmed ids.
     * @param pubmedStr id fetched from database
     * @return the pubmed id
     */
    @SuppressWarnings("boxing")
    protected Integer fixPubMedId(String pubmedStr) {
        return Integer.parseInt(pubmedStr);
    }

    /**
     * Return the item identifier of the publication Item for the given pubmed id.
     * @param pubMedId the pubmed id
     * @return the publication item id
     * @throws ObjectStoreException if the item can't be stored
     */
    protected String makePublication(Integer pubMedId) throws ObjectStoreException {
        if (publications.containsKey(pubMedId)) {
            return publications.get(pubMedId);
        }
        Item publication = getChadoDBConverter().createItem("Publication");
        publication.setAttribute("pubMedId", pubMedId.toString());
        getChadoDBConverter().store(publication); // Stores Publication
        String publicationId = publication.getIdentifier();
        publications.put(pubMedId, publicationId);
        return publicationId;
    }

    /**
     * Set the publications collection of the feature with the given (chado) feature id.
     */
    private void makeFeaturePublications(Integer featureId, List<String> argPublicationIds)
        throws ObjectStoreException {
        FeatureData fdat = featureMap.get(featureId);
        if (fdat == null) {
            throw new RuntimeException("feature " + featureId + " not found in features Map");
        }
        if (argPublicationIds.size() == 0) {
            return;
        }
        List<String> publicationIds = new ArrayList<String>(argPublicationIds);
        ReferenceList referenceList = new ReferenceList();
        referenceList.setName("publications");
        referenceList.setRefIds(publicationIds);
        getChadoDBConverter().store(referenceList, fdat.getIntermineObjectId());
    }

    /**
     * Return the interesting rows from the features table.
     * This is a protected method so that it can be overridden for testing
     * @param connection the db connection
     * @return the SQL result set
     * @throws SQLException if a database problem occurs
     */
    protected ResultSet getFeatureTableResultSet(Connection connection)
        throws SQLException {
        String query = "SELECT * FROM " + tempFeatureTableName;
        LOG.info("executing getFeatureTableResultSet(): " + query);
        Statement stmt = connection.createStatement();
        ResultSet res = stmt.executeQuery(query);
        return res;
    }

    /**
     * Convert the list of features to a string to be used in a SQL query.
     * @return the list of features as a string (in SQL list format)
     */
    private String getFeaturesString(List<String> featuresList) {
        StringBuffer featureListString = new StringBuffer();
        Iterator<String> i = featuresList.iterator();
        while (i.hasNext()) {
            String item = i.next();
            featureListString.append("'" + item + "'");
            if (i.hasNext()) {
                featureListString.append(", ");
            }
        }
        return featureListString.toString();
    }

    /**
     * Return a comma-separated string containing the organism_ids which to query from chado.
     */
    private String getOrganismIdsString() {
        return StringUtil.join(getChadoDBConverter().getChadoIdToOrgDataMap().keySet(), ", ");
    }

    /**
     * Create a temporary table containing only the features that interest us.  Also create indexes
     * for the type and feature_id columns.
     * The table is used in later queries.  This is a protected method so that it can be overridden
     * for testing.
     * @param connection the Connection
     * @throws SQLException if there is a problem
     */
    protected void createFeatureTempTable(Connection connection) throws SQLException {
        List<String> featuresList = new ArrayList<String>(getFeatures()); // DESIRED_FEATURES
        String featureTypesString = getFeaturesString(featuresList);      // comma-delimited list of features
        String organismConstraint = getOrganismConstraint();              // SQL restricting to selected organism(s)
        String orgConstraintForQuery = "";
        if (!StringUtils.isEmpty(organismConstraint)) {
            orgConstraintForQuery = " AND " + organismConstraint;
        }

	// create temp table holding all desired features
        String query =
	    "CREATE TEMPORARY TABLE " + tempFeatureTableName + " AS"
	    + " SELECT feature_id, feature.name, uniquename, cvterm.name as type,"
	    + " seqlen, is_analysis, residues, md5checksum, organism_id"
	    + " FROM feature, cvterm"
	    + " WHERE cvterm.name IN (" + featureTypesString  + ")"
	    + orgConstraintForQuery
	    + " AND NOT feature.is_obsolete"
	    + " AND feature.type_id = cvterm.cvterm_id "
	    + (getExtraFeatureConstraint() != null
	       ? " AND (" + getExtraFeatureConstraint() + ")"
	       : "");
        Statement stmt = connection.createStatement();
        LOG.info("executing createFeatureTempTable(): " + query);
        stmt.execute(query);
	// create indexes and analyze on temp table
        String idIndexQuery = "CREATE INDEX "+tempFeatureTableName+"_feature_index ON "+tempFeatureTableName+" (feature_id)";
        LOG.info("executing: " + idIndexQuery);
        stmt.execute(idIndexQuery);
        String typeIndexQuery = "CREATE INDEX "+tempFeatureTableName+"_type_index ON "+tempFeatureTableName+ " (type)";
        LOG.info("executing: " + typeIndexQuery);
        stmt.execute(typeIndexQuery);
        String analyze = "ANALYZE " + tempFeatureTableName;
        LOG.info("executing: " + analyze);
        stmt.execute(analyze);
    }

    /**
     * Return some SQL that can be included in the WHERE part of query that restricts features
     * by organism.  "organism_id" must be selected.
     * @return the SQL
     */
    protected String getOrganismConstraint() {
        String organismIdsString = getOrganismIdsString();
        if (StringUtils.isEmpty(organismIdsString)) {
            return "";
        }
        return "organism_id IN (" + organismIdsString + ")";
    }

    /**
     * Return an extra constraint to be used when querying the feature table.  Any feature table
     * column or cvterm table column can be constrained.  The cvterm will match the type_id field
     * in the feature.
     * eg. "uniquename not like 'BAD_ID%'"
     * @return the constraint as SQL or nul if there is no extra constraint.
     */
    protected String getExtraFeatureConstraint() {
        // no default
        return null;
    }

    /**
     * Return an SQL query that finds the feature_ids of all rows from the feature table that we
     * need to process.
     * @return the SQL string
     */
    protected String getFeatureIdQuery() {
        return "SELECT feature_id FROM " + tempFeatureTableName;
    }


    private String getSrcFeatureIdQuery() {
        return
	    "SELECT feature_id FROM feature, cvterm"
	    + "  WHERE type_id = cvterm.cvterm_id"
	    + "    AND feature.is_obsolete = 'f' "
	    + "    AND cvterm.name IN (" + getFeaturesString(SRC_FEATURE_TYPES) + ")"
	    + (getExtraFeatureConstraint() != null
	       ? " AND (" + getExtraFeatureConstraint() + ")"
	       : "");
    }

    /**
     * Return the interesting rows from the feature_relationship table.  The feature pairs are
     * returned in both subject, object and object, subject orientations so that the relationship
     * processing can be configured in a natural way.
     * This is a protected method so that it can be overridden for testing
     * @param connection the db connection
     * @param subjectFirst if true the subject_id column from the relationship table will be before
     *   the object_id in the results, otherwise it will be after.  ie.
     *   "feature_relationship_id, subject_id as feature1_id, object_id as feature2_id, ..."
     *   vs "feature_relationship_id, object_id as feature1_id, subject_id as feature2_id, ..."
     * @return the SQL result set
     * @throws SQLException if a database problem occurs
     */
    protected ResultSet getFeatureRelationshipResultSet(Connection connection, boolean subjectFirst) throws SQLException {
        String subObjString;
        if (subjectFirst) {
            subObjString = "subject_id as feature1_id, object_id as feature2_id";
        } else {
            subObjString = "object_id as feature1_id, subject_id as feature2_id";
        }
        String extraQueryBits = "";
	if (subjectFirst) {
	    extraQueryBits = getGenesPolypeptidesQuery(); // Why only subjectFirst?
	}
        String query = "SELECT feature_relationship_id, "+subObjString
	    + ", cvterm.name AS type_name"
	    + "  FROM feature_relationship, cvterm"
	    + "  WHERE cvterm.cvterm_id = type_id"
	    + "      AND subject_id IN (" + getFeatureIdQuery() + ")" // subject must be in feature temp table
	    + "      AND object_id IN (" + getFeatureIdQuery() + ")"  // object must be in feature temp table
	    + extraQueryBits                                          // get gene-polypeptide relations via mRNA pivot
	    + " ORDER BY feature1_id";
        LOG.info("getFeatureRelationshipResultSet() executing:\n" + query);
        Statement stmt = connection.createStatement();
        ResultSet res = stmt.executeQuery(query);
        return res;
    }

    /**
     * Return extra SQL that will be added to the feature_relationships table query.  It will add
     * fake relationships to the query to make it look like there are gene <-> polypeptide relations
     * with type 'producedby'.
     * Results in a resultset like (0, gene_id, polypeptide_id, 'producedby')
     * UPDATE: keep this specific to genes that are in the feature temp table to limit to chosen organisms; query originally hit entire database for all organisms.
     * UPDATE: simplified query using explicit type_id
     * UPDATE: UNION not UNION ALL (why would we want duplicate rows?)
     */
    private String getGenesPolypeptidesQuery() {
	// TEST: switched f1.feature_id and f3.feature_id in output
        return " UNION "
	    + " SELECT 0 AS feature_relationship_id, f3.feature_id AS feature1_id, f1.feature_id AS feature2_id, 'producedby' AS type_name "
            + " FROM "+tempFeatureTableName+" f1, feature_relationship fr1,"
            + "      "+tempFeatureTableName+" f2, feature_relationship fr2,"
            + "      "+tempFeatureTableName+" f3 "

            + " WHERE f1.feature_id = fr1.object_id"  // gene is object of rel 1
            + " AND   f2.feature_id = fr1.subject_id" // mRNA is subject of rel 1
            + " AND   f2.feature_id = fr2.object_id"  // mRNA is object of rel 2
            + " AND   f3.feature_id = fr2.subject_id" // polypeptide is subject of rel 2

            + " AND fr1.type_id = 42871"  // mRNA part_of gene, part_of=42871
            + " AND fr2.type_id = 42884"  // polypeptide derives_from mRNA, derives_from=42884

            + " AND f1.type = 'gene'"           // f1 is gene
            + " AND f2.type = 'mRNA'"           // f2 is mRNA
            + " AND f3.type = 'polypeptide'";   // f3 is polypeptide
    }

    /**
     * Return the interesting rows from the featureloc table.
     * This is a protected method so that it can be overridden for testing
     * @param connection the db connection
     * @return the SQL result set
     * @throws SQLException if a database problem occurs
     */
    protected ResultSet getFeatureLocResultSet(Connection connection) throws SQLException {
        String query =
	    "SELECT featureloc_id, feature_id, srcfeature_id, fmin, is_fmin_partial, fmax, is_fmax_partial, strand"
	    + " FROM featureloc"
	    + " WHERE locgroup=0"
	    + " AND feature_id IN ("+getFeatureIdQuery()+")"
	    + " AND srcfeature_id IN ("+getFeatureIdQuery()+")"
	    + " ORDER BY featureloc_id";
	LOG.info("getFeatureLocResultSet():query="+query);
        Statement stmt = connection.createStatement();
        ResultSet res = stmt.executeQuery(query);
	LOG.info("getFeatureLocResultSet():DONE.");
        return res;
    }


    /**
     * Return the interesting matches from the featureloc and feature tables.
     * feature<->featureloc<->match_feature<->featureloc<->feature
     * This is a protected method so that it can be overridden for testing
     * @param connection the db connection
     * @return the SQL result set
     * @throws SQLException if a database problem occurs
     */
    /* doesn't seem to be used
       protected ResultSet getMatchLocResultSet(Connection connection) throws SQLException {
       String query =
       "SELECT f1loc.featureloc_id, f1loc.srcfeature_id as feature_id,"
       + " f2loc.srcfeature_id AS srcfeature_id, f2loc.fmin,"
       + " false AS is_fmin_partial, f2loc.fmax, false AS is_fmax_partial,"
       + " f2loc.strand"
       + " FROM feature match, featureloc f1loc, featureloc f2loc, cvterm mt"
       + " WHERE match.feature_id = f1loc.feature_id"
       + " AND match.feature_id = f2loc.feature_id"
       + " AND match.type_id = mt.cvterm_id AND mt.name IN ('match', 'cDNA_match')"
       + " AND f1loc.srcfeature_id <> f2loc.srcfeature_id"
       + " AND f1loc.srcfeature_id IN (" + getFeatureIdQuery() + ")"
       + " AND f2loc.srcfeature_id IN (" + getSrcFeatureIdQuery() + ")";

       // Previous query included feature table three times
       // "SELECT f1loc.featureloc_id, f1.feature_id, f2.feature_id AS srcfeature_id, f2loc.fmin,"
       //  + "     false AS is_fmin_partial, f2loc.fmax, false AS is_fmax_partial, f2loc.strand"
       //  + "   FROM feature match, feature f1, featureloc f1loc, feature f2, featureloc f2loc,"
       //  + "        cvterm mt"
       //  + "  WHERE match.feature_id = f1loc.feature_id AND match.feature_id = f2loc.feature_id"
       //  + "    AND f1loc.srcfeature_id = f1.feature_id AND f2loc.srcfeature_id = f2.feature_id"
       //  + "    AND match.type_id = mt.cvterm_id AND mt.name IN ('match', 'cDNA_match')"
       //  + "    AND f1.feature_id <> f2.feature_id"
       //  + "    AND f1.feature_id IN (" + getFeatureIdQuery() + ")"
       //  + "    AND f2.feature_id IN (" + getSrcFeatureIdQuery() + ")";
       LOG.info("executing getMatchLocResultSet(): " + query);
       Statement stmt = connection.createStatement();
       ResultSet res = stmt.executeQuery(query);
       return res;
       }
    */

    /**
     * Return the interesting rows from the dbxref table.
     * This is a protected method so that it can be overridden for testing
     * @param connection the db connection
     * @return the SQL result set
     * @throws SQLException if a database problem occurs
     */
    protected ResultSet getDbxrefResultSet(Connection connection) throws SQLException {
        String query =
	    "SELECT feature.feature_id, accession, db.name AS db_name, is_current"
	    + "  FROM dbxref, feature_dbxref, feature, db"
	    + "  WHERE feature_dbxref.dbxref_id = dbxref.dbxref_id "
	    + "    AND feature_dbxref.feature_id = feature.feature_id "
	    + "    AND feature.feature_id IN"
	    + "        (" + getFeatureIdQuery() + ")"
	    + "    AND dbxref.db_id = db.db_id";
        LOG.info("executing getDbxrefResultSet(): " + query);
        Statement stmt = connection.createStatement();
        ResultSet res = stmt.executeQuery(query);
        return res;
    }

    /**
     * Return the interesting rows from the featureprop table.
     * This is a protected method so that it can be overridden for testing
     * @param connection the db connection
     * @return the SQL result set
     * @throws SQLException if a database problem occurs
     */
    protected ResultSet getFeaturePropResultSet(Connection connection) throws SQLException {
        String query =
	    "select feature_id, value, cvterm.name AS type_name FROM featureprop, cvterm"
	    + "   WHERE featureprop.type_id = cvterm.cvterm_id"
	    + "       AND feature_id IN (" + getFeatureIdQuery() + ")";
        LOG.info("executing getFeaturePropResultSet(): " + query);
        Statement stmt = connection.createStatement();
        ResultSet res = stmt.executeQuery(query);
        return res;
    }

    //    SELECT fp.value
    //    FROM feature f, featureprop fp, cvterm cvt
    //    WHERE f.feature_id = fp.feature_id AND fp.type_id = cvt.cvterm_id AND
    //      cvt.name = 'promoted_gene_type' AND f.uniquename = 'FBgn0000011';

    /**
     * Return the interesting rows from the libraryprop table.
     * This is a protected method so that it can be overridden for testing
     * @param connection the db connection
     * @return the SQL result set
     * @throws SQLException if a database problem occurs
     */
    protected ResultSet getLibraryFeatureResultSet(Connection connection) throws SQLException {
        String query =
	    "select f.feature_id, lp.value, lp_type.name AS type_name "
	    + "FROM feature f, library_feature lf, library l, libraryprop lp,"
	    + " cvterm lp_type "
	    + " WHERE  f.feature_id=lf.feature_id "
	    + " AND lf.library_id=l.library_id "
	    + " AND l.library_id=lp.library_id "
	    + " AND lp.type_id=lp_type.cvterm_id "
	    + " AND f.feature_id IN (" + getFeatureIdQuery() + ")";
        LOG.info("executing getLibraryFeatureResultSet(): " + query);
        Statement stmt = connection.createStatement();
        ResultSet res = stmt.executeQuery(query);
        return res;
    }

    /**
     * Return the interesting rows from the librarycvterm table.
     * This is a protected method so that it can be overridden for testing
     * @param connection the db connection
     * @return the SQL result set
     * @throws SQLException if a database problem occurs
     */
    protected ResultSet getLibraryCVTermResultSet(Connection connection) throws SQLException {
        String query =
	    "select f.feature_id, d.accession AS term_identifier "
	    + " FROM feature f, library_feature lf, library l, library_cvterm lcvt,"
	    + " cvterm cvt, cv, dbxref d "
	    + " WHERE cv.name IN ('FlyBase anatomy CV','cellular_component') "
	    + " AND lf.library_id=l.library_id AND l.library_id=lcvt.library_id "
	    + " AND lcvt.cvterm_id=cvt.cvterm_id "
	    + " AND cvt.dbxref_id = d.dbxref_id "
	    + " AND f.feature_id IN (" + getFeatureIdQuery() + ")";
        LOG.info("executing getLibraryFeatureResultSet(): " + query);
        Statement stmt = connection.createStatement();
        ResultSet res = stmt.executeQuery(query);
        return res;
    }


    /**
     * Return the interesting rows from the feature_cvterm/cvterm table.  Only returns rows for
     * those features returned by getFeatureIdQuery().
     * This is a protected method so that it can be overridden for testing
     * @param connection the db connection
     * @return the SQL result set
     * @throws SQLException if a database problem occurs
     */
    protected ResultSet getFeatureCVTermResultSet(Connection connection) throws SQLException {
        String query =
	    "SELECT DISTINCT feature_id, cvterm.cvterm_id, cvterm.name AS cvterm_name,"
	    + " cv.name AS cv_name "
	    + " FROM feature_cvterm, cvterm, cv "
	    + " WHERE feature_id IN (" + getFeatureIdQuery() + ")"
	    + " AND cvterm.cvterm_id = feature_cvterm.cvterm_id "
	    + " AND cvterm.cv_id = cv.cv_id "
	    + " ORDER BY feature_id";
        LOG.info("executing getFeatureCVTermResultSet(): " + query);
        Statement stmt = connection.createStatement();
        ResultSet res = stmt.executeQuery(query);
        return res;
    }

    /**
     * Return the interesting rows from the synonym table.
     * This is a protected method so that it can be overridden for testing
     * @param connection the db connection
     * @return the SQL result set
     * @throws SQLException if a database problem occurs
     */
    protected ResultSet getSynonymResultSet(Connection connection) throws SQLException {
        String query =
	    "SELECT DISTINCT feature_id, synonym.name AS synonym_name,"
	    + " cvterm.name AS type_name, is_current"
	    + " FROM feature_synonym, synonym, cvterm"
	    + " WHERE feature_synonym.synonym_id = synonym.synonym_id"
	    + " AND synonym.type_id = cvterm.cvterm_id"
	    + " AND feature_id IN (" + getFeatureIdQuery() + ")"
	    + " ORDER BY is_current DESC";
        LOG.info("executing getSynonymResultSet(): " + query);
        Statement stmt = connection.createStatement();
        ResultSet res = stmt.executeQuery(query);
        return res;
    }

    /**
     * Return the interesting rows from the pub table.
     * This is a protected method so that it can be overridden for testing
     * @param connection the db connection
     * @return the SQL result set
     * @throws SQLException if a database problem occurs
     */
    protected ResultSet getPubResultSet(Connection connection) throws SQLException {
        String query =
	    "SELECT DISTINCT feature_pub.feature_id, dbxref.accession as pub_db_identifier"
	    + " FROM feature_pub, dbxref, db, pub, pub_dbxref"
	    + " WHERE feature_pub.pub_id = pub.pub_id"
	    + " AND pub_dbxref.dbxref_id = dbxref.dbxref_id"
	    + " AND dbxref.db_id = db.db_id"
	    + " AND pub.pub_id = pub_dbxref.pub_id"
	    + " AND db.name = 'pubmed'"
	    + " AND feature_id IN (" + getFeatureIdQuery() + ")"
	    + " ORDER BY feature_pub.feature_id";
        LOG.info("executing getPubResultSet(): " + query);
        Statement stmt = connection.createStatement();
        ResultSet res = stmt.executeQuery(query);
        return res;
    }

    /**
     * Call DataConverter.createSynonym(), store the Item then record in FeatureData that we've
     * created it.
     * @param fdat the FeatureData
     * @param identifier the identifier to store in the Synonym
     * @return the new Synonym
     * @throws ObjectStoreException if there is a problem while storing
     */
    protected Item createSynonym(FeatureData fdat, String identifier)
        throws ObjectStoreException {
        if (fdat.getExistingSynonyms().contains(identifier)) {
            String msg = "feature identifier " + identifier + " is already a synonym for: "
		+ fdat.getExistingSynonyms();
            LOG.info(msg);
            //          TODO:  why would a duplicate synonym require an exception to be thrown?
            //          throw new IllegalArgumentException(msg);
            return null;
        }
        Item returnItem = null;
        try {
            returnItem = getChadoDBConverter().createSynonym(fdat.getItemIdentifier(), identifier,
							     false);
        } catch (ObjectStoreException e) {
            throw new RuntimeException("Couldn't create synonym", e);
        }
        fdat.addExistingSynonym(identifier);
        return returnItem;
    }

    /**
     * Fetch the populated map of chado feature id to FeatureData objects.
     * @return map of feature details
     */
    protected Map<Integer, FeatureData> getFeatureMap() {
        return this.featureMap;
    }

    /**
     * Store and tally debug lines in debugLineMap.
     * @param a string to be output
     */
    private void debugToLog(String output) {
	if (debugLineMap.containsKey(output)) {
	    int count = debugLineMap.get(output).intValue();
	    count++;
	    debugLineMap.put(output, new Integer(count));
	} else {
	    debugLineMap.put(output, new Integer(1));
	}
    }

    /**
     * Dump the debugLineMap and clear it.
     */
    private void logDebugLines() {
	for (String debugLine : debugLineMap.keySet()) {
	    LOG.info(debugLine+" {"+debugLineMap.get(debugLine).intValue()+"}");
	}
	debugLineMap.clear();
    }

}
