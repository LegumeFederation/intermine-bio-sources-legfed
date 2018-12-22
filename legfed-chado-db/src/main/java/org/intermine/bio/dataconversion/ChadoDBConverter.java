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

import java.io.File;
import java.lang.reflect.Constructor;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.intermine.bio.util.OrganismData;
import org.intermine.bio.util.OrganismRepository;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.sql.Database;
import org.intermine.metadata.StringUtil;
import org.intermine.xml.full.Item;

/**
 * DataConverter to read from a Chado database into items
 *
 * @author Kim Rutherford, Sam Hokin
 */
public class ChadoDBConverter extends BioDBConverter {
    
    protected static final Logger LOG = Logger.getLogger(ChadoDBConverter.class);

    // a Map from chado organism_id to OrganismData
    private final Map<Integer, OrganismData> chadoToOrgData = new HashMap<Integer, OrganismData>();

    // a Map from chado organism_id to OrganismData for desired homologue organisms
    private final Map<Integer, OrganismData> chadoToHomologueOrgData = new HashMap<Integer, OrganismData>();

    // a Map from chado organism_id to strain name
    private final Map<Integer,String> chadoToStrainName = new HashMap<Integer,String>();

    // a Map from chado organism_id to strain name for desired homologue strains
    private final Map<Integer,String> chadoToHomologueStrainName = new HashMap<Integer,String>();

    private String processors = "";

    private String reactomeFilename = "";

    private String phytozomeVersion = "";

    private final Set<OrganismData> organismsToProcess = new HashSet<OrganismData>();
    private final Set<String> strainsToProcess = new HashSet<String>();
    
    private final Set<OrganismData> homologueOrganismsToProcess = new HashSet<OrganismData>();
    private final Set<String> homologueStrainsToProcess = new HashSet<String>();
    
    private final OrganismRepository organismRepository;

    private final List<ChadoProcessor> completedProcessors = new ArrayList<ChadoProcessor>();

    private Connection connection;

    /**
     * Create a new ChadoDBConverter object.
     * @param database the database to read from
     * @param tgtModel the Model used by the object store we will write to with the ItemWriter
     * @param writer an ItemWriter used to handle the resultant Items
     * @throws SQLException if we fail to get a database connection

     */
    public ChadoDBConverter(Database database, Model tgtModel, ItemWriter writer) throws SQLException {
        super(database, tgtModel, writer, null, null);
        organismRepository = OrganismRepository.getOrganismRepository();
        if (getDatabase() == null) {
            // no Database when testing and no connection needed
            connection = null;
        } else {
            connection = getDatabase().getConnection();
            LOG.info("Connection schema:"+connection.getSchema());
        }
    }

    /**
     * Set the taxon ids for the desired organisms.
     * from chado with these organisms will be processed.
     *
     * @param organisms a space separated list of the organism taxon ids or abbreviations to look up in the organism table.
     */
    public void setOrganisms(String organisms) {
        String[] bits = StringUtil.split(organisms, " ");
        for (String taxonId: bits) {
            LOG.info("setOrganisms:taxonId="+taxonId);
            OrganismData od = organismRepository.getOrganismDataByTaxon(taxonId);
            if (od == null) {
                throw new RuntimeException("Can't find organism for taxonId " + taxonId);
            }
            LOG.info("od="+od);
            organismsToProcess.add(od);
        }
    }

    /**
     * Set the taxon ids for the desired homologues.
     * @param organisms a space separated list of the taxon ids to look up in the organism table
     */
    public void setHomologueOrganisms(String organisms) {
        String[] bits = StringUtil.split(organisms, " ");
	LOG.info(bits.length+" homologue organisms from:"+organisms);
        for (String taxonId: bits) {
            LOG.info("setHomologueOrganisms:taxonId="+taxonId);
            OrganismData od = organismRepository.getOrganismDataByTaxon(taxonId);
            if (od == null) {
                throw new RuntimeException("Can't find organism for taxonId=" + taxonId);
            }
            homologueOrganismsToProcess.add(od);
        }
    }

    /**
     * Set the strain names for desired organisms.
     *
     * @param strains a space separated list of the strain names.
     */
    public void setStrains(String strains) {
        String[] bits = StringUtil.split(strains, " ");
        for (String strainName : bits) {
            LOG.info("setStrains:strainName="+strainName);
            strainsToProcess.add(strainName);
        }
    }

    /**
     * Set the strain names for desired homologues.
     *
     * @param strains a space separated list of the strain names.
     */
    public void setHomologueStrains(String strains) {
        String[] bits = StringUtil.split(strains, " ");
        for (String strainName : bits) {
            LOG.info("setHomologueStrains:strainName="+strainName);
            homologueStrainsToProcess.add(strainName);
        }
    }

    /**
     * Set the class names of the ChadoProcessors to run.
     * @param processors a space separated list of the fully-qualified class names of module
     * processors to run
     */
    public void setProcessors(String processors) {
        this.processors = processors;
    }

    /**
     * Set the name of the reactome file to run in ReactomeProcessor
     * @param filename the full file name of the tab-delimited reactome file
     */
    public void setReactomeFilename(String reactomeFilename) {
        this.reactomeFilename = reactomeFilename;
    }

    /**
     * Return the name of the tab-delimited reactome file
     */
    public String getReactomeFilename() {
        return reactomeFilename;
    }

    /**
     * Set the Phytozome version for querying the phylotree table.
     * @param phytozomeVersion the Phytozome version, e.g. phytozome_10_2
     */
    public void setPhytozomeVersion(String phytozomeVersion) {
        this.phytozomeVersion = phytozomeVersion;
    }

    /**
     * Return the Phytozome version
     */
    public String getPhytozomeVersion() {
        return phytozomeVersion;
    }
    
    /**
     * Return a map from chado organism_id to OrganismData object for all the organisms that we are processing.
     * @return the Map
     */
    public Map<Integer, OrganismData> getChadoIdToOrgDataMap() {
        return chadoToOrgData;
    }

    /**
     * Return a map from chado organism_id to OrganismData object for the organisms that we are processing for homology.
     * @return the Map
     */
    public Map<Integer, OrganismData> getChadoIdToHomologueOrgDataMap() {
        return chadoToHomologueOrgData;
    }

    /**
     * Return a map from chado organism_id to strain name.
     * @ return the strain name
     */
    public Map<Integer,String> getChadoIdToStrainNameMap() {
        return chadoToStrainName;
    }

    /**
     * Return a map from chado organism_id to strain name for homology organisms.
     * @ return the strain name
     */
    public Map<Integer,String> getChadoIdToHomologueStrainNameMap() {
        return chadoToHomologueStrainName;
    }

    /**
     * Get the connection to use when processing.
     * @return the Connection, or null while testing
     */
    protected Connection getConnection() {
        return connection;
    }

    /**
     * Process the data from the Database and write to the ItemWriter.
     * {@inheritDoc}
     */
    @Override
    public void process() throws Exception {

        if (StringUtils.isEmpty(processors)) {
            throw new IllegalArgumentException("processors not set in ChadoDBConverter");
        }

        Map<Integer,OrganismData> tempChadoOrgMap = getChadoOrganismIds(getConnection());
        Map<Integer,String> tempChadoStrainMap = getChadoStrainNames(getConnection());

        // DEBUG
        for (OrganismData od : organismsToProcess) {
            LOG.info("organismsToProcess: ["+od+"]");
        }
        for (Integer chadoId : tempChadoOrgMap.keySet()) {
            LOG.info("tempChadoOrgMap: key="+chadoId+" value=["+tempChadoOrgMap.get(chadoId)+"]");
        }
        for (String strainName : strainsToProcess) {
            LOG.info("strainsToProcess: "+strainName);
        }
        for (Integer chadoId : tempChadoStrainMap.keySet()) {
            LOG.info("tempChadoStrainMap: key="+chadoId+" value="+tempChadoStrainMap.get(chadoId));
        }

        // build the map of chadoId to desired organisms (and homologue organisms) to process
        for (Integer chadoId : tempChadoOrgMap.keySet()) {
            OrganismData od = tempChadoOrgMap.get(chadoId);
            if (organismsToProcess.contains(od)) {
                chadoToOrgData.put(chadoId, od);
            }
            if (homologueOrganismsToProcess.contains(od)) {
                chadoToHomologueOrgData.put(chadoId, od);
            }
        }

        // build the map of chadoId to desired strains (and homologue strains) to process
        for (Integer chadoId : tempChadoStrainMap.keySet()) {
            String strainName = tempChadoStrainMap.get(chadoId);
            if (strainsToProcess.contains(strainName)) {
                chadoToStrainName.put(chadoId, strainName);
            }
            if (homologueStrainsToProcess.contains(strainName)) {
                chadoToHomologueStrainName.put(chadoId, strainName);
            }
        }

        String[] bits = processors.trim().split("[ \\t]+");
        for (int i = 0; i < bits.length; i++) {
            String className = bits[i];
            if (!StringUtils.isEmpty(className)) {
                Class<?> cls = Class.forName(className);
                Constructor<?> constructor = cls.getDeclaredConstructor(ChadoDBConverter.class);
                ChadoProcessor currentProcessor = (ChadoProcessor) constructor.newInstance(this);
                currentProcessor.process(getConnection());
                getCompletedProcessors().add(currentProcessor);
            }
        }
    }

    /**
     * Return a map from chado organism_id to OrganismData.
     * @param conn the db connection
     * @return a Map from abbreviation to chado organism_id
     * @throws SQLException if there is a database problem
     */
    protected Map<Integer,OrganismData> getChadoOrganismIds(Connection conn) throws SQLException {
        String query = "SELECT organism_id, abbreviation, genus, species FROM organism";
        LOG.info("Executing: " + query);
        Statement stmt = conn.createStatement();
        ResultSet res = stmt.executeQuery(query);
        LOG.info("ResultSet returned.");
        Map<Integer,OrganismData> retMap = new HashMap<>();
        OrganismRepository or = OrganismRepository.getOrganismRepository();
        LOG.info("OrganismRepository gotten.");
        while (res.next()) {
            int organismId = res.getInt("organism_id");
            String abbreviation = res.getString("abbreviation");
            String genus = res.getString("genus");
            String species = res.getString("species");
            String strain = null;
            // strains are signified with underscore on species in chado, e.g. arietinum_desi, arietinum_kabuli
            if (genus!=null && species!=null && species.contains("_")) {
                String[] parts = species.split("_");
                species = parts[0];
            }
            LOG.info(organismId+":"+abbreviation+":"+genus+":"+species);
            
            // use genus and species to get Taxon ID
            OrganismData od1 = or.getOrganismDataByGenusSpecies(genus, species);
            if (od1==null) {
                throw new RuntimeException("Could not get OrganismData from genus,species:"+genus+","+species);
            }
            String taxonId = od1.getTaxonId();
            
            // get full OrganismData with taxonId
            OrganismData od = or.getOrganismDataByTaxon(taxonId);
            
            retMap.put(new Integer(organismId), od);
        }
        return retMap;
    }

    /**
     * Return a map from chado organism_id to strain name.
     * @param conn the db connection
     * @return a Map from abbreviation to chado organism_id
     * @throws SQLException if the is a database problem
     */
    protected Map<Integer,String> getChadoStrainNames(Connection conn) throws SQLException {
        String query = "SELECT organism_id, abbreviation, genus, species FROM organism";
        LOG.info("Executing: " + query);
        Statement stmt = conn.createStatement();
        ResultSet res = stmt.executeQuery(query);
        LOG.info("ResultSet returned.");
        Map<Integer,String> retMap = new HashMap<>();
        while (res.next()) {
            int organismId = res.getInt("organism_id");
            String abbreviation = res.getString("abbreviation");
            String genus = res.getString("genus");
            String species = res.getString("species");
            String strain = null;
            // strains are signified with underscore on species in chado, e.g. arietinum_desi, arietinum_kabuli
            if (genus!=null && species!=null && species.contains("_")) {
                String[] parts = species.split("_");
                species = parts[0];
                strain = parts[1];
                retMap.put(new Integer(organismId), strain);
                LOG.info(organismId+":"+abbreviation+":"+genus+":"+species+":"+strain);
            }
        }
        return retMap;
    }

    /**
     * Return the OrganismData objects for the organisms listed in the source configuration.
     * @return the organismsToProcess
     */
    public Set<OrganismData> getOrganismsToProcess() {
        return organismsToProcess;
    }

    /**
     * Return the names of strains to process.
     * @return the strains to process
     */
    public Set<String> getStrainsToProcess() {
        return strainsToProcess;
    }

    /**
     * Return the OrganismData objects for the homologue organisms listed in the source configuration.
     * @return homologueOrganismsToProcess
     */
    public Set<OrganismData> getHomologueOrganismsToProcess() {
        return homologueOrganismsToProcess;
    }

    /**
     * Look at the list of completed processors and return the processor of the given type.  If
     * there is none or more than one, throw a RuntimeException
     * @param cls the class
     * @return the ChadoProcessor
     */
    public ChadoProcessor findProcessor(Class<? extends ChadoProcessor> cls) {
        ChadoProcessor returnProcessor = null;

        for (ChadoProcessor processor: getCompletedProcessors()) {
            if (cls.isAssignableFrom(processor.getClass())) {
                if (returnProcessor == null) {
                    returnProcessor = processor;
                } else {
                    throw new RuntimeException("Completed processors list contains two objects of "
                                               + "type: " + cls.getName());
                }
            }
        }

        if (returnProcessor == null) {
            throw new RuntimeException("Can't find `" + cls.getName() + "` before `"
                                       + this.getClass().getName()
                                       + "` in the list of completed processors - must run "
                                       + cls.getName() + " first.");
        }
        return returnProcessor;
    }
    
    /**
     * Default implementation that makes a data set title based on the data source name.
     * {@inheritDoc}
     */
    @Override
    public String getDataSetTitle(String taxonId) {
        OrganismData od = organismRepository.getOrganismDataByTaxon(taxonId);
        if (od != null) {
            return getDataSourceName() + " data set for " + od.getGenus() + " " + od.getSpecies();
        }
        return getDataSourceName() + " data set";
    }

    /**
     * @return the completedProcessors
     */
    public List<ChadoProcessor> getCompletedProcessors() {
        return completedProcessors;
    }

    /**
     * @return a Set of desired chado organism_ids
     */
    public Set<Integer> getDesiredChadoOrganismIds() {
        Set<Integer> organismIds = new HashSet<>();
        for (Integer organismId : chadoToOrgData.keySet()) {
            organismIds.add(organismId);
        }
        for (Integer organismId : chadoToStrainName.keySet()) {
            organismIds.add(organismId);
        }
        return organismIds;
    }

    /**
     * @return a Set of desired chado homologue organism_ids
     */
    public Set<Integer> getDesiredChadoHomologueOrganismIds() {
        Set<Integer> organismIds = new HashSet<>();
        for (Integer organismId : chadoToHomologueOrgData.keySet()) {
            organismIds.add(organismId);
        }
        for (Integer organismId : chadoToHomologueStrainName.keySet()) {
            organismIds.add(organismId);
        }
        return organismIds;
    }

    /**
     * Provides the Organism item for the given chado organism_id.
     * @param chadoId the chado organism_id
     * @return the Organism Item
     */
    public Item getOrganismItem(Integer chadoId) {
	OrganismData od = chadoToOrgData.get(chadoId);
	String taxonId = od.getTaxonId();
	return getOrganismItem(taxonId);
    }

    /**
     * Provides the Strain item for the given chado organism_id.
     * @param chadoId the chado organism_id
     * @return the Strain Item
     */
    public Item getStrainItem(Integer chadoId) {
	OrganismData od = chadoToOrgData.get(chadoId);
	String taxonId = od.getTaxonId();
	String strainName = chadoToStrainName.get(chadoId);
	return getStrainItem(strainName, taxonId);
    }

    /**
     * Provides the Organism item for the given homologue chado organism_id.
     * @param chadoId the chado organism_id
     * @return the Organism Item
     */
    public Item getHomologueOrganismItem(Integer chadoId) {
	OrganismData od = chadoToHomologueOrgData.get(chadoId);
	String taxonId = od.getTaxonId();
	return getOrganismItem(taxonId);
    }

    /**
     * Provides the Strain item for the given homologue chado organism_id.
     * @param chadoId the chado organism_id
     * @return the Strain Item
     */
    public Item getHomologueStrainItem(Integer chadoId) {
	OrganismData od = chadoToHomologueOrgData.get(chadoId);
	String taxonId = od.getTaxonId();
	String strainName = chadoToHomologueStrainName.get(chadoId);
	return getStrainItem(strainName, taxonId);
    }

    /**
     * Provides the Strain.primaryIdentifier for the given chado organism_id.
     * @param chadoId the chado organism_id
     * @return the Strain.primaryIdentifier
     */
    public String getStrainName(Integer chadoId) {
	OrganismData od = chadoToOrgData.get(chadoId);
	return chadoToStrainName.get(chadoId);
    }

}
