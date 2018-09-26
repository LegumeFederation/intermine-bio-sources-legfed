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

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.Map;
import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;

import org.apache.log4j.Logger;

import org.intermine.bio.util.OrganismData;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Attribute;
import org.intermine.xml.full.Item;
import org.intermine.xml.full.Reference;

/**
 * Create and store Protein, ProteinMatch, ProteinHmmMatch, ConsensusRegion, ProteinDomain and mRNA objects by querying the chado feature and featureloc tables.
 * (The mRNAs are fully loaded as sequence features in SequenceProcessor, so here we simply associate them with proteins.)
 *
 * Since this processor deals only with chado data, Items are stored in maps with Integer keys equal to
 * the chado feature.feature_id.
 *
 * project.xml parameters:
 *   taxonid_variety e.g. "3920_IT97K-499-35"
 *
 * @author Sam Hokin, NCGR
 */
public class ProteinProcessor extends ChadoProcessor {
	
    private static final Logger LOG = Logger.getLogger(ProteinProcessor.class);

    /**
     * Create a new ProteinProcessor
     *
     * @param chadoDBConverter the ChadoDBConverter that is controlling this processor
     */
    public ProteinProcessor(ChadoDBConverter chadoDBConverter) {
        super(chadoDBConverter);
    }

    /**
     * {@inheritDoc}
     * We process the chado database by reading the feature and featureloc tables.
     */
    @Override
    public void process(Connection connection) throws SQLException, ObjectStoreException {
        
        // initialize our DB stuff
        Statement stmt1 = connection.createStatement();
        Statement stmt2 = connection.createStatement();
        Statement stmt3 = connection.createStatement();
        ResultSet rs1;
        ResultSet rs2;
        ResultSet rs3;
        
        // get and store the desired organisms from the supplied source taxon IDs
        Map<Integer,Item> organismMap = new HashMap<Integer,Item>();
        Map<Integer,OrganismData> chadoToOrgData = getChadoDBConverter().getChadoIdToOrgDataMap();
        for (Integer organismId : chadoToOrgData.keySet()) {
            OrganismData organismData = chadoToOrgData.get(organismId);
            String taxonId = organismData.getTaxonId();
            String variety = organismData.getVariety();
            Item organism = getChadoDBConverter().createItem("Organism");
            organism.setAttribute("taxonId", String.valueOf(taxonId));
            organism.setAttribute("variety", variety); // required
            store(organism);
	    organismMap.put(organismId, organism);
        }
        LOG.info("Created "+organismMap.size()+" organism Items.");
        if (organismMap.size()==0) {
            throw new RuntimeException("Property organisms must contain at least one taxonID_variety in project.xml.");
        }

        // CV term IDs
        int proteinTypeId = getCVTermId(stmt1, "polypeptide");
        int proteinDomainTypeId = getCVTermId(stmt1, "polypeptide_domain");
        int proteinMatchTypeId = getCVTermId(stmt1, "protein_match");
        int proteinHmmMatchTypeId = getCVTermId(stmt1, "protein_hmm_match");
        int consensusRegionTypeId = getCVTermId(stmt1, "consensus_region");
        int mRNATypeId = getCVTermId(stmt1, "mRNA");

        // consensus regions are organism-independent and are associated with multiple proteins / protein HMM matches
        Map<Integer,Item> crMap = new HashMap<Integer,Item>();

        // protein domains can be associated with multiple proteins
        Map<Integer,Item> proteinDomainMap = new HashMap<Integer,Item>();

        // cycle through organisms
        for (Integer organismId : organismMap.keySet()) {
            Item organism = organismMap.get(organismId);
            // query the proteins for this organism
            rs1 = stmt1.executeQuery("SELECT * FROM feature WHERE organism_id="+organismId+" AND type_id="+proteinTypeId);
            while (rs1.next()) {

                // create the Protein item
                int proteinId = rs1.getInt("feature_id");
                String proteinUniqueName = rs1.getString("uniquename");
                String proteinName = rs1.getString("name");
                String proteinResidues = rs1.getString("residues");
                int proteinLength = rs1.getInt("seqlen");
                String proteinMd5Checksum = rs1.getString("md5checksum");
                Item protein = getChadoDBConverter().createItem("Protein");
                protein.setAttribute("primaryIdentifier", proteinUniqueName);
                protein.setAttribute("secondaryIdentifier", proteinName);
                protein.setAttribute("length", String.valueOf(proteinLength));
                protein.setAttribute("chadoFeatureId", String.valueOf(proteinId));
                protein.setAttribute("chadoUniqueName", proteinUniqueName);
                protein.setAttribute("chadoName", proteinName);
                if (proteinMd5Checksum!=null) protein.setAttribute("md5checksum", proteinMd5Checksum);
                protein.setReference("organism", organism);

                // create and store the protein Sequence item
                Item proteinSequence = getChadoDBConverter().createItem("Sequence");
                proteinSequence.setAttribute("length", String.valueOf(proteinLength));
                if (proteinMd5Checksum!=null) proteinSequence.setAttribute("md5checksum", proteinMd5Checksum);
                proteinSequence.setAttribute("residues", proteinResidues);
                store(proteinSequence);
                protein.setReference("sequence", proteinSequence);

                // query, create and store the consensus region(s) (without their sequences, which we'll add in GeneFamilyProcessor)
                rs2 = stmt2.executeQuery("SELECT feature.* FROM feature,featureloc WHERE type_id="+consensusRegionTypeId+
                                         " AND feature.feature_id=featureloc.srcfeature_id AND featureloc.feature_id="+proteinId);
                while (rs2.next()) {
                    int crId = rs2.getInt("feature_id");
                    String crUniqueName = rs2.getString("uniquename");
                    String crName = rs2.getString("name");
                    String crResidues = rs2.getString("residues");
                    int crLength = rs2.getInt("seqlen");
                    String crMd5Checksum = rs2.getString("md5checksum");
                    Item cr;
                    if (crMap.containsKey(crId)) {
                        cr = crMap.get(crId);
                    } else {
                        cr = getChadoDBConverter().createItem("ConsensusRegion");
                        cr.setAttribute("primaryIdentifier", crUniqueName);
                        cr.setAttribute("secondaryIdentifier", crName);
                        cr.setAttribute("length", String.valueOf(crLength));
                        cr.setAttribute("chadoFeatureId", String.valueOf(crId));
                        cr.setAttribute("chadoUniqueName", crUniqueName);
                        cr.setAttribute("chadoName", crName);
                        if (crMd5Checksum!=null) cr.setAttribute("md5checksum", crMd5Checksum);
                        // create and store the consensus region Sequence item
                        Item crSequence = getChadoDBConverter().createItem("Sequence");
                        crSequence.setAttribute("length", String.valueOf(crLength));
                        if (crMd5Checksum!=null) crSequence.setAttribute("md5checksum", crMd5Checksum);
                        crSequence.setAttribute("residues", crResidues);
                        store(crSequence);
                        cr.setReference("sequence", crSequence);
                        store(cr);
                        crMap.put(crId, cr);
                    }
                    protein.addToCollection("consensusRegions", cr);
                }
                rs2.close();

                // query, create and store the protein matches to this protein
                rs2 = stmt2.executeQuery("SELECT feature.*,fmin,fmax FROM feature,featureloc WHERE type_id="+proteinMatchTypeId+
                                         " AND feature.feature_id=featureloc.feature_id AND featureloc.srcfeature_id="+proteinId);
                while (rs2.next()) {
                    int proteinMatchId = rs2.getInt("feature_id");
                    String proteinMatchUniqueName = rs2.getString("uniquename");
                    String proteinMatchName = rs2.getString("name");
                    int start = rs2.getInt("fmin") + 1; // zero-based
                    int end = rs2.getInt("fmax");
                    Item proteinMatch = getChadoDBConverter().createItem("ProteinMatch");
                    proteinMatch.setAttribute("primaryIdentifier", proteinMatchUniqueName);
                    proteinMatch.setAttribute("secondaryIdentifier", proteinMatchName);
                    proteinMatch.setAttribute("chadoFeatureId", String.valueOf(proteinMatchId));
                    proteinMatch.setAttribute("chadoUniqueName", proteinMatchUniqueName);
                    proteinMatch.setAttribute("chadoName", proteinMatchName);
                    proteinMatch.setReference("organism", organism);
                    proteinMatch.setReference("protein", protein);
                    Item proteinLocation = getChadoDBConverter().createItem("Location");
                    proteinLocation.setAttribute("start", String.valueOf(start));
                    proteinLocation.setAttribute("end", String.valueOf(end));
                    proteinLocation.setReference("locatedOn", protein);
                    proteinLocation.setReference("feature", proteinMatch);
                    store(proteinLocation);
                    proteinMatch.setReference("proteinLocation", proteinLocation);
                    store(proteinMatch);
                }
                rs2.close();
                
                // query, create and store the protein HMM matches to this protein
                rs2 = stmt2.executeQuery("SELECT feature.*,fmin,fmax FROM feature,featureloc WHERE type_id="+proteinHmmMatchTypeId+
                                         " AND feature.feature_id=featureloc.feature_id AND featureloc.srcfeature_id="+proteinId);
                while (rs2.next()) {
                    int proteinHmmMatchId = rs2.getInt("feature_id");
                    String proteinHmmMatchUniqueName = rs2.getString("uniquename");
                    String proteinHmmMatchName = rs2.getString("name");
                    int pStart = rs2.getInt("fmin") + 1; // zero-based
                    int pEnd = rs2.getInt("fmax");
                    Item proteinHmmMatch = getChadoDBConverter().createItem("ProteinHmmMatch");
                    proteinHmmMatch.setAttribute("primaryIdentifier", proteinHmmMatchUniqueName);
                    proteinHmmMatch.setAttribute("secondaryIdentifier", proteinHmmMatchName);
                    proteinHmmMatch.setAttribute("chadoFeatureId", String.valueOf(proteinHmmMatchId));
                    proteinHmmMatch.setAttribute("chadoUniqueName", proteinHmmMatchUniqueName);
                    proteinHmmMatch.setAttribute("chadoName", proteinHmmMatchName);
                    proteinHmmMatch.setReference("organism", organism);
                    proteinHmmMatch.setReference("protein", protein);
                    // locate the HMM on the protein
                    Item proteinLocation = getChadoDBConverter().createItem("Location");
                    proteinLocation.setAttribute("start", String.valueOf(pStart));
                    proteinLocation.setAttribute("end", String.valueOf(pEnd));
                    proteinLocation.setReference("feature", proteinHmmMatch);
                    proteinLocation.setReference("locatedOn", protein);
                    store(proteinLocation);
                    proteinHmmMatch.setReference("proteinLocation", proteinLocation);
                    // query, create and store the ProteinDomain associated with this ProteinHmmMatch (and therefore this Protein)
                    rs3 = stmt3.executeQuery("SELECT feature.*,fmin,fmax FROM feature,featureloc WHERE type_id="+proteinDomainTypeId+
                                             " AND feature.feature_id=featureloc.srcfeature_id AND featureloc.feature_id="+proteinHmmMatchId);
                    while (rs3.next()) {
                        int proteinDomainId = rs3.getInt("feature_id");
                        String proteinDomainUniqueName = rs3.getString("uniquename");
                        String proteinDomainName = rs3.getString("name");
                        int pdStart = rs3.getInt("fmin") + 1; // zero-based
                        int pdEnd = rs3.getInt("fmax");
                        Item proteinDomain;
                        if (proteinDomainMap.containsKey(proteinDomainId)) {
                            proteinDomain = proteinDomainMap.get(proteinDomainId);
                        } else {
                            proteinDomain = getChadoDBConverter().createItem("ProteinDomain");
                            proteinDomain.setAttribute("primaryIdentifier", proteinDomainUniqueName);
                            proteinDomain.setAttribute("secondaryIdentifier", proteinDomainName);
                            proteinDomain.setAttribute("chadoFeatureId", String.valueOf(proteinDomainId));
                            proteinDomain.setAttribute("chadoUniqueName", proteinDomainUniqueName);
                            proteinDomain.setAttribute("chadoName", proteinDomainName);
                            store(proteinDomain);
                            proteinDomainMap.put(proteinDomainId, proteinDomain);
                        }
                        proteinHmmMatch.setReference("proteinDomain", proteinDomain);
                        // locate the HMM on the protein domain
                        Item proteinDomainLocation = getChadoDBConverter().createItem("Location");
                        proteinDomainLocation.setAttribute("start", String.valueOf(pdStart));
                        proteinDomainLocation.setAttribute("end", String.valueOf(pdEnd));
                        proteinDomainLocation.setReference("feature", proteinHmmMatch);
                        proteinDomainLocation.setReference("locatedOn", proteinDomain);
                        store(proteinDomainLocation);
                        proteinHmmMatch.setReference("proteinDomainLocation", proteinDomainLocation);
                        protein.addToCollection("proteinDomains", proteinDomain);
                    }
                    rs3.close();
                    store(proteinHmmMatch);
                }
                rs2.close();

                // store the protein
                store(protein);
            }
            rs1.close();

        }

    }

    /**
     * Store the item.
     * @param item the Item
     * @return the database id of the new Item
     * @throws ObjectStoreException if an error occurs while storing
     */
    protected Integer store(Item item) throws ObjectStoreException {
        return getChadoDBConverter().store(item);
    }

    /**
     * Get the CVTerm ID for a given CVTerm name.
     * @param stmt the database connection statement, initialized to the chado database
     * @param name the desired CV term name
     * @return the CV term id
     * @throws SQLException
     */
     protected int getCVTermId(Statement stmt, String name) throws SQLException {
        ResultSet rs = stmt.executeQuery("SELECT cvterm_id FROM cvterm WHERE name='"+name+"'");
        if (rs.next()) {
            int cvtermId = rs.getInt("cvterm_id");
            rs.close();
            return cvtermId;
        } else {
            throw new RuntimeException("Could not determine CV term id for '"+name+"'.");
        }
     }
    
}
