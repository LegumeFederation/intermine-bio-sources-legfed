package org.intermine.bio.postprocess;

/*
 * Copyright (C) 2002-2016 FlyMine, Legume Federation
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */

import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import org.intermine.bio.util.Constants;
import org.intermine.metadata.ConstraintOp;
import org.intermine.model.bio.Chromosome;
import org.intermine.model.bio.Exon;
import org.intermine.model.bio.Gene;
import org.intermine.model.bio.GeneticMarker;
import org.intermine.model.bio.Location;
import org.intermine.model.bio.MRNA;
import org.intermine.model.bio.Polypeptide;
import org.intermine.model.bio.QTL;
import org.intermine.objectstore.ObjectStore;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.objectstore.ObjectStoreWriter;
import org.intermine.objectstore.intermine.ObjectStoreInterMineImpl;
import org.intermine.objectstore.query.ConstraintSet;
import org.intermine.objectstore.query.ContainsConstraint;
import org.intermine.objectstore.query.Query;
import org.intermine.objectstore.query.QueryClass;
import org.intermine.objectstore.query.QueryCollectionReference;
import org.intermine.objectstore.query.QueryObjectReference;
import org.intermine.objectstore.query.Results;
import org.intermine.objectstore.query.ResultsRow;
import org.intermine.postprocess.PostProcessor;

import org.apache.log4j.Logger;

/**
 * 0. Check that exons and mRNAs all have geneid - this detects when merging has replaced genes.
 * 1. Find "spanned genes" for each range of markers for a QTL. Only include QTLs that have at least two associated markers (usually three: nearest, flanking left and right).
 * 2. Relate polypeptides to genes via polypeptide.mRNA:mRNA.gene.
 *
 * Since the QTLs and genetic markers can be loaded from a variety of sources (flat files, chado), it makes sense to do this in post-processing when the 
 * QTLs, markers and genes exist in the database in a standard format. 
 *
 * The polypeptide-gene connection used to be in SequenceProcessor, but it makes more sense to do it as a post-processing task since it doesn't involve importing new data.
 *
 * @author Sam Hokin
 */
public class LegfedPostProcess extends PostProcessor {

    private static final Logger LOG = Logger.getLogger(LegfedPostProcess.class);

    /**
     * Create a new instance of LegfedPostProcess
     *
     * @param osw object store writer
-     */
    public LegfedPostProcess(ObjectStoreWriter osw) {
        super(osw);
    }

    /**
     * {@inheritDoc}
     * Main post-processing routine. Find the "spanned genes" for the QTL+markers.
     * @throws ObjectStoreException if the objectstore throws an exception
     */
    public void postProcess() throws ObjectStoreException {

        // ----------------------------------------------------------------------------------------------
        // 0 - Integrity checks
        // ----------------------------------------------------------------------------------------------

        LOG.info("Checking that all exons and mRNAs have gene IDs...");

        // query all exons
        Query qExon = new Query();
        qExon.setDistinct(false);
        QueryClass qcExon = new QueryClass(Exon.class);
        qExon.addFrom(qcExon);
        qExon.addToSelect(qcExon);
        qExon.addToOrderBy(qcExon);

        // execute the query
        Results exonResults = osw.getObjectStore().execute(qExon);
        Iterator<?> exonIter = exonResults.iterator();
        while (exonIter.hasNext()) {
            try {
                ResultsRow<?> rr = (ResultsRow<?>) exonIter.next();
                Exon exon = (Exon) rr.get(0);
                String primaryIdentifier = (String) exon.getFieldValue("primaryIdentifier");
                Gene gene = (Gene) exon.getFieldValue("gene");
                if (gene==null) LOG.error("Exon "+primaryIdentifier+" has no related gene.");
            } catch (IllegalAccessException ex) {
                throw new RuntimeException(ex);
            }
        }

        // query all mRNAs
        Query qMRNA = new Query();
        qMRNA.setDistinct(false);
        QueryClass qcMRNA = new QueryClass(MRNA.class);
        qMRNA.addFrom(qcMRNA);
        qMRNA.addToSelect(qcMRNA);
        qMRNA.addToOrderBy(qcMRNA);

        // execute the query
        Results mRNAResults = osw.getObjectStore().execute(qMRNA);
        Iterator<?> mRNAIter = mRNAResults.iterator();
        while (mRNAIter.hasNext()) {
            try {
                ResultsRow<?> rr = (ResultsRow<?>) mRNAIter.next();
                MRNA mRNA = (MRNA) rr.get(0);
                String primaryIdentifier = (String) mRNA.getFieldValue("primaryIdentifier");
                Gene gene = (Gene) mRNA.getFieldValue("gene");
                if (gene==null) LOG.error("mRNA "+primaryIdentifier+" has no related gene.");
            } catch (IllegalAccessException ex) {
                throw new RuntimeException(ex);
            }
        }
        
        LOG.info("Integrity check done.");
        
        // ----------------------------------------------------------------------------------------------
        // 1 - First section: accumulate the QTLs and their genomic spans, from their associated markers
        // ----------------------------------------------------------------------------------------------

        LOG.info("Accumulate QTLs and genomic spans from associated markers...");
        
        Query qQTL = new Query();
        qQTL.setDistinct(false);
        ConstraintSet csQTL = new ConstraintSet(ConstraintOp.AND);

        // 0 QTL
        QueryClass qcQTL = new QueryClass(QTL.class);
        qQTL.addFrom(qcQTL);
        qQTL.addToSelect(qcQTL);
        qQTL.addToOrderBy(qcQTL);

        // 1 QTL.associatedGeneticMarkers
        QueryClass qcGeneticMarker = new QueryClass(GeneticMarker.class);
        qQTL.addFrom(qcGeneticMarker);
        qQTL.addToSelect(qcGeneticMarker);
        qQTL.addToOrderBy(qcGeneticMarker);
        QueryCollectionReference qtlGeneticMarkers = new QueryCollectionReference(qcQTL, "associatedGeneticMarkers");
        csQTL.addConstraint(new ContainsConstraint(qtlGeneticMarkers, ConstraintOp.CONTAINS, qcGeneticMarker));

        // 2 GeneticMarker.chromosome
        QueryClass qcChromosome = new QueryClass(Chromosome.class);
        qQTL.addFrom(qcChromosome);
        qQTL.addToSelect(qcChromosome);
        qQTL.addToOrderBy(qcChromosome);
        QueryObjectReference gmChromosome = new QueryObjectReference(qcGeneticMarker, "chromosome");
        csQTL.addConstraint(new ContainsConstraint(gmChromosome, ConstraintOp.CONTAINS, qcChromosome));

        // 3 GeneticMarker.chromosomeLocation
        QueryClass qcLocation = new QueryClass(Location.class);
        qQTL.addFrom(qcLocation);
        qQTL.addToSelect(qcLocation);
        qQTL.addToOrderBy(qcLocation);
        QueryObjectReference gmLocation = new QueryObjectReference(qcGeneticMarker, "chromosomeLocation");
        csQTL.addConstraint(new ContainsConstraint(gmLocation, ConstraintOp.CONTAINS, qcLocation));

        // set the constraints
        qQTL.setConstraint(csQTL);

        // execute the query
        Results qtlResults = osw.getObjectStore().execute(qQTL);
        Iterator<?> qtlIter = qtlResults.iterator();

        // we'll store our QTLs and genomic span in a set for comparison when we drill through the genes
        Set<QTLSpan> qtlSpanSet = new HashSet<QTLSpan>();
        
        QTL lastQTL = null;
        String lastQTLId = "";
        String lastChrId = null;
        String startMarkerId = null;
        String endMarkerId = null;
        int minStart = 100000000;
        int maxEnd = 0;
        int count = 0;
        
        while (qtlIter.hasNext()) {
            try {
                ResultsRow<?> rr = (ResultsRow<?>) qtlIter.next();
                // objects
                QTL qtl = (QTL) rr.get(0);
                GeneticMarker gm = (GeneticMarker) rr.get(1);
                Chromosome chr = (Chromosome) rr.get(2);
                Location loc = (Location) rr.get(3);
                // field values
                String qtlId = (String) qtl.getFieldValue("primaryIdentifier");
                String chrId = (String) chr.getFieldValue("primaryIdentifier");
                String gmId = (String) gm.getFieldValue("primaryIdentifier");
                int start = ((Integer) loc.getFieldValue("start")).intValue();
                int end = ((Integer) loc.getFieldValue("end")).intValue();
                // logic
                if (qtlId.equals(lastQTLId)) {
                    count++;
                    if (start<minStart) {
                        minStart = start;
                        startMarkerId = gmId;
                    }
                    if (end>maxEnd) {
                        maxEnd = end;
                        endMarkerId = gmId;
                    }
                } else {
                    // only store QTLs with more than one marker, so there is an actual genomic range defined
                    if (maxEnd>0 && count>1) {
                        LOG.info("QTL "+lastQTLId+" has "+count+" markers; "+startMarkerId+" to "+endMarkerId+" span ("+minStart+","+maxEnd+") on "+lastChrId);
                        qtlSpanSet.add(new QTLSpan(lastQTL, lastChrId, minStart, maxEnd));
                    }
                    lastQTL = qtl;
                    lastQTLId = qtlId;
                    lastChrId = chrId;
                    minStart = start;
                    maxEnd = end;
                    startMarkerId = gmId;
                    endMarkerId = gmId;
                    count = 1;
                }
            } catch (IllegalAccessException ex) {
                throw new RuntimeException(ex);
            }
        }
        // last one
        qtlSpanSet.add(new QTLSpan(lastQTL, lastChrId, minStart, maxEnd));

        // ------------------------------------------------------------------------------------------------------------------------------
        // 1 - Second section: spin through the genes, comparing their genomic range to the QTL genomic spans, associate them if overlapping
        // ------------------------------------------------------------------------------------------------------------------------------

        LOG.info("Spin through genes, associate with QTL if spanned by QTL genomic range...");

        Query qGene = new Query();
        qGene.setDistinct(false);
        ConstraintSet csGene = new ConstraintSet(ConstraintOp.AND);

        // 0 Gene
        QueryClass qcGene = new QueryClass(Gene.class);
        qGene.addFrom(qcGene);
        qGene.addToSelect(qcGene);
        qGene.addToOrderBy(qcGene);

        // 1 Gene.chromosome
        QueryClass qcGeneChromosome = new QueryClass(Chromosome.class);
        qGene.addFrom(qcGeneChromosome);
        qGene.addToSelect(qcGeneChromosome);
        qGene.addToOrderBy(qcGeneChromosome);
        QueryObjectReference geneChromosome = new QueryObjectReference(qcGene, "chromosome");
        csGene.addConstraint(new ContainsConstraint(geneChromosome, ConstraintOp.CONTAINS, qcGeneChromosome));

        // 2 Gene.chromosomeLocation
        QueryClass qcGeneLocation = new QueryClass(Location.class);
        qGene.addFrom(qcGeneLocation);
        qGene.addToSelect(qcGeneLocation);
        qGene.addToOrderBy(qcGeneLocation);
        QueryObjectReference geneLocation = new QueryObjectReference(qcGene, "chromosomeLocation");
        csGene.addConstraint(new ContainsConstraint(geneLocation, ConstraintOp.CONTAINS, qcGeneLocation));

        // set the overall constraint
        qGene.setConstraint(csGene);

        // execute the query
        Results geneResults = osw.getObjectStore().execute(qGene);
        Iterator<?> geneIter = geneResults.iterator();

        // begin transaction
        osw.beginTransaction();
        
        while (geneIter.hasNext()) {
            try {
                ResultsRow<?> rr = (ResultsRow<?>) geneIter.next();
                // objects
                Gene gene = (Gene) rr.get(0);
                Chromosome chr = (Chromosome) rr.get(1);
                Location loc = (Location) rr.get(2);
                // field values
                String geneId = (String) gene.getFieldValue("primaryIdentifier");
                String chrId = (String) chr.getFieldValue("primaryIdentifier");
                int start = ((Integer) loc.getFieldValue("start")).intValue();
                int end = ((Integer) loc.getFieldValue("end")).intValue();
                // loop through QTLSpans, adding QTLs to a set when gene is spanned by QTL
                Set<QTL> qtlCollection = new HashSet<QTL>();
                for (QTLSpan qtlSpan : qtlSpanSet) {
                    if (chrId.equals(qtlSpan.chromosomeId) && start<=qtlSpan.end && end>=qtlSpan.start) {
                        String qtlId = (String) qtlSpan.qtl.getFieldValue("primaryIdentifier");
                        LOG.info("gene "+geneId+" is spanned by QTL "+qtlId);
                        qtlCollection.add(qtlSpan.qtl);
                    }
                }
                // now add the collection to the gene if it contains QTLs
                if (qtlCollection.size()>0) {
                    Gene tempGene = PostProcessUtil.cloneInterMineObject(gene);
                    tempGene.setFieldValue("spanningQTLs", qtlCollection);
                    osw.store(tempGene);
                }
            } catch (IllegalAccessException ex) {
                throw new RuntimeException(ex);
            }
        }
        
        // close transaction
        osw.commitTransaction();


        // ------------------------------------------------------------------------------------------------------------------------------
        // 2 - Spin through the polypeptides, mRNAs and genes, relating the polypeptides to the genes.
        // ------------------------------------------------------------------------------------------------------------------------------

        LOG.info("Spin through polypeptides, associate with the gene associated with the mRNA...");

        Query qPolypeptide = new Query();
        qPolypeptide.setDistinct(false);
        ConstraintSet csPolypeptide = new ConstraintSet(ConstraintOp.AND);

        // 0 Polypeptide
        QueryClass qcPolypeptide = new QueryClass(Polypeptide.class);
        qPolypeptide.addFrom(qcPolypeptide);
        qPolypeptide.addToSelect(qcPolypeptide);
        qPolypeptide.addToOrderBy(qcPolypeptide);

        // 1 Polypeptide.MRNA
        QueryClass qcPolypeptideMRNA = new QueryClass(MRNA.class);
        qPolypeptide.addFrom(qcPolypeptideMRNA);
        qPolypeptide.addToSelect(qcPolypeptideMRNA);
        qPolypeptide.addToOrderBy(qcPolypeptideMRNA);
        QueryObjectReference polypeptideMRNA = new QueryObjectReference(qcPolypeptide, "MRNA");
        csPolypeptide.addConstraint(new ContainsConstraint(polypeptideMRNA, ConstraintOp.CONTAINS, qcPolypeptideMRNA));

        // 2 MRNA.gene
        QueryClass qcMRNAGene = new QueryClass(Gene.class);
        qPolypeptide.addFrom(qcMRNAGene);
        qPolypeptide.addToSelect(qcMRNAGene);
        qPolypeptide.addToOrderBy(qcMRNAGene);
        QueryObjectReference mRNAGene = new QueryObjectReference(qcPolypeptideMRNA, "gene");
        csPolypeptide.addConstraint(new ContainsConstraint(mRNAGene, ConstraintOp.CONTAINS, qcMRNAGene));
        
        // set the overall constraint
        qPolypeptide.setConstraint(csPolypeptide);

        // execute the query
        Results polypeptideResults = osw.getObjectStore().execute(qPolypeptide);
        Iterator<?> polypeptideIter = polypeptideResults.iterator();

        // begin transaction
        osw.beginTransaction();
        
        while (polypeptideIter.hasNext()) {
            try {
                ResultsRow<?> rr = (ResultsRow<?>) polypeptideIter.next();
                // objects
                Polypeptide polypeptide = (Polypeptide) rr.get(0);
                MRNA mRNA = (MRNA) rr.get(1);
                Gene gene = (Gene) rr.get(2);
                // add the gene to the polypeptide and store it
                Polypeptide tempPoly = PostProcessUtil.cloneInterMineObject(polypeptide);
                tempPoly.setFieldValue("gene", gene);
                osw.store(tempPoly);
            } catch (IllegalAccessException ex) {
                throw new RuntimeException(ex);
            }
        }
        
        // close transaction
        osw.commitTransaction();

    }

    /**
     * Encapsulate a QTL and the genomic span: chromosome, start and end. Only need chromosome ID.
     */
    class QTLSpan {

        QTL qtl;
        String chromosomeId;
        int start;
        int end;

        QTLSpan(QTL qtl, String chromosomeId, int start, int end) {
            this.qtl = qtl;
            this.chromosomeId = chromosomeId;
            this.start = start;
            this.end = end;
        }

    }
        
}
