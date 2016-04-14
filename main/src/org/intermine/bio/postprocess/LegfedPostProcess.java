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
import org.intermine.model.bio.Gene;
import org.intermine.model.bio.GeneticMarker;
import org.intermine.model.bio.Location;
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
 * Find "spanned genes" for each QTL+markers.
 *
 * Since the QTLs and genetic markers can be loaded from a variety of sources (flat files, chado), it makes sense to do this in post-processing when the 
 * QTLs, markers and genes exist in the database in a standard format. 
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
        
        // ------------------------------------------------------------------------------------------
        // First section: accumulate the QTLs and their genomic spans, from their associated markers
        // ------------------------------------------------------------------------------------------

        Query qQTL = new Query();
        qQTL.setDistinct(false);
        ConstraintSet csQTL = new ConstraintSet(ConstraintOp.AND);

        // 0
        QueryClass qcQTL = new QueryClass(QTL.class);
        qQTL.addFrom(qcQTL);
        qQTL.addToSelect(qcQTL);
        qQTL.addToOrderBy(qcQTL);

        // 1
        QueryClass qcGeneticMarker = new QueryClass(GeneticMarker.class);
        qQTL.addFrom(qcGeneticMarker);
        qQTL.addToSelect(qcGeneticMarker);
        qQTL.addToOrderBy(qcGeneticMarker);

        QueryCollectionReference qtlGeneticMarkers = new QueryCollectionReference(qcQTL, "associatedGeneticMarkers");
        csQTL.addConstraint(new ContainsConstraint(qtlGeneticMarkers, ConstraintOp.CONTAINS, qcGeneticMarker));

        // 2
        QueryClass qcChromosome = new QueryClass(Chromosome.class);
        qQTL.addFrom(qcChromosome);
        qQTL.addToSelect(qcChromosome);
        qQTL.addToOrderBy(qcChromosome);

        QueryObjectReference gmChromosome = new QueryObjectReference(qcGeneticMarker, "chromosome");
        csQTL.addConstraint(new ContainsConstraint(gmChromosome, ConstraintOp.CONTAINS, qcChromosome));

        // 3
        QueryClass qcLocation = new QueryClass(Location.class);
        qQTL.addFrom(qcLocation);
        qQTL.addToSelect(qcLocation);
        qQTL.addToOrderBy(qcLocation);

        QueryObjectReference gmLocation = new QueryObjectReference(qcGeneticMarker, "chromosomeLocation");
        csQTL.addConstraint(new ContainsConstraint(gmLocation, ConstraintOp.CONTAINS, qcLocation));

        qQTL.setConstraint(csQTL);
        ObjectStore osQTL = osw.getObjectStore();

        ((ObjectStoreInterMineImpl) osQTL).precompute(qQTL, Constants.PRECOMPUTE_CATEGORY);
        Results qtlResults = osQTL.execute(qQTL, 500, true, true, true);

        osw.beginTransaction();

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
        
        Iterator<?> qtlIter = qtlResults.iterator();
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
                    if (maxEnd>0) {
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
            } catch (Exception ex) {
                LOG.error(ex.toString());
            }
        }
        // last one
        qtlSpanSet.add(new QTLSpan(lastQTL, lastChrId, minStart, maxEnd));

        // close transaction and start another one
        osw.batchCommitTransaction();

        // ------------------------------------------------------------------------------------------------------------------------------
        // Second section: spin through the genes, comparing their genomic range to the QTL genomic spans, associate them if overlapping
        // ------------------------------------------------------------------------------------------------------------------------------

        Query qGene = new Query();
        qGene.setDistinct(false);
        ConstraintSet csGene = new ConstraintSet(ConstraintOp.AND);

        // 0
        QueryClass qcGene = new QueryClass(Gene.class);
        qGene.addFrom(qcGene);
        qGene.addToSelect(qcGene);
        qGene.addToOrderBy(qcGene);

        // 1
        QueryClass qcGeneChromosome = new QueryClass(Chromosome.class);
        qGene.addFrom(qcGeneChromosome);
        qGene.addToSelect(qcGeneChromosome);
        qGene.addToOrderBy(qcGeneChromosome);

        QueryObjectReference geneChromosome = new QueryObjectReference(qcGene, "chromosome");
        csGene.addConstraint(new ContainsConstraint(geneChromosome, ConstraintOp.CONTAINS, qcGeneChromosome));

        // 2
        QueryClass qcGeneLocation = new QueryClass(Location.class);
        qGene.addFrom(qcGeneLocation);
        qGene.addToSelect(qcGeneLocation);
        qGene.addToOrderBy(qcGeneLocation);

        QueryObjectReference geneLocation = new QueryObjectReference(qcGene, "chromosomeLocation");
        csGene.addConstraint(new ContainsConstraint(geneLocation, ConstraintOp.CONTAINS, qcGeneLocation));

        qGene.setConstraint(csGene);
        ObjectStore osGene = osw.getObjectStore();

        ((ObjectStoreInterMineImpl) osGene).precompute(qGene, Constants.PRECOMPUTE_CATEGORY);
        Results geneResults = osGene.execute(qGene, 500, true, true, true);

        Iterator<?> geneIter = geneResults.iterator();
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
                    Gene tempGene;
                    tempGene = PostProcessUtil.cloneInterMineObject(gene);
                    tempGene.setFieldValue("spanningQTLs", qtlCollection);
                    osw.store(tempGene);
                }
            } catch (Exception ex) {
                LOG.error(ex.toString());
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
