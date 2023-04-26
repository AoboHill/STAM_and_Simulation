package STAM;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.util.Log;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.Nucleotide;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

@Description("Pre-process the alignment data to the allele frequency data (four-dimensional)")

public class STAMData extends Alignment{
    public Input<beast.evolution.alignment.Alignment> m_rawData = new Input<beast.evolution.alignment.Alignment>("rawdata","raw binary sequences");
    public Input<List<TaxonSet>> m_taxonsets = new Input<List<TaxonSet>>("taxonset","set of taxons that group a number of SNP sequences into a single sequence",
            new ArrayList<TaxonSet>());

    public List<List<Integer>> nrOfLineages;
    public int [][]  m_nPatternLineageCounts;





    public Map<String, int[][]> countsData; // key: taxon, int[eachSite][counts for AGCT]

    List<String> taxonNames;

    int[][][] counts;


    // The main data used for calculation
    // [pattern][eachsite][counts for AGCT]
    int[][][] countsPatterns;




    public STAMData() {
        sequenceInput.setRule(Validate.OPTIONAL);
    }



    public void initAndValidate() {

        // guess taxon set if no sequences and no taxonsets are known
        if (/*m_pSequences.get().size() == 0 && */m_taxonsets.get().size() == 0 && m_rawData.get() != null) {
            while (sequenceInput.get().size() > 0) {
                sequenceInput.get().remove(0);
            }
            // by last separator
            int nIgnored = guessTaxonSets("^(.+)[-_\\. ](.*)$", 1);
            if (nIgnored > 0) {
                // by first separator
                nIgnored = guessTaxonSets("^([^-_\\. ]+)[-_\\. ](.*)$", 1);
            }
            if (nIgnored > 0) {
                // by taxon name
                nIgnored = guessTaxonSets("^(.*)$", 0);
            }
        }

        // amalgamate binary sequences into count sequences by taxon sets
        if (m_taxonsets.get().size() > 0) {

            // Get the raw data
            List<Sequence> sequences = m_rawData.get().sequenceInput.get();

            DataType rawDataType = m_rawData.get().getDataType();

            if (rawDataType instanceof Nucleotide) {

                int numOfTaxon = m_rawData.get().getTaxonCount();

                int alignmentLength = sequences.get(0).dataInput.get().length();

                Map<String, int[][]> myData = new HashMap<>();
                //new int[numOfTaxon][alignmentLength][4];
                //initialize the data

                for(TaxonSet set : m_taxonsets.get()) {
                    myData.put(set.getID(),new int[alignmentLength][4]);
                    for (Taxon taxon : set.taxonsetInput.get()) {
                        boolean bFound = false;
                        for (int i = 0; i < sequences.size() && !bFound; i++) {
                            if (sequences.get(i).taxonInput.get().equals(taxon.getID())) {
                                char[] seqChar = sequences.get(i).dataInput.get().toCharArray();
                                myData.put(set.getID(),countThisSeq(myData, set.getID(), seqChar));


                                bFound = true;
                            }

                        }
                        if (!bFound) {
                            String seq = m_rawData.get().defaultInput.get().get(taxon.getID());
                            if (seq != null) {
                                char[] seqChar = seq.toCharArray();
                                myData.put(set.getID(),countThisSeq(myData, set.getID(), seqChar));
                            } else {
                                throw new IllegalArgumentException("Could not find taxon " + taxon.getID() + " in alignment");
                            }
                        }
                    }

                }
                this.countsData = myData;
            } else {
                throw new RuntimeException("Cannot handle data of type " + rawDataType.getTypeDescription() +
                        ". Use binary or nucleotide data instead.");
            }


        }


        String oldSiteWeights = siteWeightsInput.get();
        if (oldSiteWeights == null && m_rawData.get() != null) {
            siteWeightsInput.setValue(m_rawData.get().siteWeightsInput.get(), this);
        }

        //System.out.println("Sequence input");
        //System.out.println(sequenceInput.get());


        sequenceInput.setValue(m_rawData.get().sequenceInput.get(),this);
        findDataTypes();
        super.initAndValidate();


        // Set the taxa names for the xml testing
        taxaNames.clear();
        for (TaxonSet i : m_taxonsets.get()){
            taxaNames.add(i.getID());
        }

        if (true){
            calcPatterns();
            if (m_rawData.get() != null) {
                sequenceInput.get().clear();
                siteWeightsInput.setValue(oldSiteWeights, this);
            }

        }

        Log.info.println(toString(false));

    } // initAndValidate


    // This is for counting the number of A,G,C,T in sequences
    public int[][] countThisSeq(Map<String, int[][]> myData,String taxonSetID, char[] seqChar){
        int index = 0;
        int[][] waitForUpdate = myData.get(taxonSetID);
        for (char i : seqChar){
            if (i == 'A'){
                waitForUpdate[index][0]++;
            } else if (i == 'G') {
                waitForUpdate[index][1]++;
            } else if(i == 'C'){
                waitForUpdate[index][2]++;
            } else{
                waitForUpdate[index][3]++;
            }
            index++;
        }
        return waitForUpdate;
    }

    public int guessTaxonSets(String sRegexp, int nMinSize) {
        m_taxonsets.get().clear();
        List<Taxon> taxa = new ArrayList<Taxon>();
        for (String name : m_rawData.get().getTaxaNames()) { //Sequence sequence : m_rawData.get().sequenceInput.get()) {
            Taxon taxon = new Taxon();
            // ensure sequence and taxon do not get same ID
            //if (sequence.getID() == null || sequence.getID().equals(sequence.taxonInput.get())) {
            //	sequence.setID("_"+sequence.getID());
            //}
            taxon.setID(name); //sequence.taxonInput.get());
            taxa.add(taxon);
        }
        HashMap<String, TaxonSet> map = new HashMap<String, TaxonSet>();
        Pattern m_pattern = Pattern.compile(sRegexp);
        int nIgnored = 0;
        for (Taxon taxon : taxa) {
            Matcher matcher = m_pattern.matcher(taxon.getID());
            if (matcher.find()) {
                String sMatch = matcher.group(1);
                try {
                    if (map.containsKey(sMatch)) {

                        TaxonSet set = map.get(sMatch);
                        set.taxonsetInput.setValue(taxon, set);

                    } else {
                        TaxonSet set = new TaxonSet();
                        if (sMatch.equals(taxon.getID())) {
                            set.setID(sMatch + "_");
                        } else {
                            set.setID(sMatch);
                        }
                        set.taxonsetInput.setValue(taxon, set);
                        map.put(sMatch, set);
                    }
                } catch (Exception ex) {
                    ex.printStackTrace();
                }
            } else {
                nIgnored++;
            }
        }
        // add taxon sets
        for (TaxonSet set : map.values()) {
            if (set.taxonsetInput.get().size() > nMinSize) {
                m_taxonsets.setValue(set, this);
            } else {
                nIgnored += set.taxonsetInput.get().size();
            }
        }
        return nIgnored;
    }

    // Get the weights of the patterns
    public int getPatternWeight(int id) {
        if (id < patternWeight.length) {
            return patternWeight[id];
        }
        return 0;
    }


    // Get the counts pattern data, used in likelihood calculation
    public int[][][] getCountsPatterns(){
        return countsPatterns;
    }

    public int getCountsPatternsCounts(){
        return countsPatterns.length;
    }


    /** calculate patterns from sequence data
     * The difference with standard sequence data is that constant sites
     * are removed + 2 patterns are added at the end representing these
     * constant sites, but with zero weight. The likelihood calculator
     * deals with these different sites.
     *
     * **/
    // for the counts data, we design a new method to calculate the unique pattern
    @Override
    public void calcPatterns() {
        // find unique patterns
        Set<String> taxa = countsData.keySet();

        int nTaxa = taxa.size();
        // find unique patterns and find the weights at the meantime
        int nSeq = m_rawData.get().sequenceInput.get().get(0).dataInput.get().length();
        int[] weights = new int[nSeq];
        Boolean[] isVisited = new Boolean[nSeq];

        for (int i = 0; i < weights.length; i++) {
            int j = 0;
            for (j = 0; j < i; j++) {
                if (isEqual(i,j)) {
                    break;
                }
            }
            weights[j]++;
        }
        //System.out.println(Arrays.toString(weights));

        // count nr of patterns
        int nPatterns = 0;
        for (int i = 0; i < weights.length; i++) {
            if (weights[i]>0) {
                nPatterns++;
            }
        }
        //System.out.println(nPatterns);

        patternWeight = new int[nPatterns];
        countsPatterns = new int[nPatterns][nTaxa][4];
        patternIndex = new int[nSeq];

        nPatterns = 0;
        int iSite = 0;


        // instantiate patterns
        for (int i = 0; i < nSeq; i++) {
            if (weights[i]>0) {
                patternWeight[nPatterns] = weights[i];
                List<String> taxonList = new ArrayList<>(countsData.keySet());
                for (int j = 0; j < nTaxa; j++) {
                    countsPatterns[nPatterns][j] = countsData.get(taxonList.get(j))[i];

                }
                for (int k = 0; k < weights[i]; k++) {
                    patternIndex[iSite++] = nPatterns;
                }
                nPatterns++;
            }
        }




        // pattern
        List<String> taxonList = new ArrayList<>(countsData.keySet());
        Arrays.fill(patternIndex, -1);
        for (int i = 0; i < nSeq; i++) {
            int[][] siteCounts = new int[nTaxa][4];
            for (int j = 0; j < nTaxa; j++) {
                siteCounts[j] = countsData.get(taxonList.get(j))[i];
            }
            for (int j = 0; j < nPatterns; j++) {
                boolean found = true;
                for (int k = 0; k < nTaxa; k++) {
                    if (!Arrays.equals(siteCounts[k], countsPatterns[j][k])) {
                        found = false;
                        break;
                    }
                }
                if (found) {
                    patternIndex[i] = j;
                    j = nPatterns;
                }
            }
        }



    }


    @Override
    public int getPatternCount() {
        return patternWeight.length;
    }


    // check if it is the unique pattern
    protected boolean isEqual ( int iSite1, int iSite2){
        for (String taxon: countsData.keySet()){
            //System.out.println(Arrays.toString(countsData.get(taxon)[iSite1]));
            //System.out.println(Arrays.toString(countsData.get(taxon)[iSite2]));
            if (!Arrays.equals(countsData.get(taxon)[iSite1], countsData.get(taxon)[iSite2])){
                return false;
            }
        }
        return true;
    }



    private long getTotalWeight() {
        long totalWeight = 0;
        for (int weight : patternWeight) {
            totalWeight += weight;
        }
        return totalWeight;
    }

    @Override
    public String toString(boolean singleLine) {
        long totalWeight = getTotalWeight();
        StringBuilder builder = new StringBuilder();
        builder.append(getClass().getSimpleName() + "(" + getID() + ")");

        if (singleLine) {
            builder.append(": [taxa, patterns, sites] = [" + getTaxonCount() + ", " + getPatternCount());
            builder.append(", " + getTotalWeight() + "]");
        } else {

            long siteCount = getSiteCount();

            builder.append('\n');
            builder.append("  " + getTaxonCount() + " taxa");
            builder.append('\n');
            builder.append("  " + siteCount + (siteCount == 1 ? " site" : " sites") + (totalWeight == getSiteCount() ? "" : " with weight " + totalWeight + ""));
            builder.append('\n');
            if (siteCount > 1) {
                builder.append("  " + getPatternCount() + " patterns");
                builder.append('\n');
            }
        }
        return builder.toString();
    }

}
