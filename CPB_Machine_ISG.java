
package cpb_machine_isg;

import java.io.*;
import java.util.HashMap;
import java.util.Random;

public class CPB_Machine_ISG {

    public static String[] expOrder;
    public static String[][] aaCode, aminoAcids, aminoCodons;
    public static int[] aminoCodonCounts, codCounts, aaCounts;
    public static int[][] codPairCounts, aaPairCounts;
    
    public static HashMap<String,Integer> aaLookUp, codLookUp;

    public final static int aminoNum=21;
    public static int count=0, a=0, c=0, g=0, t=0, n=0, cg=0, ua=0, totCod=0, nullCod=0, totCP=0, stops=0, totLength=0, totSeqs=0, cpbCount;
    public static double gc=0, at=0, cpg=0, upa=0;
    
    public static int aaN=0, acN=0, agN=0, atN=0;
    public static double apa=0, apc=0, apg=0, apt=0;
    public static int caN=0, ccN=0, cgN=0, ctN=0;
    public static double cpa=0, cpc=0, cpt=0;//cpg=0
    public static int gaN=0, gcN=0, ggN=0, gtN=0;
    public static double gpa=0, gpc=0, gpg=0, gpt=0;
    public static int taN=0, tcN=0, tgN=0, ttN=0;
    public static double tpa=0, tpc=0, tpg=0, tpt=0;//upa=0
    
    public static int totBr=0, totNonBr=0, totBrN=0, totNonBrN=0;
    public static int aBr=0, cBr=0, gBr=0, tBr=0, nBr=0,aNonBr=0, cNonBr=0, gNonBr=0, tNonBr=0, nNonBr=0;
    public static int aaBrN=0, acBrN=0, agBrN=0, atBrN=0, caBrN=0, ccBrN=0, cgBrN=0, ctBrN=0, gaBrN=0, gcBrN=0, ggBrN=0, gtBrN=0, taBrN=0, tcBrN=0, tgBrN=0, ttBrN=0;
    public static int aaNonBrN=0, acNonBrN=0, agNonBrN=0, atNonBrN=0, caNonBrN=0, ccNonBrN=0, cgNonBrN=0, ctNonBrN=0, gaNonBrN=0, gcNonBrN=0, ggNonBrN=0, gtNonBrN=0, taNonBrN=0, tcNonBrN=0, tgNonBrN=0, ttNonBrN=0;
    public static double apaBr=0, apcBr=0, apgBr=0, aptBr=0,cpaBr=0, cpcBr=0, cpgBr=0, cptBr=0, gpaBr=0, gpcBr=0, gpgBr=0, gptBr=0, tpaBr=0, tpcBr=0, tpgBr=0, tptBr=0;
    public static double apaNonBr=0, apcNonBr=0, apgNonBr=0, aptNonBr=0,cpaNonBr=0, cpcNonBr=0, cpgNonBr=0, cptNonBr=0, gpaNonBr=0, gpcNonBr=0, gpgNonBr=0, gptNonBr=0, tpaNonBr=0, tpcNonBr=0, tpgNonBr=0, tptNonBr=0;
    
    public static double codBias[], aaBias[], cpBias[][], cpb[][];
    public static double cpbAv=0, cpbSum=0, trueCpbMin, cpbMin;
    public static String filename="", selFilename="", outname="", seq="", seqName="", prevName;
    
    public static int seqCount=0;
    
    public static boolean average=true, outType=true;
    
    public static int badCount=0, bad=0, notGood=0, notCount=0;
  
    public static String nameHeader="", selData[][];
    
    public static int polCol=0;
    
    public static double rawAA=0, rawAC=0, rawAG=0, rawAT=0, rawCA=0, rawCC=0, rawCG=0, rawCT=0, rawGA=0, rawGC=0, rawGG=0, rawGT=0, rawTA=0, rawTC=0, rawTG=0, rawTT=0;
    
    public static void main(String[] args) {

        if(args.length==1) {
            //FASTA seqs - single line format essential
            filename=args[0];
        }
        else if(args.length==2) {
            //FASTA seqs - single line format essential
            filename=args[0];
            selFilename=args[1];
        }
        else {
            System.out.println("Error - incorrect parameters "+args.length);
            //System.exit(0);
            filename="/Users/rjorton/Desktop/RichardRefSeq/seqs_new_acc.txt";
            filename="/Users/rjorton/Desktop/RichardRefSeq/bunya_extra_sort.fasta";
            filename="/Users/rjorton/Desktop/RichardRefSeq/Bridge/toga_new.fasta";
            filename="/Users/rjorton/Desktop/RichardRefSeq/arena_new_sort.fasta";
            filename="/Users/rjorton/Desktop/RichardRefSeq/May2017/bunya_new_acc_sort_new.fasta";
            filename="/Users/rjorton/Dropbox/CPBmachine/Gifford/raw sequence data/V2/Parvovirinae_new.fasta";
            filename="/Users/rjorton/Dropbox/CPBcamels/astro_new.fasta";
            filename="/Users/rjorton/Desktop/RichardRefSeq/Bridge/crap/corona_new.fasta";
            filename="/Users/rjorton/Downloads/ebola_all_coding_new.fasta";
            
            //To join to other data file containing meta-data - match accessions
            selFilename="/Users/rjorton/Desktop/RichardRefSeq/May2017/virus_accessionsV3B.txt";
        }

        outname=filename.substring(0, filename.lastIndexOf("."))+"_cpb_dat.txt";
        
        System.out.println("In Seq File= "+filename);
        System.out.println("Out Dat File = "+outname);
        
        
        File inFile = new File(selFilename);
        int count=0;
        try {
            BufferedReader input =  new BufferedReader(new FileReader(inFile));

            try {
                String line = null;

                while (( line = input.readLine()) != null) {
                    count++;
                }
              }
            finally {
                input.close();
            }
        }
        catch (IOException ex) {
            ex.printStackTrace();
        }
        selData=new String[count][4];
        count=0;
        try {
            BufferedReader input =  new BufferedReader(new FileReader(inFile));

            try {
                String line = null;

                while (( line = input.readLine()) != null) {
                    String splits[]=line.split("\t");
                    selData[count][0]=splits[0];//accession
                    selData[count][1]=splits[1];//taxID
                    selData[count][2]=splits[2];//virus name
                    selData[count][3]=splits[3];//family
                    count++;
                }
              }
            finally {
                input.close();
            }
        }
        catch (IOException ex) {
            ex.printStackTrace();
        }
        count=0;
        
        createCode();
        
        codCounts=new int[64];
        codBias=new double[64];
        
        aaCounts=new int[21];
        aaBias=new double[21];
        
        codPairCounts=new int[64][64];
        cpBias=new double[64][64];
        cpb=new double[64][64];
        
        aaPairCounts=new int[21][21];
        
        try{
            FileWriter fstream = new FileWriter(outname);
            BufferedWriter out = new BufferedWriter(fstream);
        
            out.write(genHeader());
            
            try {
                File seqFile = new File(filename);
                BufferedReader input =  new BufferedReader(new FileReader(seqFile));

                try {
                    String line = null;
                    String prevLine=line;
                    //String acc="";
                    seqName="START";
                    count=0;

                    while (( line = input.readLine()) != null) {
                        
                        if(line.indexOf(">")==0) {
                            int from=5;
                            int to=line.indexOf("_cds_");
 
                       
                            from=0;
                            to=line.length();
                            
                            prevName=seqName;
                            seqName=line.substring(from,to);
                            
                            if(count==0)
                                prevName=seqName;
                            
                            seqCount++;
                        }
                        else {   
                            //if the current sequence name/accesion does not equal the previous
                            //means new sequence - output the existing the data and then clear it
                            if(!seqName.equals(prevName)) {
                                out.write(genOutput());
                                clearDat();
                            }
                            
                            seq=line.toUpperCase();
                            
                            int trunc=(int)((double)seq.length()/(double)3);
                            trunc=trunc*3;
                            seq=seq.substring(0, trunc);
                                
                            if(seq.length()%3!=0)
                                System.out.println("Error - seq length not divisible by 3 "+prevName+" "+seq.length()+" "+prevLine);
                            //Not really needed anymore
                            else if(seq.length()<3)
                                System.out.println("Error - seq length less than 3! "+seq.length());
                            else
                                analyseSeq(seq);
                        }
                        
                        prevLine=line;   
                        count++;
                    }
                }
                finally {
                    input.close();
                }

            }
            catch (IOException ex) {
                ex.printStackTrace();
            }
        
            //catch last
            prevName=seqName;
            out.write(genOutput());
            out.close();
        }

        catch (Exception e) {
            e.printStackTrace();
            System.err.println("Error: " + e.getMessage());
        }
        

    }
    
    public static void clearDat() {
        a=c=g=t=n=cg=ua=totCod=nullCod=totCP=stops=0;
        aaN=acN=agN=atN=caN=ccN=cgN=ctN=gaN=gcN=ggN=gtN=taN=tcN=tgN=ttN=0;
        apa=apc=apg=apt=cpa=cpc=cpg=cpt=gpa=gpc=gpg=gpt=tpa=tpc=tpg=tpt=0;
        
        aBr=cBr=gBr=tBr=nBr=aNonBr=cNonBr=gNonBr=tNonBr=nNonBr=0;
        aaBrN=acBrN=agBrN=atBrN=caBrN=ccBrN=cgBrN=ctBrN=gaBrN=gcBrN=ggBrN=gtBrN=taBrN=tcBrN=tgBrN=ttBrN=0;
        aaNonBrN=acNonBrN=agNonBrN=atNonBrN=caNonBrN=ccNonBrN=cgNonBrN=ctNonBrN=gaNonBrN=gcNonBrN=ggNonBrN=gtNonBrN=taNonBrN=tcNonBrN=tgNonBrN=ttNonBrN=0;
        apaBr=apcBr=apgBr=aptBr=cpaBr=cpcBr=cpgBr=cptBr=gpaBr=gpcBr=gpgBr=gptBr=tpaBr=tpcBr=tpgBr=tptBr=0;
        apaNonBr=apcNonBr=apgNonBr=aptNonBr=cpaNonBr=cpcNonBr=cpgNonBr=cptNonBr=gpaNonBr=gpcNonBr=gpgNonBr=gptNonBr=tpaNonBr=tpcNonBr=tpgNonBr=tptNonBr=0;
        totBr=totNonBr=totBrN=totNonBrN=0;
        
        gc=at=0;
        cpbSum=cpbAv=0;
        cpbCount=0;
        totLength=totSeqs=0;
        cpg=upa=0;
        cpbMin=1000000;
        trueCpbMin=1000000;
        
        for(int i=0;i<codCounts.length;i++)
            codCounts[i]=0;
        
        for(int i=0;i<codBias.length;i++)
            codBias[i]=0;
        
        for(int i=0;i<aaCounts.length;i++)
            aaCounts[i]=0;
        
        for(int i=0;i<aaBias.length;i++)
            aaBias[i]=0;
        
        for(int i=0;i<aaPairCounts.length;i++)
            for(int j=0;j<aaPairCounts.length;j++)
                aaPairCounts[i][j]=0;
        
        for(int i=0;i<codPairCounts.length;i++)
            for(int j=0;j<codPairCounts[i].length;j++)
                codPairCounts[i][j]=0;
        
        for(int i=0;i<cpBias.length;i++)
            for(int j=0;j<cpBias.length;j++)
                cpBias[i][j]=0;
        
        for(int i=0;i<cpb.length;i++)
            for(int j=0;j<cpb.length;j++)
                cpb[i][j]=0;
        
        
        rawAA=rawAC=rawAG=rawAT=0;
        rawCA=rawCC=rawCG=rawCT=0;
        rawGA=rawGC=rawGG=rawGT=0;
        rawTA=rawTC=rawTG=rawTT=0;
    }
    
    public static String genHeader() {
        
        String head="TaxID\tSpecies\tSeqName\tGood\tComplete\tSeqs\tSeqLength\tCodons\tBadCodons\tCodonPairs\tStops\t";
        head+="A\tC\tG\tT\tN\t";
        head+="GC\tAT\tCpG\tUpA\t";
//
//        for(int j=0;j<aaCounts.length;j++)
//            head+=aminoAcids[0][j]+"-Count\t";

        for(int j=0;j<aaBias.length;j++)
            head+=aminoAcids[0][j]+"-Bias\t";
//  
//        for(int j=0;j<codCounts.length;j++)
//            head+=aaCode[0][j]+"-Counts\t";

        for(int j=0;j<codBias.length;j++)
            head+=aaCode[0][j]+"-Bias\t";
//
//        for(int j=0;j<aaPairCounts.length;j++)
//            for(int k=0;k<aaPairCounts[j].length;k++)
//                    head+=aminoAcids[0][j]+"-"+aminoAcids[0][k]+"\t";
        

        for(int j=0;j<codPairCounts.length;j++)
            for(int k=0;k<codPairCounts[j].length;k++)
                    head+=aaCode[0][j]+"("+aaCode[1][j]+")-"+aaCode[0][k]+"("+aaCode[1][k]+")\t";
   
        head+="CpsMin*2\tCpsAv";
        //standard - without cpg and upa
        head+="\tApA\tApC\tApG\tApU\tCpA\tCpC\tCpU\tGpA\tGpC\tGpG\tGpU\tUpC\tUpG\tUpT";
        //bridge
        head+="\tbrApA\tbrApC\tbrApG\tbrApU\tbrCpA\tbrCpC\tbrCpG\tbrCpU\tbrGpA\tbrGpC\tbrGpG\tbrGpU\tbrUpA\tbrUpC\tbrUpG\tbrUpU";
        //non-bridge
        head+="\tNonBrApA\tNonBrApC\tNonBrApG\tNonBrApU\tNonBrCpA\tNonBrCpC\tNonBrCpG\tNonBrCpU\tNonBrGpA\tNonBrGpC\tNonBrGpG\tNonBrGpU\tNonBrUpA\tNonBrUpC\tNonBrUpG\tNonBrUpU";
        
        head+="\tRawAA\tRawAC\tRawAG\tRawAT";
        head+="\tRawCA\tRawCC\tRawCG\tRawCT";
        head+="\tRawGA\tRawGC\tRawGG\tRawGT";
        head+="\tRawTA\tRawTC\tRawTG\tRawTT";
        
        head+="\n";
        
        return head;
    }
    
    public static String genOutput() {
        
        String dat="";
        boolean test2=false;
        
        for(int i=0;i<selData.length;i++) {
            if(selData[i][0].equals(prevName)) {
               dat+=selData[i][1]+"\t"+selData[i][2]+"\t"; 
               test2=true;
               break;
            }
        }
        
        if(!test2) {
            dat+="?\t?\t";
        }
        
        dat+=prevName+"\t"+bad+"\t"+notGood+"\t"+totSeqs+"\t"+totLength+"\t"+totCod+"\t"+nullCod+"\t"+totCP+"\t"+" "+stops+"\t";
        dat+=((double)a/(double)totLength)+"\t"+((double)c/(double)totLength)+"\t"+((double)g/(double)totLength)+"\t"+((double)t/(double)totLength)+"\t"+((double)n/(double)totLength)+"\t";
        dat+=gc+"\t"+at+"\t"+cpg+"\t"+upa+"\t";

//        for(int j=0;j<aaCounts.length;j++)
//            dat+=aaCounts[j]+"\t";

        for(int j=0;j<aaBias.length;j++)
            dat+=aaBias[j]+"\t";

//        for(int j=0;j<codCounts.length;j++)
//            dat+=codCounts[j]+"\t";

        for(int j=0;j<codBias.length;j++)
            dat+=codBias[j]+"\t";
    
//        for(int j=0;j<aaPairCounts.length;j++)
//            for(int k=0;k<aaPairCounts[j].length;k++)
//                    dat+=aaPairCounts[j][k]+"\t";
   
//        for(int j=0;j<codPairCounts.length;j++)
//            for(int k=0;k<codPairCounts[j].length;k++)
//                    dat+=codPairCounts[j][k]+"\t";
        
        if(outType) {
            for(int j=0;j<cpb.length;j++)
                for(int k=0;k<cpb[j].length;k++)
                        dat+=cpb[j][k]+"\t";
        }
        else {
            for(int j=0;j<cpBias.length;j++)
                for(int k=0;k<cpBias[j].length;k++)
                        dat+=cpBias[j][k]+"\t";
        }
        
        dat+=cpbMin+"\t"+cpbAv+"\t";
        
        dat+=apa+"\t"+apc+"\t"+apg+"\t"+apt+"\t";
        dat+=cpa+"\t"+cpc+"\t"+cpt+"\t";//removed cpg
        dat+=gpa+"\t"+gpc+"\t"+gpg+"\t"+gpt+"\t";
        dat+=tpc+"\t"+tpg+"\t"+tpt+"\t";//removed tpa
        
        //bridge
        dat+=apaBr+"\t"+apcBr+"\t"+apgBr+"\t"+aptBr+"\t";
        dat+=cpaBr+"\t"+cpcBr+"\t"+cpgBr+"\t"+cptBr+"\t";
        dat+=gpaBr+"\t"+gpcBr+"\t"+gpgBr+"\t"+gptBr+"\t";
        dat+=tpaBr+"\t"+tpcBr+"\t"+tpgBr+"\t"+tptBr+"\t";
        
        //NonBridge
        dat+=apaNonBr+"\t"+apcNonBr+"\t"+apgNonBr+"\t"+aptNonBr+"\t";
        dat+=cpaNonBr+"\t"+cpcNonBr+"\t"+cpgNonBr+"\t"+cptNonBr+"\t";
        dat+=gpaNonBr+"\t"+gpcNonBr+"\t"+gpgNonBr+"\t"+gptNonBr+"\t";
        dat+=tpaNonBr+"\t"+tpcNonBr+"\t"+tpgNonBr+"\t"+tptNonBr+"\t";
        
        //Raw
        dat+=rawAA+"\t"+rawAC+"\t"+rawAG+"\t"+rawAT+"\t";
        dat+=rawCA+"\t"+rawCC+"\t"+rawCG+"\t"+rawCT+"\t";
        dat+=rawGA+"\t"+rawGC+"\t"+rawGG+"\t"+rawGT+"\t";
        dat+=rawTA+"\t"+rawTC+"\t"+rawTG+"\t"+rawTT+"\t";
        
        dat+="\n";
        
        //dat+="CodCount="+codCounts[0]+"\n";
        //dat+="CodPairCount="+codPairCounts[0][0]+"\n";
        //dat+="CodPairCountLookup="+codPairCounts[codLookUp.get("AAA")][codLookUp.get("AAA")]+"\n";
        //dat+="AAcount="+aaCounts[9]+"\n";
        //dat+="AAcountLookup="+aaCounts[aaLookUp.get(aaCode[1][0])]+"\n";
        //dat+="AApairCount="+aaPairCounts[9][9]+"\n";
        //dat+="AApairCountLookup="+aaPairCounts[aaLookUp.get(aaCode[1][codLookUp.get("AAA")])][aaLookUp.get(aaCode[1][codLookUp.get("AAA")])]+"\n";
        //dat+="AApairCountLookup2="+aaPairCounts[aaLookUp.get(aaCode[1][0])][aaLookUp.get(aaCode[1][0])]+"\n";
        //dat+="CPB="+cpb[0][0]+"\n";
        //dat+=Math.log((double)codPairCounts[0][0]/((((double)codCounts[0]*(double)codCounts[0])/(aaCounts[aaLookUp.get(aaCode[1][0])]*aaCounts[aaLookUp.get(aaCode[1][0])]))/aaPairCounts[aaLookUp.get(aaCode[1][0])][aaLookUp.get(aaCode[1][0])]))+"\n";

        //dat+=((double)codPairCounts[0][0]/((((double)codCounts[0]*(double)codCounts[0])/(aaCounts[aaLookUp.get(aaCode[1][0])]*aaCounts[aaLookUp.get(aaCode[1][0])]))/aaPairCounts[aaLookUp.get(aaCode[1][0])][aaLookUp.get(aaCode[1][0])]))+"\n";
                    
            
        
        return dat;

    }
    
    public static void analyseSeq(String seq) {

        //as coding seqs of same genome are merged together keep a running total of all the data
        //and calculate the biases as don't know if this will be the last sequence
        
        totLength+=seq.length();
        totSeqs++;
        
        //Nucleotides
        for(int i=0;i<seq.length();i++) {
            
            //all
            if(seq.charAt(i)=='A')
                a++;
            else if(seq.charAt(i)=='C')
                c++;
            else if(seq.charAt(i)=='G')
                g++;
            else if(seq.charAt(i)=='T')
                t++;
            else
                n++;
            
            //bridge- 3,6,9 etc correspond to codon pos 1 [seq starts at 0]
            if(i%3==0 & i>0) {   
                totBr++;
                
                if(seq.charAt(i-1)=='A')
                    aBr++;
                else if(seq.charAt(i-1)=='C')
                    cBr++;
                else if(seq.charAt(i-1)=='G')
                    gBr++;
                else if(seq.charAt(i-1)=='T')
                    tBr++;
                else
                    nBr++;
                
                totBr++;
                
                if(seq.charAt(i)=='A')
                    aBr++;
                else if(seq.charAt(i)=='C')
                    cBr++;
                else if(seq.charAt(i)=='G')
                    gBr++;
                else if(seq.charAt(i)=='T')
                    tBr++;
                else
                    nBr++;
            }
        }
        
        //Non-bridge actually includes all nucls 123-123 bridge is 3-1, non bridge is 1-2,2-3,1-2,2-3
        aNonBr=a;
        cNonBr=c;
        gNonBr=g;
        tNonBr=t;
        nNonBr=n;
        
        //totLength should be same as a,c,g,t,n
        totNonBr=totLength;
        
        //standard gc/at content - nothing to do with dinucs
        gc=(double)(g+c)/(double)totLength;
        at=(double)(a+t)/(double)totLength;
        
        //Count the dinucleotides - start at 1 as doing current and previous base
        for(int i=1;i<seq.length();i++) {
            
            //originally only did these two cg, ua dinucs so separate
            if(seq.charAt(i-1)=='C' & seq.charAt(i)=='G')
                cg++;
            if(seq.charAt(i-1)=='T' & seq.charAt(i)=='A')
                ua++;
            
            if(seq.charAt(i-1)=='A' & seq.charAt(i)=='A')
                aaN++;
            if(seq.charAt(i-1)=='A' & seq.charAt(i)=='C')
                acN++;
            if(seq.charAt(i-1)=='A' & seq.charAt(i)=='G')
                agN++;
            if(seq.charAt(i-1)=='A' & seq.charAt(i)=='T')
                atN++;
            if(seq.charAt(i-1)=='C' & seq.charAt(i)=='A')
                caN++;
            if(seq.charAt(i-1)=='C' & seq.charAt(i)=='C')
                ccN++;
            if(seq.charAt(i-1)=='C' & seq.charAt(i)=='G')
                cgN++;
            if(seq.charAt(i-1)=='C' & seq.charAt(i)=='T')
                ctN++;
            if(seq.charAt(i-1)=='G' & seq.charAt(i)=='A')
                gaN++;
            if(seq.charAt(i-1)=='G' & seq.charAt(i)=='C')
                gcN++;
            if(seq.charAt(i-1)=='G' & seq.charAt(i)=='G')
                ggN++;
            if(seq.charAt(i-1)=='G' & seq.charAt(i)=='T')
                gtN++;
            if(seq.charAt(i-1)=='T' & seq.charAt(i)=='A')
                taN++;
            if(seq.charAt(i-1)=='T' & seq.charAt(i)=='C')
                tcN++;
            if(seq.charAt(i-1)=='T' & seq.charAt(i)=='G')
                tgN++;
            if(seq.charAt(i-1)=='T' & seq.charAt(i)=='T')
                ttN++;
            
            //BRIDGE-seq start at 0, so 3,6,9 etc are codon pos 1
            if(i%3==0) {
                totBrN++;
                
                if(seq.charAt(i-1)=='A' & seq.charAt(i)=='A')
                    aaBrN++;
                if(seq.charAt(i-1)=='A' & seq.charAt(i)=='C')
                    acBrN++;
                if(seq.charAt(i-1)=='A' & seq.charAt(i)=='G')
                    agBrN++;
                if(seq.charAt(i-1)=='A' & seq.charAt(i)=='T')
                    atBrN++;
                if(seq.charAt(i-1)=='C' & seq.charAt(i)=='A')
                    caBrN++;
                if(seq.charAt(i-1)=='C' & seq.charAt(i)=='C')
                    ccBrN++;
                if(seq.charAt(i-1)=='C' & seq.charAt(i)=='G')
                    cgBrN++;
                if(seq.charAt(i-1)=='C' & seq.charAt(i)=='T')
                    ctBrN++;
                if(seq.charAt(i-1)=='G' & seq.charAt(i)=='A')
                    gaBrN++;
                if(seq.charAt(i-1)=='G' & seq.charAt(i)=='C')
                    gcBrN++;
                if(seq.charAt(i-1)=='G' & seq.charAt(i)=='G')
                    ggBrN++;
                if(seq.charAt(i-1)=='G' & seq.charAt(i)=='T')
                    gtBrN++;
                if(seq.charAt(i-1)=='T' & seq.charAt(i)=='A')
                    taBrN++;
                if(seq.charAt(i-1)=='T' & seq.charAt(i)=='C')
                    tcBrN++;
                if(seq.charAt(i-1)=='T' & seq.charAt(i)=='G')
                    tgBrN++;
                if(seq.charAt(i-1)=='T' & seq.charAt(i)=='T')
                    ttBrN++;
            }
        }
        
        //Nonbridge is total minus bridge
        aaNonBrN=aaN-aaBrN;
        acNonBrN=acN-acBrN;
        agNonBrN=agN-agBrN;
        atNonBrN=atN-atBrN;
        
        caNonBrN=caN-caBrN;
        ccNonBrN=ccN-ccBrN;
        cgNonBrN=cgN-cgBrN;
        ctNonBrN=ctN-ctBrN;
        
        gaNonBrN=gaN-gaBrN;
        gcNonBrN=gcN-gcBrN;
        ggNonBrN=ggN-ggBrN;
        gtNonBrN=gtN-gtBrN;
        
        taNonBrN=taN-taBrN;
        tcNonBrN=tcN-tcBrN;
        tgNonBrN=tgN-tgBrN;
        ttNonBrN=ttN-ttBrN;
        
        //totLength-totSeqs is the total dinucls (length-1 for each seq)
        //tot non bridge dinculs is total dinucs minus bridge dinucs
        totNonBrN=totLength-totSeqs-totBrN;
        
        //calculate dinucle biases
        cpg=((double)cg/((double)totLength-(double)totSeqs))/(((double)c/(double)totLength)*((double)g/(double)totLength));
        upa=((double)ua/((double)totLength-(double)totSeqs))/(((double)t/(double)totLength)*((double)a/(double)totLength));
       
        apa=((double)aaN/((double)totLength-(double)totSeqs))/(((double)a/(double)totLength)*((double)a/(double)totLength));
        apc=((double)acN/((double)totLength-(double)totSeqs))/(((double)a/(double)totLength)*((double)c/(double)totLength));
        apg=((double)agN/((double)totLength-(double)totSeqs))/(((double)a/(double)totLength)*((double)g/(double)totLength));
        apt=((double)atN/((double)totLength-(double)totSeqs))/(((double)a/(double)totLength)*((double)t/(double)totLength));
        
        cpa=((double)caN/((double)totLength-(double)totSeqs))/(((double)c/(double)totLength)*((double)a/(double)totLength));
        cpc=((double)ccN/((double)totLength-(double)totSeqs))/(((double)c/(double)totLength)*((double)c/(double)totLength));
        //cpg=((double)cgN/((double)totLength-(double)totSeqs))/(((double)c/(double)totLength)*((double)g/(double)totLength));
        cpt=((double)ctN/((double)totLength-(double)totSeqs))/(((double)c/(double)totLength)*((double)t/(double)totLength));
        
        gpa=((double)gaN/((double)totLength-(double)totSeqs))/(((double)g/(double)totLength)*((double)a/(double)totLength));
        gpc=((double)gcN/((double)totLength-(double)totSeqs))/(((double)g/(double)totLength)*((double)c/(double)totLength));
        gpg=((double)ggN/((double)totLength-(double)totSeqs))/(((double)g/(double)totLength)*((double)g/(double)totLength));
        gpt=((double)gtN/((double)totLength-(double)totSeqs))/(((double)g/(double)totLength)*((double)t/(double)totLength));
        
        tpa=((double)taN/((double)totLength-(double)totSeqs))/(((double)t/(double)totLength)*((double)a/(double)totLength));
        tpc=((double)tcN/((double)totLength-(double)totSeqs))/(((double)t/(double)totLength)*((double)c/(double)totLength));
        tpg=((double)tgN/((double)totLength-(double)totSeqs))/(((double)t/(double)totLength)*((double)g/(double)totLength));
        tpt=((double)ttN/((double)totLength-(double)totSeqs))/(((double)t/(double)totLength)*((double)t/(double)totLength));
        
        //bridge
        apaBr=((double)aaBrN/((double)totBrN))/(((double)aBr/(double)totBr)*((double)aBr/(double)totBr));
        apcBr=((double)acBrN/((double)totBrN))/(((double)aBr/(double)totBr)*((double)cBr/(double)totBr));
        apgBr=((double)agBrN/((double)totBrN))/(((double)aBr/(double)totBr)*((double)gBr/(double)totBr));
        aptBr=((double)atBrN/((double)totBrN))/(((double)aBr/(double)totBr)*((double)tBr/(double)totBr));
        
        cpaBr=((double)caBrN/((double)totBrN))/(((double)cBr/(double)totBr)*((double)aBr/(double)totBr));
        cpcBr=((double)ccBrN/((double)totBrN))/(((double)cBr/(double)totBr)*((double)cBr/(double)totBr));
        cpgBr=((double)cgBrN/((double)totBrN))/(((double)cBr/(double)totBr)*((double)gBr/(double)totBr));
        cptBr=((double)ctBrN/((double)totBrN))/(((double)cBr/(double)totBr)*((double)tBr/(double)totBr));
        
        gpaBr=((double)gaBrN/((double)totBrN))/(((double)gBr/(double)totBr)*((double)aBr/(double)totBr));
        gpcBr=((double)gcBrN/((double)totBrN))/(((double)gBr/(double)totBr)*((double)cBr/(double)totBr));
        gpgBr=((double)ggBrN/((double)totBrN))/(((double)gBr/(double)totBr)*((double)gBr/(double)totBr));
        gptBr=((double)gtBrN/((double)totBrN))/(((double)gBr/(double)totBr)*((double)tBr/(double)totBr));
        
        tpaBr=((double)taBrN/((double)totBrN))/(((double)tBr/(double)totBr)*((double)aBr/(double)totBr));
        tpcBr=((double)tcBrN/((double)totBrN))/(((double)tBr/(double)totBr)*((double)cBr/(double)totBr));
        tpgBr=((double)tgBrN/((double)totBrN))/(((double)tBr/(double)totBr)*((double)gBr/(double)totBr));
        tptBr=((double)ttBrN/((double)totBrN))/(((double)tBr/(double)totBr)*((double)tBr/(double)totBr));
        
        //NonBridge
        apaNonBr=((double)aaNonBrN/((double)totNonBrN))/(((double)aNonBr/(double)totNonBr)*((double)aNonBr/(double)totNonBr));
        apcNonBr=((double)acNonBrN/((double)totNonBrN))/(((double)aNonBr/(double)totNonBr)*((double)cNonBr/(double)totNonBr));
        apgNonBr=((double)agNonBrN/((double)totNonBrN))/(((double)aNonBr/(double)totNonBr)*((double)gNonBr/(double)totNonBr));
        aptNonBr=((double)atNonBrN/((double)totNonBrN))/(((double)aNonBr/(double)totNonBr)*((double)tNonBr/(double)totNonBr));
        
        cpaNonBr=((double)caNonBrN/((double)totNonBrN))/(((double)cNonBr/(double)totNonBr)*((double)aNonBr/(double)totNonBr));
        cpcNonBr=((double)ccNonBrN/((double)totNonBrN))/(((double)cNonBr/(double)totNonBr)*((double)cNonBr/(double)totNonBr));
        cpgNonBr=((double)cgNonBrN/((double)totNonBrN))/(((double)cNonBr/(double)totNonBr)*((double)gNonBr/(double)totNonBr));
        cptNonBr=((double)ctNonBrN/((double)totNonBrN))/(((double)cNonBr/(double)totNonBr)*((double)tNonBr/(double)totNonBr));
        
        gpaNonBr=((double)gaNonBrN/((double)totNonBrN))/(((double)gNonBr/(double)totNonBr)*((double)aNonBr/(double)totNonBr));
        gpcNonBr=((double)gcNonBrN/((double)totNonBrN))/(((double)gNonBr/(double)totNonBr)*((double)cNonBr/(double)totNonBr));
        gpgNonBr=((double)ggNonBrN/((double)totNonBrN))/(((double)gNonBr/(double)totNonBr)*((double)gNonBr/(double)totNonBr));
        gptNonBr=((double)gtNonBrN/((double)totNonBrN))/(((double)gNonBr/(double)totNonBr)*((double)tNonBr/(double)totNonBr));
        
        tpaNonBr=((double)taNonBrN/((double)totNonBrN))/(((double)tNonBr/(double)totNonBr)*((double)aNonBr/(double)totNonBr));
        tpcNonBr=((double)tcNonBrN/((double)totNonBrN))/(((double)tNonBr/(double)totNonBr)*((double)cNonBr/(double)totNonBr));
        tpgNonBr=((double)tgNonBrN/((double)totNonBrN))/(((double)tNonBr/(double)totNonBr)*((double)gNonBr/(double)totNonBr));
        tptNonBr=((double)ttNonBrN/((double)totNonBrN))/(((double)tNonBr/(double)totNonBr)*((double)tNonBr/(double)totNonBr));
        
        if(Double.isNaN(apaBr))
            apaBr=0;
        if(Double.isNaN(apcBr))
            apcBr=0;
        if(Double.isNaN(apgBr))
            apgBr=0;
        if(Double.isNaN(aptBr))
            aptBr=0;
        if(Double.isNaN(cpaBr))
            cpaBr=0;
        if(Double.isNaN(cpcBr))
            cpcBr=0;
        if(Double.isNaN(cpgBr))
            cpgBr=0;
        if(Double.isNaN(cptBr))
            cptBr=0;
        if(Double.isNaN(gpaBr))
            gpaBr=0;
        if(Double.isNaN(gpcBr))
            gpcBr=0;
        if(Double.isNaN(gpgBr))
            gpgBr=0;
        if(Double.isNaN(gptBr))
            gptBr=0;
        if(Double.isNaN(tpaBr))
            tpaBr=0;
        if(Double.isNaN(tpcBr))
            tpcBr=0;
        if(Double.isNaN(tpgBr))
            tpgBr=0;
        if(Double.isNaN(tptBr))
            tptBr=0;
        if(Double.isNaN(apaNonBr))
            apaNonBr=0;
        if(Double.isNaN(apcNonBr))
            apcNonBr=0;
        if(Double.isNaN(apgNonBr))
            apgNonBr=0;
        if(Double.isNaN(aptNonBr))
            aptNonBr=0;
        if(Double.isNaN(cpaNonBr))
            cpaNonBr=0;
        if(Double.isNaN(cpcNonBr))
            cpcNonBr=0;
        if(Double.isNaN(cpgNonBr))
            cpgNonBr=0;
        if(Double.isNaN(cptNonBr))
            cptNonBr=0;
        if(Double.isNaN(gpaNonBr))
            gpaNonBr=0;
        if(Double.isNaN(gpcNonBr))
            gpcNonBr=0;
        if(Double.isNaN(gpgNonBr))
            gpgNonBr=0;
        if(Double.isNaN(gptNonBr))
            gptNonBr=0;
        if(Double.isNaN(tpaNonBr))
            tpaNonBr=0;
        if(Double.isNaN(tpcNonBr))
            tpcNonBr=0;
        if(Double.isNaN(tpgNonBr))
            tpgNonBr=0;
        if(Double.isNaN(tptNonBr))
            tptNonBr=0;
        
        if(Double.isNaN(apa))
            apa=0;
        if(Double.isNaN(apc))
            apc=0;
        if(Double.isNaN(apg))
            apg=0;
        if(Double.isNaN(apt))
            apt=0;
        if(Double.isNaN(cpa))
            cpa=0;
        if(Double.isNaN(cpc))
            cpc=0;
        if(Double.isNaN(cpg))
            cpg=0;
        if(Double.isNaN(cpt))
            cpt=0;
        if(Double.isNaN(gpa))
            gpa=0;
        if(Double.isNaN(gpc))
            gpc=0;
        if(Double.isNaN(gpg))
            gpg=0;
        if(Double.isNaN(gpt))
            gpt=0;
        if(Double.isNaN(tpa))
            tpa=0;
        if(Double.isNaN(tpc))
            tpc=0;
        if(Double.isNaN(tpg))
            tpg=0;
        if(Double.isNaN(tpt))
            tpt=0;
        if(Double.isNaN(upa))
            upa=0;
        
        rawAA=(double)aaN/((double)totLength-(double)totSeqs);
        rawAC=(double)acN/((double)totLength-(double)totSeqs);
        rawAG=(double)agN/((double)totLength-(double)totSeqs);
        rawAT=(double)atN/((double)totLength-(double)totSeqs);
        rawCA=(double)caN/((double)totLength-(double)totSeqs);
        rawCC=(double)ccN/((double)totLength-(double)totSeqs);
        rawCG=(double)cgN/((double)totLength-(double)totSeqs);
        rawCT=(double)ctN/((double)totLength-(double)totSeqs);
        rawGA=(double)gaN/((double)totLength-(double)totSeqs);
        rawGC=(double)gcN/((double)totLength-(double)totSeqs);
        rawGG=(double)ggN/((double)totLength-(double)totSeqs);
        rawGT=(double)gtN/((double)totLength-(double)totSeqs);
        rawTA=(double)taN/((double)totLength-(double)totSeqs);
        rawTC=(double)tcN/((double)totLength-(double)totSeqs);
        rawTG=(double)tgN/((double)totLength-(double)totSeqs);
        rawTT=(double)ttN/((double)totLength-(double)totSeqs);
                
        //Count codons and AAs
        for(int i=0;i<seq.length();i+=3) {
            //not really needed as seqs not multiple of 3 are kicked out above
            if((i+3)>seq.length())
                break;
            
            String cod=""+seq.charAt(i)+seq.charAt(i+1)+seq.charAt(i+2);
            
            if(codLookUp.get(cod)==null)
                nullCod++;
            else {
                codCounts[codLookUp.get(cod)]++;
                aaCounts[aaLookUp.get(aaCode[1][codLookUp.get(cod)])]++;
            }
            
            totCod++;
        }
        
        //AA bias
        for(int i=0;i<aaCounts.length;i++) {
            aaBias[i]=(double)aaCounts[i]/(double)totCod;
        }
        
        //Codon Bias
        for(int i=0;i<codCounts.length;i++) {
            if(codCounts[i]>0 & aaCounts[aaLookUp.get(aaCode[1][i])]==0)
                System.out.println("Error "+i+" codCounts>0 but AA=0 codCount="+codCounts[i]+" "+aaLookUp.get(aaCode[1][i]));
            
            if(aaCounts[aaLookUp.get(aaCode[1][i])]>0)
                codBias[i]=(double)codCounts[i]/(double)aaCounts[aaLookUp.get(aaCode[1][i])];
        }
        
        //Codon Pair & AA Pair counts
        for(int i=3;i<seq.length();i+=3) {
            if((i+3)>seq.length())
                break;
            
            String cod1=""+seq.charAt(i-3)+seq.charAt(i-2)+seq.charAt(i-1);
            String cod2=""+seq.charAt(i)+seq.charAt(i+1)+seq.charAt(i+2);
            
            //ignore stops at codon 1 - even if readthrough
            if(cod1.equals("TAG")|cod1.equals("TAA")|cod1.equals("TGA"))
                stops++;
            
            //ignore codon pairs with Ns or ambiguitites
            if(codLookUp.get(cod1)!=null & codLookUp.get(cod2)!=null) {
                codPairCounts[codLookUp.get(cod1)][codLookUp.get(cod2)]++;
                aaPairCounts[aaLookUp.get(aaCode[1][codLookUp.get(cod1)])][aaLookUp.get(aaCode[1][codLookUp.get(cod2)])]++;
            }
            
            totCP++;
        }
        
        //Codon Pair Bias raw
        for(int i=0;i<cpBias.length;i++)
            for(int j=0;j<cpBias[i].length;j++)
                if(codPairCounts[i][j]>0)
                    cpBias[i][j]=(double)codPairCounts[i][j]/(double)aaPairCounts[aaLookUp.get(aaCode[1][i])][aaLookUp.get(aaCode[1][j])];
                else
                    cpBias[i][j]=0;
        
        
        //CPB = ln [codonPairCount / [[(codon1Count x codon2Count) / (aa1Count * aa2Count)] x aaPairCount]]
        //Codon Pair Bias coleman formula
        for(int i=0;i<cpb.length;i++) {
            for(int j=0;j<cpb[i].length;j++) {
                if(codPairCounts[i][j]>0) {
                    cpb[i][j]=Math.log((double)codPairCounts[i][j]/((((double)codCounts[i]*(double)codCounts[j])/(aaCounts[aaLookUp.get(aaCode[1][i])]*aaCounts[aaLookUp.get(aaCode[1][j])]))*aaPairCounts[aaLookUp.get(aaCode[1][i])][aaLookUp.get(aaCode[1][j])]));
                    cpbSum+=cpb[i][j];
                    cpbCount++;
                }
            }
        }
        
        cpbAv=cpbSum/(double)cpbCount;
        
        for(int i=0;i<cpb.length;i++) {
            for(int j=0;j<cpb[i].length;j++) {
                if(cpb[i][j]<trueCpbMin & cpb[i][j]!=cpbMin)
                    trueCpbMin=cpb[i][j];
            }
        }
        
        //not used anymore
        if(trueCpbMin<0)
            cpbMin=trueCpbMin*2;
        else
            cpbMin=trueCpbMin/2;
        
        for(int i=0;i<cpb.length;i++) {
            for(int j=0;j<cpb[i].length;j++) {
                //if(cpb[i][j]==0)
                    //cpb[i][j]=cpbAv;
                
                if(codPairCounts[i][j]==0) {
                    //if codon pair count is 0, and the AA pair count is also zero - set to the cpbAv
                    if(aaPairCounts[aaLookUp.get(aaCode[1][i])][aaLookUp.get(aaCode[1][j])]==0) {
                        //cpb[i][j]=-100;
                        cpb[i][j]=cpbAv;
                    }
                    //if codon pair count is 0, but there are AA pair counts - set to -9999
                    else {
                        cpb[i][j]=-9999;
                        //cpb[i][j]=cpbMin;
                    }
                }
            }
        }
        
        
    }
    public static void createCode() {
        
        aminoAcids=new String[2][aminoNum];//stores the AA codes and names
        aminoAcids[0][0]="L";
        aminoAcids[1][0]="leu";
        
        aminoAcids[0][1]="P";
        aminoAcids[1][1]="pro";
        
        aminoAcids[0][2]="H";
        aminoAcids[1][2]="his";
        
        aminoAcids[0][3]="Q";
        aminoAcids[1][3]="gln";
        
        aminoAcids[0][4]="R";
        aminoAcids[1][4]="arg";
        
        aminoAcids[0][5]="I";
        aminoAcids[1][5]="ile";

        aminoAcids[0][6]="M";
        aminoAcids[1][6]="met";//start
        
        aminoAcids[0][7]="T";
        aminoAcids[1][7]="thr";
        
        aminoAcids[0][8]="N";
        aminoAcids[1][8]="asn";
        
        aminoAcids[0][9]="K";
        aminoAcids[1][9]="lys";
        
        aminoAcids[0][10]="S";
        aminoAcids[1][10]="ser";
      
        aminoAcids[0][11]="V";
        aminoAcids[1][11]="val";
        
        aminoAcids[0][12]="A";
        aminoAcids[1][12]="ala";
        
        aminoAcids[0][13]="D";
        aminoAcids[1][13]="asp";
        
        aminoAcids[0][14]="E";
        aminoAcids[1][14]="glu";
        
        aminoAcids[0][15]="G";
        aminoAcids[1][15]="gly";
        
        aminoAcids[0][16]="F";
        aminoAcids[1][16]="phe";
        
        aminoAcids[0][17]="Y";
        aminoAcids[1][17]="tyr";

        aminoAcids[0][18]="C";
        aminoAcids[1][18]="cys";
        
        aminoAcids[0][19]="W";
        aminoAcids[1][19]="trp";
        
        aminoAcids[0][20]="X";
        aminoAcids[1][20]="stp";//stop
        
        //check for duplicate AAs - as list above was manually defined
        for(int i=0;i<aminoAcids[0].length;i++) {
            for(int j=i+1;j<aminoAcids[0].length;j++) {
                if(aminoAcids[0][i].equalsIgnoreCase(aminoAcids[0][j])) {
                    System.out.println("Duplicate AAs - "+i+ " "+j+" - "+aminoAcids[0][i]+" and "+aminoAcids[0][j]);
                }
                
                if(aminoAcids[1][i].equalsIgnoreCase(aminoAcids[1][j])) {
                    System.out.println("Duplicate AAs - "+i+ " "+j+" - "+aminoAcids[1][i]+" and "+aminoAcids[1][j]);
                }
            }
        }
        
        aaCode=new String[2][64];
        
        aaCode[0][0]="AAA";
        aaCode[1][0]="K";
        aaCode[0][1]="AAC";
        aaCode[1][1]="N";
        aaCode[0][2]="AAG";
        aaCode[1][2]="K";
        aaCode[0][3]="AAT";
        aaCode[1][3]="N";
        
        aaCode[0][4]="ACA";
        aaCode[1][4]="T";
        aaCode[0][5]="ACC";
        aaCode[1][5]="T";
        aaCode[0][6]="ACG";
        aaCode[1][6]="T";
        aaCode[0][7]="ACT";
        aaCode[1][7]="T";
        
        aaCode[0][8]="AGA";
        aaCode[1][8]="R";
        aaCode[0][9]="AGC";
        aaCode[1][9]="S";
        aaCode[0][10]="AGG";
        aaCode[1][10]="R";
        aaCode[0][11]="AGT";
        aaCode[1][11]="S";
        
        aaCode[0][12]="ATA";
        aaCode[1][12]="I";
        aaCode[0][13]="ATC";
        aaCode[1][13]="I";
        aaCode[0][14]="ATG";
        aaCode[1][14]="M";//start
        aaCode[0][15]="ATT";
        aaCode[1][15]="I";
        
        aaCode[0][16]="CAA";
        aaCode[1][16]="Q";
        aaCode[0][17]="CAC";
        aaCode[1][17]="H";
        aaCode[0][18]="CAG";
        aaCode[1][18]="Q";
        aaCode[0][19]="CAT";
        aaCode[1][19]="H";
        
        aaCode[0][20]="CCA";
        aaCode[1][20]="P";
        aaCode[0][21]="CCC";
        aaCode[1][21]="P";
        aaCode[0][22]="CCG";
        aaCode[1][22]="P";
        aaCode[0][23]="CCT";
        aaCode[1][23]="P";
        
        aaCode[0][24]="CGA";
        aaCode[1][24]="R";
        aaCode[0][25]="CGC";
        aaCode[1][25]="R";
        aaCode[0][26]="CGG";
        aaCode[1][26]="R";
        aaCode[0][27]="CGT";
        aaCode[1][27]="R";
        
        aaCode[0][28]="CTA";
        aaCode[1][28]="L";
        aaCode[0][29]="CTC";
        aaCode[1][29]="L";
        aaCode[0][30]="CTG";
        aaCode[1][30]="L";
        aaCode[0][31]="CTT";
        aaCode[1][31]="L";
        
        aaCode[0][32]="GAA";
        aaCode[1][32]="E";
        aaCode[0][33]="GAC";
        aaCode[1][33]="D";
        aaCode[0][34]="GAG";
        aaCode[1][34]="E";
        aaCode[0][35]="GAT";
        aaCode[1][35]="D";
        
        aaCode[0][36]="GCA";
        aaCode[1][36]="A";
        aaCode[0][37]="GCC";
        aaCode[1][37]="A";
        aaCode[0][38]="GCG";
        aaCode[1][38]="A";
        aaCode[0][39]="GCT";
        aaCode[1][39]="A";
        
        aaCode[0][40]="GGA";
        aaCode[1][40]="G";
        aaCode[0][41]="GGC";
        aaCode[1][41]="G";
        aaCode[0][42]="GGG";
        aaCode[1][42]="G";
        aaCode[0][43]="GGT";
        aaCode[1][43]="G";
        
        aaCode[0][44]="GTA";
        aaCode[1][44]="V";
        aaCode[0][45]="GTC";
        aaCode[1][45]="V";
        aaCode[0][46]="GTG";
        aaCode[1][46]="V";
        aaCode[0][47]="GTT";
        aaCode[1][47]="V";
        
        aaCode[0][48]="TAA";
        aaCode[1][48]="X";//stop
        aaCode[0][49]="TAC";
        aaCode[1][49]="Y";
        aaCode[0][50]="TAG";
        aaCode[1][50]="X";//stop
        aaCode[0][51]="TAT";
        aaCode[1][51]="Y";
        
        aaCode[0][52]="TCA";
        aaCode[1][52]="S";
        aaCode[0][53]="TCC";
        aaCode[1][53]="S";
        aaCode[0][54]="TCG";
        aaCode[1][54]="S";
        aaCode[0][55]="TCT";
        aaCode[1][55]="S";
        
        aaCode[0][56]="TGA";
        aaCode[1][56]="X";//stop
        aaCode[0][57]="TGC";
        aaCode[1][57]="C";
        aaCode[0][58]="TGG";
        aaCode[1][58]="W";//tryp - only codon for tryp - not essential then?
        aaCode[0][59]="TGT";
        aaCode[1][59]="C";
        
        aaCode[0][60]="TTA";
        aaCode[1][60]="L";
        aaCode[0][61]="TTC";
        aaCode[1][61]="F";
        aaCode[0][62]="TTG";
        aaCode[1][62]="L";
        aaCode[0][63]="TTT";
        aaCode[1][63]="F";
        
        for(int i=0;i<aaCode[0].length;i++) {
            
            //check for blanks
            if((aaCode[0][i].equalsIgnoreCase(""))) {
                System.out.println("Blank codon - "+i+" - "+aaCode[0][i]+" "+aaCode[1][i]);
            }
            
            boolean aaTest=false;
            //check that AAcodes are correct
            for(int j=0;j<aminoAcids[0].length;j++) {
                if(aaCode[1][i].equalsIgnoreCase(aminoAcids[0][j])) {
                    aaTest=true;
                    break;
                }
            }
            if(!aaTest)
                System.out.println("Codon AA not found in aminoAcids- "+i+" - "+aaCode[0][i]+" "+aaCode[1][i]);
            
            //check codons are unique
            for(int j=i+1;j<aaCode[0].length;j++) {
                if(aaCode[0][i].equalsIgnoreCase(aaCode[0][j])) {
                    System.out.println("Duplicate codons - "+i+" "+j+" - "+aaCode[0][i]+" "+aaCode[0][j]);
                }
            }
        }
        
        
        //check each amino acid appears in codons
        for(int i=0;i<aminoAcids[0].length;i++) {
            boolean aaTest=false;
            
            for(int j=0;j<aaCode[1].length;j++) {
                if(aminoAcids[0][i].equalsIgnoreCase(aaCode[1][j])) {
                    aaTest=true;
                    break;
                }
            }
            
            if(!aaTest) {
                System.out.println("aminoAcid not found in codons - "+i+" - "+aminoAcids[0][i]+" "+aminoAcids[1][i]);
            }
        }
        
        //create AA index - where are they located in the arrays
        aaLookUp=new HashMap<String,Integer>();
        for(int i=0;i<aminoAcids[0].length;i++) {
            aaLookUp.put(aminoAcids[0][i], i);
        }
        
        codLookUp=new HashMap<String,Integer>();
        for(int i=0;i<aaCode[0].length;i++) {
            codLookUp.put(aaCode[0][i], i);
        }
        
        //For each AA - see how many codons code it - and populate their data
        aminoCodons=new String[aminoNum][10];//ten is arbitrary and excess - maximum number of codons that an individual AA has
        for(int i=0;i<aminoCodons.length;i++) {
            for(int j=0;j<aminoCodons[i].length;j++) {
                aminoCodons[i][j]="";
            }
        }
        aminoCodonCounts=new int[aminoNum];//stores the number of codons each AA has
        for(int i=0;i<aminoAcids[0].length;i++) {

            count=0;
            
            for(int j=0;j<aaCode[1].length;j++) {
                if(aminoAcids[0][i].equalsIgnoreCase(aaCode[1][j])) {
                    aminoCodons[i][count]=aaCode[0][j];
                    count++;
                }
            }
            
            aminoCodonCounts[i]=count;
        }
 
    }   
    
}
