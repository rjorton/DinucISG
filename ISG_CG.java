
package isg_cg;

import java.io.*;
import java.util.HashMap;
import java.util.Random;

public class ISG_CG {

    public static String filename, outname;
    public static String seqName, prevName, seq, gene, trans;
    public static int count=0, miss=0;
    public static int a,c,g,t,n;
    public static int aaN,acN,agN,atN,caN,ccN,cgN,ctN,gaN,gcN,ggN,gtN,taN,tcN,tgN,ttN;
    public static double apa,apc,apg,apt,cpa,cpc,cpg,cpt,gpa,gpc,gpg,gpt,tpa,tpc,tpg,tpt;
    public static double gc, at;
    public static int totSeqs, totLength;
    
    public static double aaCon,acCon,agCon,atCon,caCon,ccCon,cgCon,ctCon;
    public static double gaCon,gcCon,ggCon,gtCon,taCon,tcCon,tgCon,ttCon;
    
    public static void main(String[] args) {

        if(args.length==1) {
            filename=args[0];
        }
        else {
            System.out.println("Error - incorrect parameters "+args.length);
            filename="/Users/rjorton/Zap/coding_human_cds_new.fasta";
            filename="/Users/rjorton/Dropbox/SGM/flu/PB2_new.fasta";
            filename="/Users/rjorton/Dropbox/human_cdna_new.fasta";
        }

        outname=filename.substring(0, filename.lastIndexOf("."))+"_dat.txt";
   
        try{
            FileWriter fstream = new FileWriter(outname);
            BufferedWriter out = new BufferedWriter(fstream);
        
            out.write(genHeader());
            
            try {
                File seqFile = new File(filename);
                BufferedReader input =  new BufferedReader(new FileReader(seqFile));

                try {
                    String line = null;
      
                    seqName="START";
                    prevName="";
                    count=0;
                    miss=0;
                    while (( line = input.readLine()) != null) {
                        
                        if(line.indexOf(">")==0) {
                            //>ENSG00000002587|ENST00000002596
                            //Added in check for Salmon as no |
                            if(line.indexOf("|")>0) {
                                gene=line.substring(1,line.indexOf("|"));
                                trans=line.substring(line.indexOf("|")+1,line.length());
                            }
                            else {
                                gene=line.substring(1,line.length());
                                trans=gene;
                            }
                            
                            seqName=line;
                            
                            if(prevName.equals(seqName))
                                System.out.println("Duplicate = "+seqName);
                            
                            prevName=seqName;
                        }
                        else {   
                            seq=line.toUpperCase();
                            
                            boolean test=true;
                            if(seq.equalsIgnoreCase("sequence unavailable")) {
                                seq="";
                                test=false;
                            }
                                
                            if(seq.length()%3!=0)
                                System.out.println(totSeqs+" Error - seq length not divisible by 3 "+prevName+" "+seq.length());
                            
                            clearDat();
                            analyseSeq(seq);
                            if(test)
                                out.write(genOutput()+"\n");
                            else {
                                System.out.println("Skipping sequence unavailble: "+seqName);
                                miss++;
                            }
                            
                            
                        }
                        
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
        
            out.close();
        }

        catch (Exception e) {
            e.printStackTrace();
            System.err.println("Error: " + e.getMessage());
        }
        
        System.out.println("TotSeqs="+totSeqs);
        System.out.println("Skipped="+miss);

    }
    
    public static void clearDat() {
        a=c=g=t=n=0;
        aaN=acN=agN=atN=caN=ccN=cgN=ctN=gaN=gcN=ggN=gtN=taN=tcN=tgN=ttN=0;
        apa=apc=apg=apt=cpa=cpc=cpg=cpt=gpa=gpc=gpg=gpt=tpa=tpc=tpg=tpt=0;
        cgCon=taCon=0;
        
        aaCon=acCon=agCon=atCon=caCon=ccCon=cgCon=ctCon=0;
        gaCon=gcCon=ggCon=gtCon=taCon=tcCon=tgCon=ttCon=0;
        
        gc=at=0;
        totLength=0;
        
    }
    public static String genHeader() {
        
        String head="";
        head+="Num\tName\tGene\tTranscript\ttotLength\tCpG\tCG%\tUpA\tUA%\ta\tc\tg\tt\tn";
        head+="\tApA\tApC\tApG\tApT\tCpA\tCpC\tCpT\tGpA\tGpC\tGpG\tGpT\tUpC\tUpG\tUpT";
        head+="\tGCcont\tATcont\tAA%\tAC%\tAG%\tAT%\tCA%\tCC%\tCT%\tGA%\tGC%\tGG%\tGT%\tUC%\tUG%\tUT%\n";
        
        return head;
    }
    
    public static String genOutput() {
        
        String dat="";
 
        dat=totSeqs+"\t"+seqName+"\t"+gene+"\t"+trans+"\t"+totLength+"\t"+cpg+"\t"+cgCon+"\t"+tpa+"\t"+taCon;
        dat+="\t"+a+"\t"+c+"\t"+g+"\t"+t+"\t"+n;
        dat+="\t"+apa+"\t"+apc+"\t"+apg+"\t"+apt;
        dat+="\t"+cpa+"\t"+cpc+"\t"+cpt;
        dat+="\t"+gpa+"\t"+gpc+"\t"+gpg+"\t"+gpt;
        dat+="\t"+tpc+"\t"+tpg+"\t"+tpt;
        dat+="\t"+gc+"\t"+at;
        dat+="\t"+aaCon+"\t"+acCon+"\t"+agCon+"\t"+atCon;
        dat+="\t"+caCon+"\t"+ccCon+"\t"+ctCon;
        dat+="\t"+gaCon+"\t"+gcCon+"\t"+ggCon+"\t"+gtCon;
        dat+="\t"+tcCon+"\t"+tgCon+"\t"+ttCon;

        return dat;

    }
    
    public static void analyseSeq(String seq) {

        totLength+=seq.length();
        totSeqs++;
        
        for(int i=0;i<seq.length();i++) {
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
        }
        
        gc=(double)(g+c)/(double)totLength;
        at=(double)(a+t)/(double)totLength;
        
        for(int i=1;i<seq.length();i++) {
   
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
        }
        
        //cgCon=(double)(cgN)/((double)totLength-1);
        //taCon=(double)(taN)/((double)totLength-1);
        
        aaCon=(double)(aaN)/((double)totLength-1);
        acCon=(double)(acN)/((double)totLength-1);
        agCon=(double)(agN)/((double)totLength-1);
        atCon=(double)(atN)/((double)totLength-1);
        caCon=(double)(caN)/((double)totLength-1);
        ccCon=(double)(ccN)/((double)totLength-1);
        cgCon=(double)(cgN)/((double)totLength-1);
        ctCon=(double)(ctN)/((double)totLength-1);
        gaCon=(double)(gaN)/((double)totLength-1);
        gcCon=(double)(gcN)/((double)totLength-1);
        ggCon=(double)(ggN)/((double)totLength-1);
        gtCon=(double)(gtN)/((double)totLength-1);
        taCon=(double)(taN)/((double)totLength-1);
        tcCon=(double)(tcN)/((double)totLength-1);
        tgCon=(double)(tgN)/((double)totLength-1);
        ttCon=(double)(ttN)/((double)totLength-1);
        
        
        //ERROR - SHOULD be totLength-1!!!!!!!! all these below were totLength
        apa=((double)aaN/((double)totLength-1))/(((double)a/(double)totLength)*((double)a/(double)totLength));
        apc=((double)acN/((double)totLength-1))/(((double)a/(double)totLength)*((double)c/(double)totLength));
        apg=((double)agN/((double)totLength-1))/(((double)a/(double)totLength)*((double)g/(double)totLength));
        apt=((double)atN/((double)totLength-1))/(((double)a/(double)totLength)*((double)t/(double)totLength));
        
        cpa=((double)caN/((double)totLength-1))/(((double)c/(double)totLength)*((double)a/(double)totLength));
        cpc=((double)ccN/((double)totLength-1))/(((double)c/(double)totLength)*((double)c/(double)totLength));
        cpg=((double)cgN/((double)totLength-1))/(((double)c/(double)totLength)*((double)g/(double)totLength));
        cpt=((double)ctN/((double)totLength-1))/(((double)c/(double)totLength)*((double)t/(double)totLength));
        
        gpa=((double)gaN/((double)totLength-1))/(((double)g/(double)totLength)*((double)a/(double)totLength));
        gpc=((double)gcN/((double)totLength-1))/(((double)g/(double)totLength)*((double)c/(double)totLength));
        gpg=((double)ggN/((double)totLength-1))/(((double)g/(double)totLength)*((double)g/(double)totLength));
        gpt=((double)gtN/((double)totLength-1))/(((double)g/(double)totLength)*((double)t/(double)totLength));
        
        tpa=((double)taN/((double)totLength-1))/(((double)t/(double)totLength)*((double)a/(double)totLength));
        tpc=((double)tcN/((double)totLength-1))/(((double)t/(double)totLength)*((double)c/(double)totLength));
        tpg=((double)tgN/((double)totLength-1))/(((double)t/(double)totLength)*((double)g/(double)totLength));
        tpt=((double)ttN/((double)totLength-1))/(((double)t/(double)totLength)*((double)t/(double)totLength));

    }
        
}
