################################################
#File Name: run.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Tue 15 Feb 2022 08:33:09 PM CST
################################################

#!/bin/sh 


<<'!'
### test ###
    echo "G7N8T9" | while read id;do
    [ -f protein_fasta/$id.fasta ] && rm -rf protein_fasta/$id.fasta
    wget https://www.uniprot.org/uniprot/$id.fasta -P protein_fasta/ 
    grep -v '^>' protein_fasta/$id.fasta > protein_fasta/$id.rmtil.fasta
    done

    blat -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 -q=prot -t=dnax  \
    /home/user/data/lit/database/public/genome/rheMac10/rheMac10.fa.gz \
    protein_fasta/G7N8T9.fasta output/G7N8T9.psl

    blat -q=prot -t=dnax  \
    /home/user/data/lit/database/public/genome/rheMac10/rheMac10.fa.gz \
    protein_fasta/G7N8T9.fasta output/G7N8T9.default.psl

    blat -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 -q=prot -t=dnax  \
    /home/user/data/lit/database/public/genome/rheMac10/mrna.fa.gz \
    protein_fasta/G7N8T9.fasta output/G7N8T9.mrna.psl

    # default
    # -minScore=N Default is 30.
    # -minIdentity=N  Default is 90 for nucleotide searches, 25 for protein or translated protein searches.
    blat -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 -q=prot -t=dnax  \
    /home/user/data/lit/database/public/genome/hg38/hg38.fa.gz \
    protein_fasta/G7N8T9.fasta output/G7N8T9.hg38.psl 

    # default
    # -minScore=N Default is 30.
    # -minIdentity=N  Default is 90 for nucleotide searches, 25 for protein or translated protein searches.
    blat -q=prot -t=dnax  \
    /home/user/data/lit/database/public/genome/hg38/hg38.fa.gz \
    protein_fasta/G7N8T9.fasta output/G7N8T9.hg38.psl 

    pslScore.pl output/G7N8T9.default.psl  1> output/G7N8T9.default.psl.score.txt 2> /dev/null
    less output/G7N8T9.default.psl.score.txt | grep '^chr' | sort -k6,6nr -k5,5nr | awk -v OFS='\t' '($6==100 && NR <= 3){print $1,$2,$3}' | \
    bedtools intersect -wa -wb -a /home/user/data/lit/database/public/annotation/gene_and_gene_predictions/rheMac10.refGene.clean.bed \
    -b -

    pslScore.pl output/G7N8T9.hg38.psl  1> output/G7N8T9.hg38.psl.score.txt 2> /dev/null
    less output/G7N8T9.hg38.psl.score.txt | grep '^chr' | sort -k6,6nr -k5,5nr | awk -v OFS='\t' '($6>=80){print $1,$2,$3,$5,$6}' | \
    bedtools intersect -wo -a /home/user/data/lit/database/public/annotation/gene_and_gene_predictions/hg38.refGene.clean.bed \
    -b - | awk -v OFS='\t' 'BEGIN{print "ID\tgeneName\tscore\tidentity\tintersect_number"} {print "'$id'",$4,$8,$9,$10}'

    # pslScore.pl output/G7N8T9.psl 1> output/G7N8T9.psl.score.txt 2> /dev/null
    # less output/G7N8T9.psl.score.txt | grep '^chr' | sort -k6,6nr -k5,5nr | awk -v OFS='\t' '($6>=80 && NR <= 3){print $1,$2,$3}' | \
    # bedtools intersect -wa -wb -a /home/user/data/lit/database/public/annotation/gene_and_gene_predictions/rheMac10.refGene.clean.bed \
    # -b -

    less /home/user/data/lit/database/public/annotation/gene_and_gene_predictions/hg38.refGene.bed.gz | \
    awk -v OFS='\t' 'NR>1 {print $3,$5,$6,$13}' > /home/user/data/lit/database/public/annotation/gene_and_gene_predictions/hg38.refGene.clean.bed

    less /home/user/data/lit/database/public/annotation/gene_and_gene_predictions/rheMac10.refGene.bed | \
    awk -v OFS='\t' 'NR>1 {print $3,$5,$6,$13}' > /home/user/data/lit/database/public/annotation/gene_and_gene_predictions/rheMac10.refGene.clean.bed

    less rheMac10.ensembl.gene.stable.csv | sed '1d' | \
    awk  -v OFS='\t' -v FS=',' '$5=="" {print "chr"$1,$2,$3,$4,"NULL"} $5!=""{print "chr"$1,$2,$3,$4,$5}' > \
    /home/user/data/lit/database/public/annotation/gene_and_gene_predictions/rheMac10.ensembl.gene.stable.bed

    for spe in hg38 rheMac10;do
        less $spe.ensembl.gene.stable.csv | sed '1d' | \
        awk  -v OFS='\t' -v FS=',' '{print "chr"$1,$2,$3,$4}' > \
        /home/user/data/lit/database/public/annotation/gene_and_gene_predictions/$spe.ensembl.gene.stable.bed
    done

    for spe in hg38 rheMac10;do
        less $spe.ensembl.gene.stable.csv | sed '1d' | \
        awk  -v OFS='\t' -v FS=',' '$5=="" {print "chr"$1,$2,$3,$4":""NULL"} $5!=""{print "chr"$1,$2,$3,$4":"$5}' > \
        /home/user/data/lit/database/public/annotation/gene_and_gene_predictions/$spe.ensembl.gene.stable.bed
    done

### test ###
!

<<'!'
### test: parasol multi-task ###

    cat EGK.txt | \
    head -n 3 | \
    while read id;do
        dir=/home/user/data2/lit/project/help/proteinID2geneName
        echo "blat -q=prot -t=dnax  \
        /home/user/data/lit/database/public/genome/hg38/hg38.fa.gz \
        $dir/protein_fasta/$id.fasta $dir/output/psl/$id.hg38.psl 
        pslScore.pl $dir/output/psl/$id.hg38.psl  1> $dir/output/score/$id.hg38.psl.score.txt 2> /dev/null

        less $dir/output/score/$id.hg38.psl.score.txt | grep '^chr' | sort -k6,6nr -k5,5nr | awk -v OFS='\t' '(\$6>=80){print \$1,\$2,\$3,\$5,\$6}' | \
        bedtools intersect -wo -a /home/user/data/lit/database/public/annotation/gene_and_gene_predictions/hg38.ensembl.gene.stable.bed \
        -b - | awk -v OFS='\t' '{print "'$id'",\$4,\$5,\$9,\$10,\$11}' >> $dir/res.hg38.txt

        blat -q=prot -t=dnax  \
        /home/user/data/lit/database/public/genome/rheMac10/rheMac10.fa.gz \
        $dir/protein_fasta/$id.fasta $dir/output/psl/$id.psl 
        pslScore.pl $dir/output/psl/$id.psl  1> $dir/output/score/$id.psl.score.txt 2> /dev/null

        less $dir/output/score/$id.psl.score.txt | grep '^chr' | sort -k6,6nr -k5,5nr | awk -v OFS='\t' '(\$6>=80){print \$1,\$2,\$3,\$5,\$6}' | \
        bedtools intersect -wo -a /home/user/data/lit/database/public/annotation/gene_and_gene_predictions/rheMac10.ensembl.gene.stable.bed \
        -b - | awk -v OFS='\t' '{print "'$id'",\$4,\$5,\$9,\$10,\$11}' >> $dir/res.rheMac10.txt" >> proteinID2geneName.Joblst
    done

    cat EGK.txt | \
    head -n 3 | \
    while read id;do
        dir=/home/user/data2/lit/project/help/proteinID2geneName
        echo "/home/user/data/lit/Tools/ucsc/bin/blat/blat -q=prot -t=dnax /home/user/data/lit/database/public/genome/hg38/hg38.fa.gz $dir/protein_fasta/$id.fasta $dir/output/psl/$id.psl" >> test.Joblst
    done

    para stop
    para flushResults
    para resetCounts
    para freeBatch

    parasol list batches

    para create jobList   # Create job tracking database

    para try              # Run ten jobs

    para check            # Check up on jobs

    para push             # Push up to 100,000 jobs

    para check            # Check up on jobs

    para push             # Retry any crashed jobs so far

    para time             # Collect timing information.

### test:  parasol multi-tastk ###
!


# run start here

<<'!' 
### download protein fasta ###
    nohup cat EGK.txt | \
        # head -n 3 | \
        while read id;do
        [ -f protein_fasta/$id.fasta ] && rm -rf protein_fasta/$id.fasta
        echo "***download $id"
        wget https://www.uniprot.org/uniprot/$id.fasta -P protein_fasta/ 
    done >log/download.log 2>&1 &

    wget https://www.uniprot.org/uniprot/G7MFL7.fasta -P protein_fasta/ 
    wget https://www.uniprot.org/uniprot/G7MIN2.fasta -P protein_fasta/ 
    wget https://www.uniprot.org/uniprot/G7N5V9.fasta -P protein_fasta/ 

### download protein fasta ###
!

# several failed proteins failed to download(unnable to connect with network) are manually downloaded

<<'!' 
### download DNA database for blat ###
    echo -e "rheMac10\nhg38" | \
    # head -n1 | \
    while read spe;do 
        data_dir=/home/user/data/lit/database/public/genome/$spe
        [ -d $spe ] || mkdir $spe
        [ -d $spe ] && rm -rf $spe/.[a-z]*
        time (rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/$spe/bigZips/$spe.fa.gz $spe && \
        cp $spe/$spe.fa.gz $data_dir && rm -rf $spe/$spe.fa.gz) 
        time (rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/$spe/bigZips/mrna.fa.gz $spe && \
        cp $spe/mrna.fa.gz $data_dir && rm -rf $spe/mrna.fa.gz) 
        time (rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/$spe/bigZips/$spe.chrom.sizes $spe && \
        cp $spe/$spe.chrom.sizes $data_dir && rm -rf $spe/$spe.chrom.sizes) 
    done

    rm -rf hg38 rheMac10
### download DNA database for blat ###
!

<<'!'
# bedtools getfasta for blat database ###
    # $spe.ensembl.gene.stable.csv is downloaded from ensembl archive(version 105)
    for spe in hg38 rheMac10;do
    less $spe.ensembl.gene.stable.csv | sed '1d' | \
    awk  -v OFS='\t' -v FS=',' '$5=="" {print "chr"$1,$2,$3,$4":""NULL"} $5!=""{print "chr"$1,$2,$3,$4":"$5}' > \
    /home/user/data/lit/database/public/annotation/gene_and_gene_predictions/$spe.ensembl.gene.stable.bed
    done

    for spe in hg38 rheMac10;do
        bedtools getfasta -name -fi /home/user/data/lit/database/public/genome/$spe/$spe.fa -bed \
        /home/user/data/lit/database/public/annotation/gene_and_gene_predictions/$spe.ensembl.gene.stable.bed > $spe.ensembl.gene.stable.fasta
    done

    less hg38.ensembl.gene.stable.fasta | grep -v WARNING > hg38.ensembl.gene.stable.clean.fasta

    less rheMac10.ensembl.gene.stable.fasta | egrep -v 'WARNING|index' > rheMac10.ensembl.gene.stable.clean.fasta

# bedtools getfasta dor blat database ###
!


### blat ###
    [ -d output/psl/ ] || mkdir -p output/psl/
    [ -d output/score/ ] || mkdir -p output/score/

    for spe in hg38 rheMac10;do
        [ -f res.$spe.txt ] && rm -rf res.$spe.txt
        awk -v OFS='\t' 'BEGIN{print "ID\tgene\tscore\tidentity"}' > res.$spe.txt
    done

    cat EGK.txt | \
    # head -n 3 | \
    while read id;do
        for spe in hg38 rheMac10;do
            echo "*** `date "+%Y-%m-%d %H:%M:%S"` processing $id ref on $spe"

            time ( blat -q=prot -t=dnax  \
            $spe.ensembl.gene.stable.clean.fasta \
            protein_fasta/$id.fasta output/psl/$id.$spe.psl 
            pslScore.pl output/psl/$id.$spe.psl  1> output/score/$id.$spe.psl.score.txt 2> /dev/null )


            if [ $spe == "hg38" ];then
                    less output/score/$id.$spe.psl.score.txt | grep '^ENSG' | awk -v OFS='\t' '{print "'$id'",$1,$5,$6}' >> res.$spe.txt
                else
                    less output/score/$id.$spe.psl.score.txt | grep '^ENSMMUG' | awk -v OFS='\t' '{print "'$id'",$1,$5,$6}' >> res.$spe.txt
            fi
        
        done
    done
### blat ###
       
# nohup bash run.sh > log/run.log 2>&1 &

# subset for the best blat results
for spe in hg38 rheMac10;do
    less res.${spe}.txt | sed '1d' | awk -v OFS='\t' -v FS=':' '{print $1,$2,$3,$4,$5}' | awk -v OFS='\t' '{print $1,$2,$3,$6,$7}' | \
    sort -k1,1 -k4,4nr -k5,5nr > res.${spe}.reo.txt

    [ -f res.${spe}.unique.txt ] && rm -rf res.${spe}.unique.txt
    cat EGK.txt | \
    # head -n 3 | \
    while read id;do
        less res.${spe}.reo.txt | grep $id | head -n1 >> res.${spe}.unique.txt
    done
done

# sort
sort -k1,1 res.hg38.unique.txt > res.hg38.unique.sorted.txt
sort -k1,1 res.rheMac10.unique.txt > res.rheMac10.unique.sorted.txt

# join by protein id
join -a1 -a2 res.hg38.unique.sorted.txt res.rheMac10.unique.sorted.txt > res.join.txt

# protein id;ensembl id;geneName;geneName(complement)
less res.join.txt | awk '{if ($7=="NULL") print $1,$6,$7,$3; 
else if($7=="") print $1,"unmapped","unmapped",$3;
else print $1,$6,$7,$7}' > res.final.txt