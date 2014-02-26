#!bash

python strip_fastas.py -i ../fosmid_qc/data/EndSequenceDB/[A-Z][A-Z][A-Z][A-Z].fasta -o workspace/
for f in workspace/*txt
do
  echo perl id2plate3a.pl $f > $f.id2plate.txt
done

for file in workspace/*id2plate.txt
do 
   python process_id2plate3a.py -i $file -o workspace/ 
done

for file in ../fosmid_qc/data/EndSequenceDB/[A-Z][A-Z][A-Z][A-Z].fasta
do 
   file2=${file##\.\.\/fosmid_qc\/data\/EndSequenceDB\/}
   python fasta_remapper.py -i $file -m workspace/${file2}.names.txt.id2plate.txt.map.txt -o workspace;
done
