## TRAINING CLASSIFIER ##
# trim reads
qiime feature-classifier extract-reads \
  --i-sequences /datasets/silva_ref_files/silva-138-99-seqs.qza \
  --p-f-primer GTGYCAGCMGCCGCGGTAA \
  --p-r-primer GGACTACNVGGGTWTCTAAT \
  --p-trunc-len 150 \
  --o-reads ref-seqs-trimmed.qza

# match trimmed reads with known taxonomies
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ref-seqs-trimmed.qza \
  --i-reference-taxonomy /datasets/silva_ref_files/silva-138-99-tax.qza \
  --o-classifier classifier.qza

# apply classifier to dataset
# Use the trained classifier to assign taxonomy to your reads (rep-seqs.qza)
qiime feature-classifier classify-sklearn \
  --i-classifier classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

# copying visualization file to local computer 
# cd: local computer (Desktop/MICB_305/Group_Project)
scp root@10.34.36.91:/work/manuscript/taxa-bar-plots.qzv .
scp root@10.34.36.91:/work/manuscript/table.qzv .
