## INITIAL STEPS: DENMULTIPLEX AND DENOISE ##
# cd: local computer (Desktop/MICB_305/Group_Project)
# import dataset from computer onto server 
scp Desktop/MICB_305/Group_Project/metadata_filtered.csv root@10.34.36.91:/work/manuscript

# cd: remote server (/work/manuscript)
# #importing manifest file as well as demultiplexing
qiime tools import \
  --type "SampleData[SequencesWithQuality]" \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path /datasets/project_2/depression/depression_manifest.tsv \
  --output-path ./demux_seqs.qza
 
# converting demultiplexed file to visualization format
qiime demux summarize \
  --i-data demux_seqs.qza \
  --o-visualization demux.qzv

# copying visualization file to local computer 
# cd: local computer (Desktop/MICB_305/Group_Project)
scp root@10.34.36.91:/work/manuscript/demux.qzv .

# denoising
# deterine ASVs
# cd: remote server (/work/manuscript)
qiime dada2 denoise-single \
  --i-demultiplexed-seqs demux_seqs.qza \
  --p-trim-left 0 \
  --p-trunc-len 151 \
  --o-representative-sequences rep-seqs.qza \
  --o-table table.qza \
  --o-denoising-stats stats.qza

# visualize DADA2 outputs 
qiime metadata tabulate \
  --m-input-file stats.qza \
  --o-visualization stats.qzv

qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file metadata_filtered.tsv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv
