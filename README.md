# soa - Systematic Operon Analysis

### Overview
A pipeline for systematic operon analysis (soa) in aims to elucidate genome wide transcription regulatory networks by using a comparative genomics approach to discover motifs that are implicated in the regulation of one or more operon(s). The distance between two operons in "regulon space" is determined as the pairwise distance between their regulatory motifs. A graphical network where the operons are nodes and the distance between any two nodes (i.e. operons) is defined by the distance of their motifs provides insight to the global regulatory network.

This core comparative genomics of this pipeline is a modified, high throughput, version of the pipeline designed for operon_conserve.

### Pipeline
<< ADD DIAGRAM OF PIPELINE >>

### Module/Class Definitions

* **GenomeFeature (soa_features.py)**
  *  A GenomeFeature object refers to a specific feature that is pulled from the Genebank record of interest.
  *  Memeber variables:
     * `genome_accession` - `str` 
     * `protein_accession` - `str`
     * `locus_tag` - `str`
     * `genome_fragment_name` - `str`
     * `coding_start` - `int`
     * `coding_end` - `int`
     * `strand` - `+ or -`
     * `five_end` - `int`
     * `three_end` - `int`
     * `aa_sequence` - `str`
     * `req_limit` - `int`
     * `sleep_time` - `int`
   * Member Functions:
     * `get_intergenic_distance(self, other)` - returns `int`
     * `__str__(self)`
     * `__eq__(self, other)`


* **AnnotatedHit(GenomeFeature) (soa_features.py)**
  * A child class of the GenomeFeature class. An instance of the AnnotatedHit class is made for each BLAST hit that is returned. In addition to all the information in the GenomeFeature class, the AnnotatedHit class contains additional information pertaining to the BLAST hit.
  * Additional Member Variables:
    * `query_accession` - `str`
    * `operon_id` - `str`
    * `align_start` - `int`
    * `align_end` - `int`
    * `feature_found` - `bool`
    * `percent_identity` - `float`
    * `alignment_seq` - `str`
  * Additional Functions:
    * `fetch_feature(self, record, coverage_cutoff=0.25)`
    * `__str__(self)`
    * `__eq__(self, other)`

* **Operon (soa_operon.py)**
  * An instance of the Operon class is made for each operon that is reassembeled/discovered in each GenomeFragment.
  * Member variables:
    * `features` - `GenomeFeature[]`
    * `operon_id` - `str`
    * `cluster_id` - `str`
    * `genome_fragment_name` - `str`
    * `genome_fragment_accession` - `str`
    * `genome_features` - `GenomeFeature[]` 
    * `strand` - `+ or -`
    * `promoter` - `str`
   * Member Functions:
     * `add_feature(self, feature)` 

* **GenomeFragment (soa_genome_fragment.py)**
  * An instance of the GenomeFragment class is made for each nucleotide record that is analyzed.
  * Member variables: 
    * `cache_directory` - `str`
    * `hits` - `AnnotatedHits[]`
    * `sorted__hits` - `AnnotatedHits{operon_id}`
    * `all_features` - `GenomeFeature[]`
    * `operons` - `Operon[]`
    * `name` - `str`
    * `genome_accession` - `str`
    * `assembly_accession` - `str`
    * `req_limit` - `int`
    * `sleep_time` - `int`
    * `full_record` - `[GenomeFeature]`
   * Relevant Member Functions:
     * `fetch_features(self)`
     * `fetch_record(self)`
     * `fetch_hit_features(self)`
     * `purge_hits(self)`
     * `add_hit(self, a_hit)`
     * `sort_hits(self)`
     * `assemble_operons(self, operon_id, feature_limit=3, intergenic_limit=1500)`
     * `get_promoters(self)`
 
* **soa_motif_analysis.py**
  * This module contains the functions that are used for the motif analysis part of the pipeline.
  * Relevant Member Functions:
    * `run_meme(input_file, output_dir, num_motifs, width_min, width_max, mode, pal, meme_exec_path)`
    * ` get_motifs(meme_data_dir, e_val_threshold)`
    * `def calculate_motif_distance(motif, other, distance_function, offset, padded, add_psuedocounts, psuedocount_val, scaling_factor)`
    * `def find_pattern(target_motif, self_score_ratio_threshold, kmer_pair_score_ratio_threshold, spacer_score_ratio_threshold)`


