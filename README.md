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
    * Functions:
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
   * Relevant Member functions:
     * `fetch_features(self)`
     * `fetch_record(self)`
     * `fetch_hit_features(self)`
     * `purge_hits(self)`
     * `add_hit(self, a_hit)`
     * `sort_hits(self)`
     * `assemble_operons(self, operon_id, feature_limit=3, intergenic_limit=1500)`
     * `get_promoters(self)`
     * 


