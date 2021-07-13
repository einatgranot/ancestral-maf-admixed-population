cwlVersion: v1.2
class: Workflow
label: Ancestry MAF global proportion wf
doc: |-
  This workflow calculates the ancestry-specific allele frequencies of bi-allelic genetic variants in admixed populations, based on global proportion ancestries.   

  Author: Einat Granot-Hershkovitz
$namespaces:
  sbg: https://sevenbridges.com

requirements:
- class: ScatterFeatureRequirement
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement

inputs:
- id: plist
  type: File
  sbg:fileTypes: CSV
  sbg:x: 0
  sbg:y: 107.12501525878906
- id: chromosome_gds
  type: File
  sbg:x: 0
  sbg:y: 321.3750305175781
- id: global_ancestry
  label: Global  ancestry file
  doc: Proportion of global ancestries SOL
  type: File?
  sbg:x: 0
  sbg:y: 214.25003051757812
- id: sol_unrel
  doc: List of third_degree unrelated SOL individuals
  type: File?
  sbg:fileTypes: .txt
  sbg:x: 0
  sbg:y: 0
- id: annotation
  type: File?
  sbg:x: 0
  sbg:y: 428.5000305175781
- id: output_prefix
  type: string
  sbg:exposed: true
- id: cpu_multiplier
  type: int
  sbg:exposed: true

outputs:
- id: concatenated_csv
  type: File?
  outputSource:
  - csv_combiner/concatenated_csv
  sbg:x: 1053.5733642578125
  sbg:y: 214.25001525878906

steps:
- id: block_creator
  label: block_creator
  in:
  - id: chromosome_gds
    source: chromosome_gds
  - id: plist
    source: plist
  run: ancestry-global-proportion.cwl.steps/block_creator.cwl
  out:
  - id: output_csv_blocks
  sbg:x: 192.39059448242188
  sbg:y: 207.25003051757812
- id: csv_combiner
  label: csv_combiner
  in:
  - id: csv_input
    source:
    - ancestry_maf_estimate_global_proportion_1/ancestry_csv
  - id: output_prefix
    source: output_prefix
  run: ancestry-global-proportion.cwl.steps/csv_combiner.cwl
  out:
  - id: concatenated_csv
  sbg:x: 812.2764892578125
  sbg:y: 214.25003051757812
- id: ancestry_maf_estimate_global_proportion_1
  label: Ancestry_MAF_estimate_global_proportion
  in:
  - id: global_ancestry
    source: global_ancestry
  - id: chromosome_gds
    source: chromosome_gds
  - id: annotation
    source: annotation
  - id: plist
    source: block_creator/output_csv_blocks
  - id: myIndex
    default: 1
  - id: sol_unrel
    source: sol_unrel
  - id: cpu
    default: 8
  - id: cpu_multiplier
    source: cpu_multiplier
  scatter:
  - plist
  run: |-
    ancestry-global-proportion.cwl.steps/ancestry_maf_estimate_global_proportion_1.cwl
  out:
  - id: ancestry_csv
  - id: workspace
  sbg:x: 523.4396362304688
  sbg:y: 211.5694580078125

hints:
- class: sbg:maxNumberOfParallelInstances
  value: '50'
- class: sbg:AWSInstanceType
  value: c5.18xlarge;ebs-gp2;100
sbg:appVersion:
- v1.2
- v1.1
sbg:content_hash: a8166c651d0dcefd0fa7ed5e3b93411c17eb8f219b8688740bd07845e139351d7
sbg:contributors:
- dave
sbg:createdBy: dave
sbg:createdOn: 1626203589
sbg:id: dave/cwl-apps/ancestry-maf-wf/1
sbg:image_url:
sbg:latestRevision: 1
sbg:modifiedBy: dave
sbg:modifiedOn: 1626203692
sbg:original_source: |-
  https://api.sb.biodatacatalyst.nhlbi.nih.gov/v2/apps/dave/cwl-apps/ancestry-maf-wf/1/raw/
sbg:project: dave/cwl-apps
sbg:projectName: CWL-Apps
sbg:publisher: sbg
sbg:revision: 1
sbg:revisionNotes: removed outdated cost info
sbg:revisionsInfo:
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1626203589
  sbg:revision: 0
  sbg:revisionNotes: Copy of dave/build-ancestry-maf-admixed-population/ancestry-maf-wf/11
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1626203692
  sbg:revision: 1
  sbg:revisionNotes: removed outdated cost info
sbg:sbgMaintained: false
sbg:toolAuthor: Einat Granot-Hershkovitz
sbg:validationErrors: []
sbg:wrapperAuthor: Einat Granot-Hershkovitz
