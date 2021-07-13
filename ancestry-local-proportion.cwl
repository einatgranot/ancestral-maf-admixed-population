cwlVersion: v1.2
class: Workflow
label: Ancestry MAF local proportion wf
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
  sbg:x: -20.45480728149414
  sbg:y: 130.83798217773438
- id: chromosome_gds
  type: File
  sbg:x: -15.356437683105469
  sbg:y: 487.8316650390625
- id: sol_unrel
  doc: List of third_degree unrelated SOL individuals
  type: File?
  sbg:fileTypes: .txt
  sbg:x: 201.62596130371094
  sbg:y: 109.57933044433594
- id: local_ancestry_european
  label: European local ancestry file
  doc: European local ancestry file SOL
  type: File?
  sbg:x: 0
  sbg:y: 213.59375
- id: local_ancestry_amerindian
  label: Amerindian local ancestry file
  doc: Amerindian local ancestry file SOL
  type: File?
  sbg:x: 0
  sbg:y: 320.390625
- id: local_ancestry_african
  label: African local ancestry file
  doc: African local ancestry file SOL
  type: File?
  sbg:x: 188.72276306152344
  sbg:y: 424.68316650390625
- id: chain_file
  label: chain file
  doc: chain file liftover from hg38 to hg19
  type: File?
  sbg:x: 0
  sbg:y: 747.5781860351562
- id: lai_annotation
  label: LAI annotation file
  type: File?
  sbg:x: 36.50493621826172
  sbg:y: 632.6138916015625
- id: custom_input
  label: RUN LOCAL
  type: boolean
  sbg:x: 369.5346374511719
  sbg:y: 826.376220703125
- id: cpu_multiplyer
  type: int
  sbg:exposed: true

outputs:
- id: concatenated_csv
  type: File?
  outputSource:
  - csv_combiner/concatenated_csv
  sbg:x: 1130.3643798828125
  sbg:y: 289.7134704589844

steps:
- id: block_creator
  label: block_creator
  in:
  - id: chromosome_gds
    source: chromosome_gds
  - id: plist
    source: plist
  run: ancestry-local-proportion.cwl.steps/block_creator.cwl
  out:
  - id: output_csv_blocks
  sbg:x: 280.3526611328125
  sbg:y: 240.0100555419922
- id: ancestry_maf_estimate_local_proportion
  label: Ancestry MAF estimate local proportion
  in:
  - id: local_ancestry_amerindian
    source: local_ancestry_amerindian
  - id: chromosome_gds
    source: chromosome_gds
  - id: lai_annotation
    source: lai_annotation
  - id: plist
    source: block_creator/output_csv_blocks
  - id: myIndex
    default: 1
  - id: sol_unrel
    source: sol_unrel
  - id: cpu
    default: 8
  - id: local_ancestry_african
    source: local_ancestry_african
  - id: local_ancestry_european
    source: local_ancestry_european
  - id: chain_file
    source: chain_file
  - id: cpu_multiplyer
    source: cpu_multiplyer
  - id: custom_input
    source: custom_input
  scatter:
  - plist
  run: ancestry-local-proportion.cwl.steps/ancestry_maf_estimate_local_proportion.cwl
  when: $(inputs.custom_input)
  out:
  - id: ancestry_csv
  - id: workspace
  sbg:x: 682.5792846679688
  sbg:y: 514.12451171875
- id: csv_combiner
  label: csv_combiner
  in:
  - id: csv_input
    source:
    - ancestry_maf_estimate_local_proportion/ancestry_csv
  - id: output_prefix
    default: local_chr_MAF
  run: ancestry-local-proportion.cwl.steps/csv_combiner.cwl
  out:
  - id: concatenated_csv
  sbg:x: 914.3644409179688
  sbg:y: 459.44232177734375

hints:
- class: sbg:AWSInstanceType
  value: c5.18xlarge;ebs-gp2;100
- class: sbg:maxNumberOfParallelInstances
  value: '50'
sbg:appVersion:
- v1.2
- v1.1
sbg:content_hash: a7b66ed52fffdcaa1e20c275dd1d406de4d194902ffd82c203b4c7487d933aa8d
sbg:contributors:
- dave
sbg:createdBy: dave
sbg:createdOn: 1626203577
sbg:id: dave/cwl-apps/ancestry-local-proportion/1
sbg:image_url:
sbg:latestRevision: 1
sbg:modifiedBy: dave
sbg:modifiedOn: 1626203642
sbg:original_source: |-
  https://api.sb.biodatacatalyst.nhlbi.nih.gov/v2/apps/dave/cwl-apps/ancestry-local-proportion/1/raw/
sbg:project: dave/cwl-apps
sbg:projectName: CWL-Apps
sbg:publisher: sbg
sbg:revision: 1
sbg:revisionNotes: removed cost info that is out dated
sbg:revisionsInfo:
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1626203577
  sbg:revision: 0
  sbg:revisionNotes: Copy of dave/build-ancestry-maf-admixed-population/ancestry-local-proportion/15
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1626203642
  sbg:revision: 1
  sbg:revisionNotes: removed cost info that is out dated
sbg:sbgMaintained: false
sbg:toolAuthor: Einat Granot-Hershkovitz
sbg:validationErrors: []
sbg:wrapperAuthor: Einat Granot-Hershkovitz
