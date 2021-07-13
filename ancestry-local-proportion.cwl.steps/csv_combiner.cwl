cwlVersion: v1.2
class: CommandLineTool
label: csv_combiner
$namespaces:
  sbg: https://sevenbridges.com

requirements:
- class: DockerRequirement
  dockerPull: |-
    images.sb.biodatacatalyst.nhlbi.nih.gov/einatgranot/ancestral_maf_admixed_population:20210211
- class: InitialWorkDirRequirement
  listing:
  - entryname: csv_combiner.R
    writable: false
    entry: |
      require(tidyverse)

      source("cwl_input.R")

      df = tibble(path = fs::dir_ls(glob = "*.csv")) %>% 
                       group_by(path) %>% 
                       do({
                         temp = .
                         csv_data = read_csv(temp$path)
                       }) %>%
                       ungroup() %>%
                       select(-path) %>%
                       write_csv(paste0(output_prefix, "_concat.csv"))
  - entryname: csv_input.json
    writable: false
    entry: "\n$(inputs.csv_input)\n"
  - entryname: cwl_input.R
    writable: false
    entry: |
      output_prefix = "$(inputs.csv_input[0].basename.split("_")[0] + "_" + inputs.output_prefix)"
- class: InlineJavascriptRequirement

inputs:
- id: csv_input
  type: File[]
  sbg:fileTypes: CSV
- id: output_prefix
  type: string

outputs:
- id: concatenated_csv
  type: File?
  outputBinding:
    glob: '*concat.csv'

baseCommand:
- Rscript csv_combiner.R

hints:
- class: sbg:SaveLogs
  value: '*.R'
id: dave/build-ancestry-maf-admixed-population/csv-combiner/1
sbg:appVersion:
- v1.2
sbg:content_hash: a1f16aa401072a6089807853859411d74b95bbad91fc5127005771193566f5749
sbg:contributors:
- dave
sbg:createdBy: dave
sbg:createdOn: 1620395854
sbg:id: dave/build-ancestry-maf-admixed-population/csv-combiner/1
sbg:image_url:
sbg:latestRevision: 1
sbg:modifiedBy: dave
sbg:modifiedOn: 1621513376
sbg:project: dave/build-ancestry-maf-admixed-population
sbg:projectName: 'BUILD: Ancestry MAF Admixed Population'
sbg:publisher: sbg
sbg:revision: 1
sbg:revisionNotes: removed inherit metadata
sbg:revisionsInfo:
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1620395854
  sbg:revision: 0
  sbg:revisionNotes: Copy of einatgranot/2020-freq-estimation-latino-ancestries/csv-combiner/14
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1621513376
  sbg:revision: 1
  sbg:revisionNotes: removed inherit metadata
sbg:sbgMaintained: false
sbg:validationErrors: []
