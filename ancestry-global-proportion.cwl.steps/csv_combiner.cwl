cwlVersion: v1.1
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
  expressionLib:
  - |2-

    var setMetadata = function(file, metadata) {
        if (!('metadata' in file)) {
            file['metadata'] = {}
        }
        for (var key in metadata) {
            file['metadata'][key] = metadata[key];
        }
        return file
    };
    var inheritMetadata = function(o1, o2) {
        var commonMetadata = {};
        if (!o2) {
            return o1;
        };
        if (!Array.isArray(o2)) {
            o2 = [o2]
        }
        for (var i = 0; i < o2.length; i++) {
            var example = o2[i]['metadata'];
            for (var key in example) {
                if (i == 0)
                    commonMetadata[key] = example[key];
                else {
                    if (!(commonMetadata[key] == example[key])) {
                        delete commonMetadata[key]
                    }
                }
            }
            for (var key in commonMetadata) {
                if (!(key in example)) {
                    delete commonMetadata[key]
                }
            }
        }
        if (!Array.isArray(o1)) {
            o1 = setMetadata(o1, commonMetadata)
            if (o1.secondaryFiles) {
                o1.secondaryFiles = inheritMetadata(o1.secondaryFiles, o2)
            }
        } else {
            for (var i = 0; i < o1.length; i++) {
                o1[i] = setMetadata(o1[i], commonMetadata)
                if (o1[i].secondaryFiles) {
                    o1[i].secondaryFiles = inheritMetadata(o1[i].secondaryFiles, o2)
                }
            }
        }
        return o1;
    };

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
    outputEval: $(inheritMetadata(self, inputs.csv_input))

baseCommand:
- Rscript csv_combiner.R

hints:
- class: sbg:SaveLogs
  value: '*.R'
id: einatgranot/2020-freq-estimation-latino-ancestries/csv-combiner/14
sbg:appVersion:
- v1.1
sbg:content_hash: a0fe6c704ade1a5667f3199009379ba36de5ae6f9e17748895c9c4c24a5ca7488
sbg:contributors:
- dave
sbg:createdBy: dave
sbg:createdOn: 1609944735
sbg:id: einatgranot/2020-freq-estimation-latino-ancestries/csv-combiner/14
sbg:image_url:
sbg:latestRevision: 14
sbg:modifiedBy: dave
sbg:modifiedOn: 1613068765
sbg:project: einatgranot/2020-freq-estimation-latino-ancestries
sbg:projectName: 2020_Freq_estimation_Latino_ancestries
sbg:publisher: sbg
sbg:revision: 14
sbg:revisionNotes: ''
sbg:revisionsInfo:
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1609944735
  sbg:revision: 0
  sbg:revisionNotes:
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1609944946
  sbg:revision: 1
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1609945059
  sbg:revision: 2
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1609945581
  sbg:revision: 3
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1609945662
  sbg:revision: 4
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1609945845
  sbg:revision: 5
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1609946771
  sbg:revision: 6
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1609947200
  sbg:revision: 7
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1609947506
  sbg:revision: 8
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1609948272
  sbg:revision: 9
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1609948794
  sbg:revision: 10
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1610029616
  sbg:revision: 11
  sbg:revisionNotes: back to *csv glob
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1610030132
  sbg:revision: 12
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1610569549
  sbg:revision: 13
  sbg:revisionNotes: added parsing of chromosome from input file name for output prefix
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1613068765
  sbg:revision: 14
  sbg:revisionNotes: ''
sbg:sbgMaintained: false
sbg:validationErrors: []
