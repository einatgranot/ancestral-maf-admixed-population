cwlVersion: v1.1
class: CommandLineTool
label: block_creator
$namespaces:
  sbg: https://sevenbridges.com

requirements:
- class: ShellCommandRequirement
- class: DockerRequirement
  dockerPull: |-
    images.sb.biodatacatalyst.nhlbi.nih.gov/einatgranot/ancestral_maf_admixed_population:20210211
- class: InitialWorkDirRequirement
  listing:
  - entryname: block_creator.R
    writable: false
    entry: |+
      require(tidyverse)

      source("cwl_input.R")

      plist = read_csv(plist_path, col_names = c("chromosome", "job_num")) %>% 
          glimpse() %>%
        mutate(chromosome = as.integer(chromosome)) %>% 
        mutate(job_num = as.integer(job_num)) %>% 
        filter(chromosome == chromosome_number) %>% 
        glimpse() %>% 
        group_by(job_num) %>% 
        do({
          temp = .
          write_csv(temp, paste0("block_", temp$job_num, ".csv"))
          
          temp = temp
          
        }) %>%
        glimpse()

  - entryname: cwl_input.R
    writable: false
    entry: |-
      plist_path="$(inputs.plist.path)" 
      chromosome_number=$(inputs.chromosome_gds.path.split("chr")[1].split("_")[0])
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
- id: chromosome_gds
  type: File
  inputBinding:
    position: 0
    shellQuote: false
- id: plist
  type: File
  inputBinding:
    position: 0
    shellQuote: false
  sbg:fileTypes: CSV

outputs:
- id: output_csv_blocks
  type: File?
  outputBinding:
    glob: '*.csv'
    outputEval: $(inheritMetadata(self, inputs.chromosome_gds))
stdout: standard.out

baseCommand:
- Rscript block_creator.R

hints:
- class: sbg:SaveLogs
  value: '*.R'
- class: sbg:SaveLogs
  value: standard.out
id: einatgranot/2020-freq-estimation-latino-ancestries/block-creator/11
sbg:appVersion:
- v1.1
sbg:content_hash: a881204a62f9587bb2e659d406181b68de3e9b6c92fc1aab096c6a97f48e4ab34
sbg:contributors:
- dave
sbg:createdBy: dave
sbg:createdOn: 1609875429
sbg:id: einatgranot/2020-freq-estimation-latino-ancestries/block-creator/11
sbg:image_url:
sbg:latestRevision: 11
sbg:modifiedBy: dave
sbg:modifiedOn: 1613068663
sbg:project: einatgranot/2020-freq-estimation-latino-ancestries
sbg:projectName: 2020_Freq_estimation_Latino_ancestries
sbg:publisher: sbg
sbg:revision: 11
sbg:revisionNotes: |-
  images.sb.biodatacatalyst.nhlbi.nih.gov/einatgranot/ancestral_maf_admixed_population
sbg:revisionsInfo:
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1609875429
  sbg:revision: 0
  sbg:revisionNotes:
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1609875666
  sbg:revision: 1
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1609875806
  sbg:revision: 2
  sbg:revisionNotes: added plist input
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1609875879
  sbg:revision: 3
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1609876413
  sbg:revision: 4
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1609880389
  sbg:revision: 5
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1609880967
  sbg:revision: 6
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1609881624
  sbg:revision: 7
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1609882176
  sbg:revision: 8
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1610031313
  sbg:revision: 9
  sbg:revisionNotes: split on chr now
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1610031900
  sbg:revision: 10
  sbg:revisionNotes: ''
- sbg:modifiedBy: dave
  sbg:modifiedOn: 1613068663
  sbg:revision: 11
  sbg:revisionNotes: |-
    images.sb.biodatacatalyst.nhlbi.nih.gov/einatgranot/ancestral_maf_admixed_population
sbg:sbgMaintained: false
sbg:validationErrors: []
