#!/usr/bin/env cwl-runner

class: CommandLineTool
id: "PanCanQC"
label: "PanCanQC"
cwlVersion: v1.0
description: |
    A Docker container for the DKFZ Bias Filter.
    ```
    Usage:
    # fetch CWL
    dockstore tool cwl --entry quay.io/jwerner_dkfz/pancanqc:latest > Dockstore.cwl
    # make a runtime JSON template and edit it
    dockstore tool convert cwl2json --cwl Dockstore.cwl > Dockstore.json
    # run it locally with the Dockstore CLI
    dockstore tool launch --entry quay.io/jwerner_dkfz/pancanqc:latest --json Dockstore.json
    ```

dct:creator:
  foaf:name: Ivo Buchhalter
  foaf:mbox: "mailto:i.buchhalter@dkfz-heidelberg.de"

requirements:
  - class: DockerRequirement
    # dockerPull: "quay.io/jwerner_dkfz/pancanqc:latest"
    dockerPull: "wernerjo/pancanqc"
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - $(inputs.input_control)
      - $(inputs.input_control_index)
      - $(inputs.input_tumor)
      - $(inputs.input_tumor_index)

inputs:
  input_control:
    type: File
    default: "/home/pcawg/data/bams/control.bam"
    doc: "Absolute filename of control bam file"
    inputBinding:
      position: 1
      valueFrom: $(self.basename)
  input_control_index:
    type: File
    default: "/home/pcawg/data/bams/control.bam.bai"
    doc: "Absolute filename of control bam index file"
  input_tumor:
    type: File
    default: "/home/pcawg/data/bams/tumor.bam"
    doc: "Absolute filename of tumor bam file"
    inputBinding:
      position: 2
      valueFrom: $(self.basename)
  input_tumor_index:
    type: File
    default: "/home/pcawg/data/bams/tumor.bam.bai"
    doc: "Absolute filename of tumor bam index file"

outputs:
  output_results_folder:
    type: Directory
    outputBinding:
      glob: results
    doc: "Output results folder"

baseCommand: ["/home/pcawg/data/scripts/PCAWG_QC.sh"]
