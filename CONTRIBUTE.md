# Nextflow Pipeline Conventions

Objectives:

- Establish a standardized file structure for Nextflow projects
- Promote modularity and maintainability across workflows
- Ensure clear conventions for configuration and parameter handling

## 0. Project Structure

- Follow a standardized directory layout for organizing workflows, configuration files, and modules.

```
<workflow_name>/
├── bin/                         (*)
│   └── <custom_script>.py       (*)
├── conf/                        (*)
│   ├── ext.config               (*)
│   └── resources.config         (*)
├── assets/                      (*)
│   └── input_schema.json        (*)
├── modules/
│   ├── config/
│   │   ├── process/
│   │   │   └── labels.config
│   │   ├── pipeline_info.config
│   │   └── profiles.config
│   ├── process/
│   │   ├── <process_name_1>/
│   │   │   └── main.nf
│   │   └── <process_name_2>/
│   │       └── main.nf
│   └── subworkflow/
│       └── <subworkflow_name_1>/
│           └── main.nf
├── process/                     (*)
├── subworkflow/                 (*)
├── data/
│   ├── test/
│   │   ├── inputs/
│   │   │   └── get_test_data.sh (*)
│   │   └── test.yml
├── nextflow.config
├── nextflow_schema.json
└── README.md

  '(*)' : optional
```

- `./bin`: Contains custom scripts or executables. Content in this directory is automatically added to the `PATH` for all tasks. Avoid or use cautiously, as explained later.
- `./conf`: Contains pipeline-specific configuration files.
- `assets/`: Contains auxiliary files such as input schemas or sample sheets.
- `./modules/config` : Contains general configuration files common to all pipelines.
- `process` and `subworkflow`: Contains pipeline-specific process/subworkflow definitions.
- `modules/process` and `modules/subworkflow`: Contains modular process/subworkflow definitions.

## 1. Help Documentation and Parameters

- Use `nf-validation@1.1.3` for managing help documentation:
  - **Include** in `main.nf` (refer to the template).
  - Create `nextflow_schema.json` based on the parameters in the `params{}` section of `nextflow.config`.

**Example** in [nextflow_schema.json](./nextflow_schema.json)

- Define resource parameters in the `resources_options` section, matching `modules/config/process/labels.config`.
- Manage tool parameterization logic in `./conf/ext.conf`.

**Example** in [.conf/ext.conf](./conf/ext.conf)

- **Note**: The plugin `nf-schema` is not compatible with EPI2ME because of its updated JSON schema.

## 2. Custom scripts and templates

- **Avoid** using `bin` in modules, as it needs to be defined at the pipeline level. 
- **Avoid** using modules binaries unless absolutely required, since it necessitates setting `nextflow.enable.moduleBinaries = true` in the main pipeline. See [Nextflow documentation on module binaries](https://www.nextflow.io/docs/latest/module.html#module-binaries)
- **Prefer** using templating in `templates`  which should be located alongside the module script `main.nf`.

**Example:**
```
  script:
  template "split_fasta_per_sequence.py"
```

## 3. Configuration Files

- Include multiple configuration files in nextflow.config. Load in the following order (from lowest to highest priority):
  - `modules/config/process/labels.config` - manages resources assigned to processes/labels.
  - `modules/config/profiles.config` - defines Nextflow profiles.
  - `modules/config/pipeline_info.config` - manage the nexflow engine reports config.
  - `conf/ext.conf` - manages `ext` passed to processes and controls tool arguments.
  - `conf/resources.conf` - manage pipeline specific resources
  
**NOTE:** Some configuration settings can be overridden by the `-c` option in Nextflow.

**Example** in `nextflow.config`:
```
includeConfig "modules/config/process/labels.config"
includeConfig "modules/config/profiles.config"
includeConfig "modules/config/pipeline_info.config"
includeConfig "conf/ext.config"
includeConfig "conf/ressources.config"
```

## 4. Modules

- Use modules to avoid duplication and maintain compatibility during updates:
- Link modules in the modules folder:

```sh
mkdir modules
git submodule add https://github.com/nexomis/nf-subworkflows.git modules/subworkflows
git submodule add https://github.com/nexomis/nf-config.git modules/config
git submodule add https://github.com/nexomis/nf-process.git modules/process
```

**NOTE:** Submodules are linked to a repository by their hash reference. To update them, you need to pull/push the changes again.

### 4.1 Process

- **Refer** to [nf-process guidelines](https://github.com/nexomis/nf-process/blob/main/README.md) rules and conventions regarding process.
- **Unique naming:** When including processes multiple times, ensure they have unique names to prevent conflicts.

**Example:**
```
include { GZ as GZ1; GZ as GZ2; GZ as GZ3 } from '../../../process/gz/main.nf'
```

### 4.2  Subworkflow

- **Integrate** most of the logic into modular subworkflows rather than the end pipeline.
- **Input Channels:**
  - Prefer a single queue channel as input for subworkflows. Merge inputs into tuples if necessary.
  - Include `null` values to handle optional inputs.

**Example in pipeline:**
```
  trimmedInputs
  | map {
    it[0].args_spades = make_spades_args(it[0])
    return [it[0].id, it]
  }
  | join(k2Inputs, by: 0, remainder: true)
  | join(refGenomeInputs, by: 0, remainder: true)
  | set {inputsForViralAssembly}
  VIRAL_ASSEMBLY(inputsForViralAssembly)
```

**Example in subworkflow:**
```
  take:
  inputs // (id, [meta, reads], [meta, k2_index], [meta, inputRefgenome])

  main:
  inputReads = inputs.map { [it[0], it[1]] }
  inputK2Index = inputs.filter { it[2] } | map { [it[0], it[2]] }
  inputRefGenome = inputs.filter { it[2] } | map { [it[0], it[3]] }
```

- **Maximize the outputs Channels:** Ensure subworkflows emit as many useful output channels as possible.
- **Avoid Global Parameters:** Do not use global parameters (`params`) directly in processes or subworkflows (e.g., `params.skip_step`).
- **Use `meta` Attributes:** Implement optional steps using `meta` attributes.

**Example 1:**
```
  inputReadsFromK3
  | map {it[1]}
  | branch {
      spades: it[0].assembler == "spades"
      no_assembly: true
    }
  | set { inputForAssembly }

  SPADES(inputForAssembly.spades)

  SPADES.out.scaffolds
  | set { scaffolds }
```

**Example 2:**
```
  finalScaffolds
  | filter { it[0].realign == "yes" }
  | set { toRealignScaffolds }

  BOWTIE2_BUILD(toRealignScaffolds)
  BOWTIE2_BUILD.out.idx
  | map { [it[0].id, it[1]] }
  | set { bwtIdx }
```
- **Optional Processing Steps:** Use strategies like `join` followed by `concat` and `unique` to insert optional processing steps.
```
joinInputForK2i1 = inputReads.join(inputK2i1, by:0)
  KRAKEN2_HOST1(joinInputForK2i1.map { it[1] }, joinInputForK2i1.map { it[2] })
  KRAKEN2_HOST1.out.unclassified_reads_fastq
  | GZ1
  | map {[it[0].id, it]}
  | concat(inputReads)
  | unique { it[0] }
  | set {inputReadsFromK1}
```
## 5. Typography Conventions

- Channel Instances: `lowerCamelCase`
- Groovy Variables (Non-Channels): `snake_case`
- Groovy Functions: `lowerCamelCase`
- Processes or Sub-Workflows: `UPPER_SNAKE_CASE`
- Process Inputs and Outputs: `snake_case`
- Sub-Workflow Inputs and Outputs: `lowerCamelCase` (Like channel instances)

## 6. Input Files/Directories by Sample

- Use sample sheets to manage input files/directories associated with samples:
- Utilize **unique IDs** in dedicated sample sheets to reference input files/directories in the main sample sheet, preventing redundant processing.
- Place sample sheets in the `assets/` directory. 

**Example** in [assets/input_schema.json](./assets/input_schema.json)

## 7. Publishing outputs

- Use workflow `publish` section to publish outputs from channels
- Subworkflows shall emit channels to be published
- Configure output directory using `-output-dir` CLI option (defaults to 'results')
- Use `output` block to customize publish targets
- Avoid using publishDir directive in processes

**Example:**
```nextflow
// main.nf
workflow {
    main:
    PRIMARY_FROM_DIR(data)

    publish:
    PRIMARY_FROM_DIR.out.trimmed >> 'fastp'
    PRIMARY_FROM_DIR.out.fastqc_trim_html >> 'fastqc_trim' 
    PRIMARY_FROM_DIR.out.fastqc_raw_html >> 'fastqc_raw'
    PRIMARY_FROM_DIR.out.multiqc_html >> 'multiqc'
}

output {
    // Basic target configuration
    'fastp' {
        enabled params.save_fastp
        mode 'copy'  // instead of symlink
    }

    // Custom path within output directory 
    'fastqc_trim' {
        path 'qc/fastqc/trimmed'
    }

    // Dynamic path based on metadata
    'fastqc_raw' {
        path { meta, html -> "qc/fastqc/raw/${meta.id}" }
    }

    // Index file to preserve metadata
    'multiqc' {
        enabled params.save_multiqc
        mode params.multiqc_mode
        index {
            path 'index.csv'
            header true
        }
    }
}
```

**Available output directives:**
- `enabled`: Enable/disable publishing for this target
- `mode`: Publishing mode ('symlink', 'copy', 'move', 'link', 'rellink', 'copyNoFollow')
- `path`: Custom path within output directory (string or closure)
- `overwrite`: Allow overwriting existing files
- `index`: Create an index file of published values (CSV or JSON)
  - `path`: Index file path (required)
  - `header`: Use first record keys as column names (CSV only)
  - `sep`: Value separator character (CSV only)
  - `mapper`: Transform values before writing

## 8. Testing

- Provide tests that can be executed with:
```
nextflow run . -params-file data/test/test.yml
```

- Set up testing parameters in `data/test/test.yml`
- Provide testing inputs in `data/test/inputs/`
  - Optionally, include a script to download input test data (e.g., `get_test_data.sh`).
- Test outputs shall be found in `results` by default, or as specified by `-output-dir` CLI option
