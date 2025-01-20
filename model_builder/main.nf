#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.genomes = "${launchDir}/data/genomes.csv"
params.data_dir = "${launchDir}/data"
params.medium= null
params.scale = 1
params.method = "gapseq"
params.threads = 12
params.gapseq_bad_score = 50
params.gapseq_good_score = 100
params.min_reactions = 100
params.simulate = false
params.memoteformat = "json"
params.db_name = "database"
params.taxversion = "gtdb220"
params.anaerobic = false
params.growth = 0.01


def helpMessage() {
    log.info"""
    ~~~ Diener Lab Metabolic Model Builder Workflow ~~~

    Usage:
    A run using all,default parameters can be started with:
    > nextflow run main.nf --resume

    A run with all parametrs set would look like:
    > nextflow run main.nf --data_dir=./data --method=gapseq --media=media.vsv

    General options:
      --data_dir [str]              The main data directory for the analysis (must contain `raw`).
      --method [str]                The algorithm to use. Either `carveme` or `gapseq`. `gapseq`
                                    requires docker or singularity.
    Growth Media:
      --medium                      A file containing growth media specification.
                                    `*.tsv` for CARVEME and `*.csv` for gapseq.
    """.stripIndent()
}

params.help = false

process init_carveme {
  cpus 1
  publishDir "${params.data_dir}", mode: "copy", overwrite: true

  output:
  path("db_stats.txt")

  script:
  """
  #!/usr/bin/env python

  import subprocess
  from os import path
  from carveme import config, project_dir
  from carveme.cli.carve import first_run_check

  diamond_db = project_dir + config.get('generated', 'diamond_db')[:-5] + ".dmnd"

  if __name__ == "__main__":
    if path.exists(diamond_db):
      subprocess.check_output(["rm", diamond_db])
    first_run_check()
    with open("db_stats.txt", "w") as out:
      res = subprocess.check_output(["diamond", "dbinfo", "-d", diamond_db])
      out.write(res.decode())
  """
}

process find_genes {
  cpus 1
  memory "2GB"
  time "1h"
  publishDir "${params.data_dir}/genes", mode: "copy", overwrite: true

  input:
  tuple val(id), val(domain), path(assembly)

  output:
  tuple val("${id}"), path("${id}.ffn"), path("${id}.faa")

  script:
  """
  prodigal -p single -i ${assembly} -o ${id}.gff -d ${id}.ffn \
           -a ${id}.faa
  """
}

process build_carveme {
  cpus 2
  memory "8GB"
  time "1h"
  publishDir "${params.data_dir}/carveme_models", mode: "copy", overwrite: true

  input:
  tuple val(id), path(genes_dna), path(genes_aa), path(db_info)

  output:
  tuple val("${id}"), path("${id}.xml.gz")

  script:
  if (params.medium)
    """
    CPX_PARAM_THREADS=${task.cpus} OMP_NUM_THREADS=${task.cpus} \
    carve ${genes_aa} -o ${id}.xml.gz --gapfill ${params.medium} \
          --diamond-args "-p ${task.cpus} --more-sensitive --top 10" \
          --fbc2 -v
    """
  else
    """
    CPX_PARAM_THREADS=${task.cpus} OMP_NUM_THREADS=${task.cpus} \
    carve ${genes_aa} -o ${id}.xml.gz --diamond-args "-p ${task.cpus}" \
      --fbc2 -v
    """
}

process build_gapseq {
  cpus 1
  memory {1.6.GB * task.attempt}
  time {10.h * task.attempt}
  maxRetries 1
  publishDir "${params.data_dir}/gapseq_draft"

  errorStrategy { task.attempt < 2 ? "retry" : "ignore" }

  input:
  tuple val(id), val(domain), path(assembly)

  output:
  tuple val("${id}"), path("${id}-draft.RDS"), path("${id}-rxnWeights.RDS"),
        path("${id}-rxnXgenes.RDS"), path("${id}-all-Pathways.tbl"),
        path("${id}-all-Reactions.tbl"), path("${id}-Transporter.tbl")

  script:
  """
    GSTMP=\$(mktemp -d -t gapseq_XXXXXXXXXX)
    trap "rm -rf \$GSTMP" EXIT

    gapseq find -O -p all -T \$GSTMP -K 1 -v 1 -t ${domain} -b ${params.gapseq_good_score} ${assembly} > ${id}.log || true
    grep "Running time:" ${id}.log || exit 1

    TMPDIR=\$GSTMP gapseq find-transport -K 1 -v 1 -b ${params.gapseq_good_score} ${assembly} > ${id}.log || true
    grep "Running time:" ${id}.log || exit 1

    (( \$(grep -c "good_blast" ${id}-all-Reactions.tbl) >= ${params.min_reactions} )) || exit 1
    (( \$(wc -l < ${id}-Transporter.tbl) >= 4 )) || exit 1

    gapseq draft \
      -r ${id}-all-Reactions.tbl \
      -t ${id}-Transporter.tbl \
      -c ${assembly} \
      -b ${domain == 'Archaea' ? 'Archaea' : 'auto'} \
      -u ${params.gapseq_good_score} \
      -l ${params.gapseq_bad_score} \
      -p ${id}-all-Pathways.tbl
  """
}

process gapfill_gapseq {
  cpus 1
  memory {3.GB * task.attempt}
  time {6.h * task.attempt}
  maxRetries 1
  errorStrategy { task.attempt < 2 ? "retry" : "ignore" }
  publishDir "${params.data_dir}/gapseq_models", mode: "copy", overwrite: true

  input:
  tuple val(id), path(draft), path(weights),
        path(rxnXgenes), path(pathways),
        path(rxns), path(transporters)

  output:
  tuple val(id), path("${id}.xml.gz")

  script:
  if (params.medium)
    """
    cp ${launchDir}/${params.medium} medium.csv
    gapseq fill -m ${draft} -n medium.csv -c ${weights} -b ${params.gapseq_bad_score} -g ${rxnXgenes} -k ${params.growth}
    gzip ${id}.xml
    """
  else
    """
    gapseq medium -m ${draft} -p ${pathways} -o medium.csv ${params.anaerobic ? "-c cpd00007:0" : ""}
    gapseq fill -m ${draft} -n medium.csv -c ${weights} -b ${params.gapseq_bad_score} -g ${rxnXgenes} -k ${params.growth}
    gzip ${id}.xml
    """
}

process merge_gapseq {
  cpus 1
  memory "6GB"
  publishDir "${params.data_dir}", mode: "copy", overwrite: true
  time "4h"

  input:
  path(files)

  output:
  tuple path("pathways.csv"), path("reactions.csv"), path("transporters.csv")

  script:
  """
  #!/usr/bin/env python3

  import pandas as pd
  import glob

  def read_gapseq(fi):
    print(f"Reading file {fi}...")
    df = pd.read_csv(fi, sep="\\t", comment="#")
    if "bitscore" in df.columns:
      df = df[df.bitscore > ${params.gapseq_bad_score}]
    if "Prediction" in df.columns:
      df = df[df.Prediction == "true"]
    id = fi.split("/")[-1].split("-")[0]
    df["sample_id"] = id
    return df


  files = {
    "pathways": glob.glob("*-all-Pathways.tbl"),
    "reactions": glob.glob("*-all-Reactions.tbl"),
    "transporters": glob.glob("*-Transporter.tbl")
  }
  for what in ["pathways", "transporters"]:
    tables = []
    for fi in files[what]:
      tables.append(read_gapseq(fi))
    df = pd.concat(tables)
    df.to_parquet(f"{what}.csv.gz", index=False)
    del tables
  """
}

process check_model {
  cpus 1
  memory "4GB"
  time "30 m"
  publishDir "${params.data_dir}/model_qualities", mode: "copy", overwrite: true

  input:
  tuple val(id), path(model)

  output:
  tuple val("${id}"), path("${id}.json.gz")

  script:
  if (params.memoteformat == "json")
    """
    CPX_PARAM_THREADS=${task.cpus} OMP_NUM_THREADS=${task.cpus} \
    memote run ${model} --filename ${id}.json || true
    gzip ${id}.json
    """
  else
    """
    CPX_PARAM_THREADS=${task.cpus} OMP_NUM_THREADS=${task.cpus} \
    memote report snapshot ${model} --filename ${id}.html
    gzip ${id}.html
    """
}

process fba {
  cpus 1
  memory "2GB"
  time "1h"

  input:
  tuple val(id), path(model)

  output:
  tuple val(id), path("${id}_exchanges.csv"), path("${id}_growth_rate.csv")

  script:
  """
  #!/usr/bin/env python3

  import cobra
  import pandas as pd
  from os import path

  model = cobra.io.read_sbml_model("${model}")

  sep = '${params.method == "gapseq" ? "," : "\\t"}'
  postfix = '${params.method == "gapseq" ? "_e0" : "_e"}'

  exids = [r.id for r in model.exchanges]
  if "${params.medium}" != "null":
    media_df = pd.read_csv("${params.medium}", sep=sep).rename(
      columns={"maxFlux": "flux", "compounds": "compound"})
    if "flux" not in media_df.columns:
      media_df["flux"] = 0.1
    if "reaction" not in media_df.columns:
      media_df["reaction"] = "EX_" + media_df["compound"] + postfix
    media_df.index = media_df.reaction
    model.medium = ${params.scale} * media_df.flux[media_df.index.isin(exids)]

  rate = pd.DataFrame({"id": "${id}", "growth_rate": model.slim_optimize()}, index=[0])
  sol = cobra.flux_analysis.pfba(model)
  ex_fluxes = sol.fluxes[
    sol.fluxes.index.isin(exids)
    & (sol.fluxes.abs() > model.tolerance)
  ]
  met_names = ex_fluxes.index.to_series().apply(
    lambda idx: model.metabolites.get_by_id(idx.replace("EX_", "")).name)
  exchanges = pd.DataFrame({
    "assembly": "${id}",
    "reaction": ex_fluxes.index,
    "flux": ex_fluxes,
    "description": met_names
  })
  rate.to_csv("${id}_growth_rate.csv", index=False)
  exchanges.to_csv("${id}_exchanges.csv", index=False)
  """
}

process summarize_fba {
  cpus 1
  memory "2GB"
  time "1h"
  publishDir "${params.data_dir}", mode: "copy", overwrite: true

  input:
  path(exchanges)
  path(rates)

  output:
  tuple path("exchanges.csv"), path("growth_rates.csv")

  script:
  """
  #!/usr/bin/env python3

  import pandas as pd
  import glob

  res = map(pd.read_csv, glob.glob("*_exchanges.csv"))
  exchanges = pd.concat(res)
  exchanges.to_csv("exchanges.csv", index=False)

  res = map(pd.read_csv, glob.glob("*_growth_rate.csv"))
  growth_rates = pd.concat(list(res))
  growth_rates.to_csv("growth_rates.csv", index=False)
  """
}

process model_db {
  cpus 8
  memory "8GB"
  time "3h"

  publishDir "${params.data_dir}", mode: "copy", overwrite: true

  input:
  path(models)
  each path(genomes)

  output:
  path("${params.db_name}_${params.taxversion}_strain_1.zip")

  script:
  """
  #!/usr/bin/env python3

  import pandas as pd
  from glob import glob
  from micom.workflows import build_database, workflow
  from cobra.io import read_sbml_model

  def gr(f):
    return (f, read_sbml_model(f).slim_optimize())

  if __name__ == "__main__":
    manifest = pd.read_csv("${genomes}")
    taxa = manifest.lineage.str.split(";", expand=True).iloc[:, 0:7]
    taxa.columns = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
    manifest = pd.concat([manifest, taxa], axis=1)
    manifest["strain"] = manifest.id.str.replace("_", " ")
    manifest["file"] = manifest.id + ".xml.gz"
    manifest = manifest[manifest.file.isin(glob("*.xml.gz"))]
    rates = workflow(
      gr,
      manifest.file,
      threads=${task.cpus}
    )
    rates = pd.Series(dict(rates))
    bad = rates[manifest.file].isna() | (rates[manifest.file] < 1e-6)
    bad = manifest[bad.values]
    if len(bad) > 0:
      print(f"The following {len(bad)} models can not grow: {', '.join(bad.id)}.")
    manifest = manifest[~manifest.file.isin(bad.file)]
    manifest["gapseq_growth_rate"] = rates[manifest.file].values

    db = build_database(
      manifest,
      "${params.db_name}_${params.taxversion}_strain_1.zip",
      rank="strain",
      threads=${params.threads}
    )
  """

}

workflow {
  // Show help message
  if (params.help) {
      helpMessage()
      exit 0
  }

  Channel
    .fromPath("${params.genomes}")
    .splitCsv(header: true)
    .map{row -> tuple(row.id, row.lineage ==~ "d__Archaea" ? "Archaea" : "Bacteria", "${params.data_dir}/raw/${row.assembly}")}
    .set{genomes}

  def models = null
  if (params.method == "carveme") {
    init_carveme()
    find_genes(genomes)
    build_carveme(find_genes.out.combine(init_carveme.out))
    models = build_carveme.out
  } else if (params.method == "gapseq") {
    build_gapseq(genomes)
    build_gapseq.out.collect{it -> tuple(it[4], it[5], it[6])} | merge_gapseq
    gapfill_gapseq(build_gapseq.out)
    models = gapfill_gapseq.out
  } else {
    error "Method must be either `carveme` or `gapseq`."
  }

  if (params.simulate) {
    models | fba
    exchanges = fba.out.map{it -> it[1]}.collect()
    rates = fba.out.map{it -> it[2]}.collect()
    summarize_fba(exchanges, rates)
  }

  check_model(models)
  models
    .map{it -> it[1]}
    .collect()
    .set{all_models}
  model_db(all_models, Channel.fromPath("${params.genomes}"))
}
