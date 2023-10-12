#!/usr/bin/env python

__author__ = "Maria J. Falaguera (mariaf@ebi.ac.uk)"
__date__ = "12 Oct 2023"

"""
GE_prior.py: Get dated GE for drug approvals.
"""

# sumbit job to gcloud machine
"""
gcloud dataproc jobs submit pyspark GE_prior.py --cluster=cf-timeseries --project=open-targets-eu-dev --region="europe-west1"
"""

import numpy as np
from pyspark.sql import functions as F
from pyspark.sql import SparkSession
from pyspark.sql import types as T

# Paths
ot_platform_version = "23.02"

# Input approvals
approvals_file = (
    # "gs://ot-team/cfalaguera/approvals/{}/2013-2022_approvals_GE_final.csv".format(
    "gs://ot-team/cfalaguera/approvals/{}/2013-2022_approvals_date_v1.csv".format(
        ot_platform_version
    )
)

# Input obtained by (1) disease ontology propagation of direct evidence (to obtain indirect evidence)
# and (2) evidence dating using literature file in the Platform
evidenceIndirectDated_file = (
    "gs://ot-team/cfalaguera/novelty/{}/evidenceIndirectDated".format(
        ot_platform_version
    )
)

# Output
approvalsPrecedence_file = (
    "gs://ot-team/cfalaguera/approvals/{}/2013-2022_approvals_GE_prec".format(
        ot_platform_version
    )
)

# Establish spark connection
spark = SparkSession.builder.getOrCreate()

# Get dated GE for drug approvals
if 1:
    # approvals
    approvals = spark.read.option("header", "true").csv(approvals_file).persist()

    # genetics
    genetics = (
        spark.read.parquet(evidenceIndirectDated_file)
        # use only genetic sources
        .filter(
            ~F.col("datasourceId").isin(
                [
                    "chembl",
                    "expression_atlas",
                    "sysbio",
                    "europepmc",
                    "phenodigm",
                    "reactome",
                    "phewas_catalog",
                ]
            )
        )
        # GE dated
        .select(
            "targetId",
            "diseaseId",
            "pmid",
            F.col("year").cast(T.IntegerType()).alias("year"),
        )
    ).persist()

    # approvals + genetics
    results = (
        # approvals
        approvals.select(
            "yearApproval",
            "drugId",
            "diseaseIds",
            "targetIds",
            "relatedIds",
            "interactorIds",
        )
        # has GE
        .filter(F.col("hasAnyGE") == "TRUE")
        # explode targets and diseases
        .withColumn("targetId", F.explode_outer(F.split("targetIds", ",")))
        .drop("targetIds")
        .withColumn("interactorId", F.explode_outer(F.split("interactorIds", ",")))
        .drop("interactorIds")
        .withColumn("diseaseId", F.explode_outer(F.split("diseaseIds", ",")))
        .drop("diseaseIds")
        .withColumn("relatedId", F.explode_outer(F.split("relatedIds", ",")))
        .drop("relatedIds")
        # hasTargetDiseaseGE
        .join(
            genetics.select(
                "targetId",
                "diseaseId",
                F.col("year").alias("yearTargetDiseaseGE"),
                F.col("pmid").alias("pmidTargetDiseaseGE"),
            ),
            ["diseaseId", "targetId"],
            "left",
        )
        .join(
            genetics.select(
                F.col("targetId").alias("interactorId"),
                "diseaseId",
                F.col("year").alias("yearInteractorDiseaseGE"),
                F.col("pmid").alias("pmidInteractorDiseaseGE"),
            ),
            ["diseaseId", "interactorId"],
            "left",
        )
        .join(
            genetics.select(
                "targetId",
                F.col("diseaseId").alias("relatedId"),
                F.col("year").alias("yearTargetRelatedGE"),
                F.col("pmid").alias("pmidTargetRelatedGE"),
            ),
            ["targetId", "relatedId"],
            "left",
        )
        .join(
            genetics.select(
                F.col("targetId").alias("interactorId"),
                F.col("diseaseId").alias("relatedId"),
                F.col("year").alias("yearInteractorRelatedGE"),
                F.col("pmid").alias("pmidInteractorRelatedGE"),
            ),
            ["interactorId", "relatedId"],
            "left",
        )
        .na.fill(float("+inf"))  # line required for min function to work later
        # get earliest GE year for every drug approval
        .groupBy(
            "yearApproval",
            "drugId",
        )
        .agg(
            # all combinations: disease vs target, related disease vs target, disease vs interacting target, related disease vs interacting target
            F.min(F.struct("yearTargetDiseaseGE", "pmidTargetDiseaseGE")).alias(
                "yearTargetDiseaseGE,pmidTargetDiseaseGE"
            ),
            F.min(F.struct("yearInteractorDiseaseGE", "pmidInteractorDiseaseGE")).alias(
                "yearInteractorDiseaseGE,pmidInteractorDiseaseGE"
            ),
            F.min(F.struct("yearTargetRelatedGE", "pmidTargetRelatedGE")).alias(
                "yearTargetRelatedGE,pmidTargetRelatedGE"
            ),
            F.min(F.struct("yearInteractorRelatedGE", "pmidInteractorRelatedGE")).alias(
                "yearInteractorRelatedGE,pmidInteractorRelatedGE"
            ),
        )
        .select(
            "yearApproval",
            "drugId",
            # target-disease
            F.col("yearTargetDiseaseGE,pmidTargetDiseaseGE.yearTargetDiseaseGE").alias(
                "yearTargetDiseaseGE"
            ),
            F.col("yearTargetDiseaseGE,pmidTargetDiseaseGE.pmidTargetDiseaseGE").alias(
                "pmidTargetDiseaseGE"
            ),
            # interactor-disease
            F.col(
                "yearInteractorDiseaseGE,pmidInteractorDiseaseGE.yearInteractorDiseaseGE"
            ).alias("yearInteractorDiseaseGE"),
            F.col(
                "yearInteractorDiseaseGE,pmidInteractorDiseaseGE.pmidInteractorDiseaseGE"
            ).alias("pmidInteractorDiseaseGE"),
            # target-related
            F.col("yearTargetRelatedGE,pmidTargetRelatedGE.yearTargetRelatedGE").alias(
                "yearTargetRelatedGE"
            ),
            F.col("yearTargetRelatedGE,pmidTargetRelatedGE.pmidTargetRelatedGE").alias(
                "pmidTargetRelatedGE"
            ),
            # interactor-related
            F.col(
                "yearInteractorRelatedGE,pmidInteractorRelatedGE.yearInteractorRelatedGE"
            ).alias("yearInteractorRelatedGE"),
            F.col(
                "yearInteractorRelatedGE,pmidInteractorRelatedGE.pmidInteractorRelatedGE"
            ).alias("pmidInteractorRelatedGE"),
        )
        # replace infinite value in year* columns with NULL
        .replace(
            2147483647,
            None,
            subset=[
                "yearTargetDiseaseGE",
                "yearTargetRelatedGE",
                "yearInteractorDiseaseGE",
                "yearInteractorRelatedGE",
            ],
        )
        # join with approvals metadata
        .join(
            approvals,
            ["drugId", "yearApproval"],
            "right",
        )
    )

    results.write.parquet(approvalsPrecedence_file)
    print(approvalsPrecedence_file)


# Count approvals with prior GE
if 1:
    print(approvalsPrecedence_file)
    approvals = spark.read.parquet(approvalsPrecedence_file)

    print("approvals:")
    approvals.groupBy("yearApproval").agg(
        F.countDistinct("drugId").alias("nDrugId")
    ).orderBy("yearApproval").show()

    print("approvals with genetic evidence:")
    approvals.filter((F.col("hasAnyGE") == "TRUE")).select(
        "yearApproval", "drugId"
    ).groupBy("yearApproval").agg(F.countDistinct("drugId").alias("nDrugId")).orderBy(
        "yearApproval"
    ).show()

    print("approvals with preceding genetic evidence:")
    approvals.filter(
        (F.col("hasAnyGE") == "TRUE")
        & (
            (F.col("yearTargetDiseaseGE") < F.col("yearApproval"))
            | (F.col("yearInteractorDiseaseGE") < F.col("yearApproval"))
            | (F.col("yearTargetRelatedGE") < F.col("yearApproval"))
            | (F.col("yearInteractorRelatedGE") < F.col("yearApproval"))
        )
    ).select("yearApproval", "drugId").groupBy("yearApproval").agg(
        F.countDistinct("drugId").alias("nDrugId")
    ).orderBy(
        "yearApproval"
    ).show()
