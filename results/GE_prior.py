#!/usr/bin/env python

__author__ = "Maria J. Falaguera"
__date__ = "04 Aug 2023"

"""
GE_prior.py: Get dated GE for drug approvals.
"""

# sumbit job to gcloud machine
"""
gcloud dataproc jobs submit pyspark GE_prior.py --cluster=cf-novelty --project=open-targets-eu-dev --region="europe-west1"
"""

import numpy as np
from pyspark.sql import functions as F
from pyspark.sql import SparkSession

# Paths
ot_platform_version = "23.02"

# Input approvals
approvals_file = (
    "gs://ot-team/cfalaguera/approvals/{}/2013-2022_approvals_GE_final.csv".format(
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
    "gs://ot-team/cfalaguera/approvals/{}/2013-2022_approvals_GE_prec.csv".format(
        ot_platform_version
    )
)

# Establish spark connection
spark = SparkSession.builder.getOrCreate()

# Get dated GE for drug approvals
if 1:
    # approvals
    approvals = spark.read.option("header", "true").csv(approvals_file).persist()
    approvals.show()

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
            # "datasourceId",
            F.col("year").alias("yearGE"),
        )
    ).persist()
    genetics.show()

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
                # F.col("datasourceId").alias("datasourceIdTargetDiseaseGE"),
                F.col("yearGE").alias("yearTargetDiseaseGE"),
            ),
            ["diseaseId", "targetId"],
            "left",
        )
        .join(
            genetics.select(
                F.col("targetId").alias("interactorId"),
                "diseaseId",
                F.col("yearGE").alias("yearInteractorDiseaseGE"),
            ),
            ["diseaseId", "interactorId"],
            "left",
        )
        .join(
            genetics.select(
                "targetId",
                F.col("diseaseId").alias("relatedId"),
                F.col("yearGE").alias("yearTargetRelatedGE"),
            ),
            ["targetId", "relatedId"],
            "left",
        )
        .join(
            genetics.select(
                F.col("targetId").alias("interactorId"),
                F.col("diseaseId").alias("relatedId"),
                F.col("yearGE").alias("yearInteractorRelatedGE"),
            ),
            ["interactorId", "relatedId"],
            "left",
        )
        # get earliest GE year for every drug approval
        .groupBy(
            "yearApproval",
            "drugId",
        )
        .agg(
            # all combinations: disease vs target, related disease vs target, disease vs interacting target, related disease vs interacting target
            F.min("yearTargetDiseaseGE").alias("yearTargetDiseaseGE"),
            F.min("yearInteractorDiseaseGE").alias("yearInteractorDiseaseGE"),
            F.min("yearTargetRelatedGE").alias("yearTargetRelatedGE"),
            F.min("yearInteractorRelatedGE").alias("yearInteractorRelatedGE"),
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
