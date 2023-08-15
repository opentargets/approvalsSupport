## Genetic Evidence for FDA-Approved Drugs Pipeline

This repository contains the pipeline for generation and analysis of genetic evidence support for FDA-approved drugs from 2013-2022, as taken from Mullard publications in Nature Reviews Drug Discovery.

![Copy of Copy of Copy of Drugs2genevidence (3)](https://github.com/opentargets/approvalsSupport/assets/122813791/809b0ba6-c689-43e4-8cf6-00b5292d9efd)



### Data Structure

- **data/2013-2022_approvals_in.csv**: Input data with FDA-approved drugs. The file provides details on the drugs mapped to ChEMBL IDs, indications manually mapped to ontology terms (EFOs, MONDo, or ORDO), and the classification of diseases.

- **data/amendMoas**: A manually curated list of mechanisms of action absent in ChEMBL.
  
- **data/amendPhenotypes**: A manually curated list of related indications.

- **data/datasourceMetadata**: List of OpenTargets data sources.

### Processing Steps

1. **code/GE_search.R**:
   - Uses Google Cloud Cluster gs://open-targets-data-releases/ via Spark for processing.
   - Uses `data/2013-2022_approvals_in.csv` as the main input file and the above-mentioned lists.
   - Generates two output files:
     - **results/2013-2022_approvals_GE_src.csv**: Organizes evidence by data source of genetic evidence support and evidence type.
     - **results/2013-2022_approvals_GE_out.csv**: Extends the input file with additional columns for target IDs, interactor IDs, and related ontology terms.

3. **code/GE_types_add.R**: 
   - Processes the output from `GE_search.R`.
   - Generates **results/2013-2022_approvals_GE.csv** with new columns indicating the type of genetic evidence found for each drug-disease pair.

4. **code/GE_prior.py**:
   - Uses Google Cloud Cluster gs://open-targets-data-releases/ via Spark for processing.
   - Uses `results/2013-2022_approvals_GE.csv` to find the date for genetic evidence support.
   - Outputs **results/2013-2022_approvals_preGE0.csv**.

6. **Manual Curation**:
   - Complex cases in the above file are manually curated to produce **results/2013-2022_approvals_preGE.csv**.

### Analysis

1. **GE_Year_plot.R**:
   - Generates a plot showing genetic evidence support for approved drugs by year.

3. **OR_plot.R**: 
   - Computes and plots odds ratios (OR) for approvals with expedited review status and for those addressing serious conditions.

4. **Venn_plot.ipynb**:
   - Produces an intersection diagram, showing the overlap between approvals with genetic evidence support, approvals with expedited review status, and approvals for serious conditions.

### How to Use

1. Ensure you have the required input data and lists.
2. Run the processing steps in the order mentioned above.
3. Conduct manual curation where necessary.
4. Run the analysis scripts to generate visualizations and insights.

### Contributing

If you would like to contribute to this project or have any queries, please open an issue or submit a pull request.

### License

[MIT License](LICENSE) 

### Acknowledgments

Thanks to Mullard publications for providing the initial dataset of FDA-approved drugs. Special thanks to contributors who helped with manual curation and data improvements. Additionally, we'd like to extend our gratitude to OpenAI for ChatGPT-4, which assisted in the disease classification process and the generation of these instructions.
