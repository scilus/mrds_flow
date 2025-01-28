# mrds_flow pipeline

This pipeline fits a multi-tensor model to diffusion-weighted imaging (DWI) data using the Multi-Resolution Discrete Search (MRDS) method. The Track Orientation Density Imaging (TODI) technique serves as a model selector, determining the optimal number of tensors to fit at each voxel.

Should you use this pipeline for your research, **please cite the following**

```
Dhollander, T., Emsell, L., Van Hecke, W., Maes, F., Sunaert, S., Suetens, P. (2014).
Track Orientation Density Imaging (TODI) and Track Orientation Distribution (TOD) based tractography.
NeuroImage, 94, 312-336. https://doi.org/10.1016/j.neuroimage.2013.12.047

Coronado-Leija, R., Ramirez-Manzanares, A., Marroquin, J.L. (2017).
Estimation of individual axon bundle properties by a Multi-Resolution Discrete-Search method.
Medical Image Analysis, 42, 26-43. https://doi.org/10.1016/j.media.2017.06.008

Hernandez-Gutierrez, E., Coronado-Leija, R., Ramirez-Manzanares, A., Barakovic, M., Magon, S., Descoteaux, M. (2023).
Improving Multi-Tensor Fitting with Global Information from Track Orientation Density Imaging.
Computational Diffusion MRI. CDMRI 2023. Lecture Notes in Computer Science, vol 14328. Springer, Cham. https://doi.org/10.1007/978-3-031-47292-3_4
```

Requirements
------------

- [Nextflow](https://www.nextflow.io)
- [scilpy](https://github.com/scilus/scilpy)

Singularity/Docker
-----------

If you are on Linux, we recommend using the Singularity to run mrds_flow pipeline.

If you have Apptainer (Singularity), build the .sif image with:
`singularity build mrds-flow_shuffle.sif docker://scilus/mrds-flow:shuffle`.

Then, launch your Nextflow command with:
`-with-singularity ABSOLUTE_PATH/mrds-flow_shuffle.sif`.

If you are on MacOS or Windows, we recommend using the Docker container to run mrds_flow pipeline.

Pull the image with:
`docker pull https://hub.docker.com/r/scilus/mrds-flow`.

Launch your Nextflow command with:
`-with-docker scilus/mrds-flow:shuffle`.

Usage
-----

See *USAGE* or run `nextflow run main.nf --help`