# hupehsuchus

This folder contains the data and R scripts used in a study that re-tested the paper by Fang et al. (2023). The published data from the original study can reproduce the published plots but contains a large amount of unauthetic data points that do not fit the morphology of respective species and mislead PCA to give the results that support the point made by Fang et al. (2023).

Fang et al., (2023) 
https://doi.org/10.1186/s12862-023-02143-9

# Primer

<ol>
    <li> Download the R and Cetacea_3D_Principal folders and place them in a single working directory. </li>
    <li> Modify the first lines of the R scripts so that your working directory is specified there. </li>
    <li> To reproduce the analysis of Fang et al. (2023), run GPA_Fang_asis.R </li>
    <li> To plot the landmarks used by Fang et al. (2023), run Plot_FangEA.R </li>
    <li> To run our preferred analysis using the 15 landmarks defined in our paper, run GPA_15.R </li>
</ol>

# List of Directories

<ul>
    <li> R      <br> Contains R scripts. </li>
    <li> Cetacea_3D_Principal  <br> Contains landmark files from 3D models of cetaceans. </li>
    <ul>
        <li> @fcsv <br> Contains 15 landmarks newly sampled from 3D models</li>
        <li> @o_fcsv <br> Contains 9 landmarks as defined by Fang et al. (2023) newly sampled from 3D models. </li>
    </ul>
</ul>

# List of R Files

- GPA_9.R   <br> GPA + PCA analysis of 9 landmarks as defined by Fang et al. (2023), newly sampled from 3D models.
- GPA_15.R  <br> GPA + PCA analysis of 15 landmarks newly sampled from 3D models.
- GPA_Fang.R    <br> GPA + PCA analysis of the data from Fang et al. (2023), with ichthyosauromorph landmarks revised.
- GPA_Fang_asis.R   <br> GPA + PCA analysis of the data from Fang et al. (2023) as is.
- GPA_Fang_asis_wo4.R   <br> GPA + PCA analysis of the data from Fang et al. (2023) with four inappropriate taxa removed.
- GPA_Fang_asis_wo19.R  <br> GPA + PCA analysis of the data from Fang et al. (2023) with four inappropriate and 15 problematic taxa removed.
- GPA_Fang_cetacea _wo4.R   <br> GPA + PCA analysis of the cetacean part of the data from Fang et al. (2023) with four inappropriate taxa removed.
- GPA_Fang_cetacea _wo19.R  <br> GPA + PCA analysis of the cetacean part of the data from Fang et al. (2023) with four inappropriate and 15 problematic taxa removed.
- GPA_Fang_cetacea.R    <br> GPA + PCA analysis of the cetacean part of the data from Fang et al. (2023), with with ichthyosauromorph landmarks revised
- GPA_Fang_cetacea_asis.R <br> GPA + PCA analysis of the cetacean part of the data from Fang et al. (2023) as is.
- GPA_Fang_wo19.R <br> GPA + PCA analysis of the data from Fang et al. (2023), with ichthyosauromorph landmarks revised and  four inappropriate and 15 problematic taxa removed. 
- Plot_FangEA.R <br> 2D Plot the landmarks from Fang et al. (2023) as well as new ones.
- Utility.R <br> Script containing functions and constants used by other script files.


# List of CSV Files

- Cb_9Bo.csv    <br> Nine landmarks as defined by Fang et al. (2023) for *Chaohusauru brevifemoralis*, assuming the posterior end of the skull is the basioccipital. 
- Cb_9St.csv    <br> Nine landmarks as defined by Fang et al. (2023) for *Chaohusauru brevifemoralis*, assuming the posterior end of the skull is the supratemporal.
- Cb_15.csv <br> Fifteen landmarks for *Chaohusauru brevifemoralis*
- Fang_EA_2023.csv  <br> Data file from Fang et al. (2023) based on a file downloaded on 5/14/2024.
- Fang_EA_2023_Hupehsuchus_Filter.csv   <br> Data file from Fang et al. (2023) based on a file downloaded on 4/12/2024.
- Hn_9Bo.csv    <br> Nine landmarks as defined by Fang et al. (2023) for *Hupehsuchus nanchangensis*, assuming the posterior end of the skull is the basioccipital.
- Hn_9St.csv    <br> Nine landmarks as defined by Fang et al. (2023) for *Hupehsuchus nanchangensis*, assuming the posterior end of the skull is the supratemporal.
- Hn_15.csv <br> Fifteen landmarks for *Hupehsuchus nanchangensis*
