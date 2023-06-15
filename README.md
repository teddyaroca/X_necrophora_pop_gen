# Population genetics of **_Xylaria necrophora_**, an emerging pathogen of soybean in the southern United States.

### Below is information regarding methods pertinent to the data in this repository. These methods are an excerpt from the manuscript submitted to .

**Title:**

# Emergence of semi and mostly clonal lineages of the soybean pathogen **_Xylaria necrophora_**.

Authors: Teddy Garcia-Aroca<sup>a\*</sup> Jonathan Richards<sup>b</sup>, Jeremy Brown<sup>b</sup>, Paul P. Price<sup>c</sup>, Cheryl P. Andam<sup>d</sup>, and Vinson P. Doyle<sup>b\*</sup>


<sup>a</sup> Department of Plant Pathology, University of Nebraska-Lincoln, Lincoln, NE 68503.

<sup>b</sup> Department of Plant Pathology and Crop Physiology, Louisiana State University, Baton Rouge, Louisiana 70803.

<sup>c</sup> LSU AgCenter, Macon Ridge Research Station, Winnsboro, Louisiana 71295.

<sup>d</sup> Department of Biological Sciences, University at Albany, State University of New York, Albany, NY 12222.


Teddy Garcia-Aroca: http://orcid.org/0000-0002-7567-4363

Paul “Trey” Price: http://orcid.org/0000-0002-1004-3616

Jonathan K. Richards: https://orcid.org/0000-0001-9342-3595

Jeremy Brown: https://orcid.org/0000-0002-1447-8633

Cheryl P. Andam: https://orcid.org/0000-0003-4428-0924

Vinson P. Doyle: https://orcid.org/0000-0002-2350-782X



| **CONTENTS**                                         |
| -----------------------------------------------------|
|												|
| 1. [DATASETS](#datasets)                        |
| &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[RStudio](#RStudio)                      |
| &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[Metadata](#Metadata)                      |
| &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[Input](#Input)                      |
| &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[Manuscript](#Manuscript)                      |
| &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[Analyses](#Scripts)                      |
| &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[Output](#Output)                      |
|																			|
| 2. [INSTRUCTIONS](#instructions)														|



# 1. DATASETS

## RStudio

This folder contains the RStudio code (markdown) and procedure for analyses of population genomics. From loading and filtering vcf files to plotting, and exporting supplementary tables.

## Metadata

In this directory, we provide all datasets and metadata per species/strain.

## Input

Contains raw data from analyses, including variant calling format (vcf) and alignment files.

## Manuscript

A single file with the latest draft version of the manuscript.

## Scripts

Contains scripts used to parse and summarise data from raw sequences (fastq files) to variant calling format files (vcf). 

## Output

Contains figures and tables presented in the manuscript submitted to .


# 2. INSTRUCTIONS 

In order to contribute, make changes, suggestions, or provide feedback to this repository, do the following:

1. [Fork](https://docs.github.com/en/get-started/quickstart/fork-a-repo) this repository from the top right.

<img width="479" alt="for_repo_image" src="https://docs.github.com/assets/cb-23088/images/help/repository/fork_button.png">

2. Then, clone this repository to your laptop/desktop computer.

`
git clone https://github.com/teddyaroca/fungalecology.git
`

3. Make changes, add, or edit current files and commit those changes to your copy of this 
repository:

| Command | Description |
| :--- | :------------------------------------- |
| `git status` | You should be able to see the changes you made (new files/folders) in red |
| `git add .` | To add all the changes to the current repository |
| `git commit -m "I changed/added x,y,z files"` | To commit those changes back to github |
| `git push` | To push the changes back to your forked repository |

Note: If you are having troubles when you try `git push`, follow [these](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token) instructions to add your personal token in the command line.

4. Once you have pushed your changes, [submit a pull request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests):

<img width="479" alt="fork_repository_example" src="https://docs.github.com/assets/cb-26570/images/help/pull_requests/pull-request-start-review-button.png">

Your changes will be reviewed by the website owner (teddyaroca) and accepted or denied depending on the accuracy/usefulness of the proposed changes.




