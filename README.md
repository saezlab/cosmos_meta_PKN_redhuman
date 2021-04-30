# COSMOS PKN generation script

## Introduction

This repository contains the scripts necessary to generate the COSMOS meta PKN, combining omnipath for PPIs, STITCHdb for allosteric regulations and redHuman metabolic reaction entwork for the metabolic part.

In short, run scripts in their numbered order to generate the cosmos_PKN.csv file. (scripts are in the script folder)

## 01: Metabolic reaction network (01alt_ocean_redhuman_to_df.R)

The first script requires the ocean package to be installed (https://github.com/saezlab/ocean). Indeed, the ocean package allows to quickly generate and process the redHuman metabolic reaction network (https://www.nature.com/articles/s41467-020-16549-2) in a sif format compatible with the other ressources (STITCHdb and omnipath).

This script also generates linkers to connect the metabolic enzyems of the reaction network with the omnipath PPI network.

## 02: Allosteric regulations (02alt_STITCH_to_SIF.R)

This scripts requires to download the expanded chemical link information file of stitch db (http://stitch.embl.de/cgi/download.pl?UserId=egbsSvpFLSc8&sessionId=q5Z8nf55HyjV&species_text=Homo+sapiens).

It first filters out textmining based interactions from the STITCH allosteric regulatio list. Then, it formats the identifiers of STITCH to be coherent with the other ressources.

## 03: omnipath PPI (03alt_clean_omnipath.R)

This scripts generates a signed directed PPI from all inteeractions available in the omnipathR package.

## 04: Joining the pieces together (04alt_FUSION.R)

This script imports the three previously generated parts of the cosmos PKN and combines them together in a single network. It also filters out from STITCHdb allosteric reactions that do not involve an genes present in omnipath or  the metabolic reaction network.
