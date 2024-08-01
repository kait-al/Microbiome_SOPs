# Burton lab process for preparing primer index plates
###### Kait Al, PhD. 2023


## Overview
This protocol will yield 5 sets of 6 x 96 well plates (576 unique index pairs), using a combination of 24 forward and reverse barcodes. See the final barcode array layout [here](../Amplicon_SOPs/V4_primers_GolayBarcodes_Plates.xlsx).

V4 barcoded primers were designed by Dr. Greg Gloor, see [here](https://github.com/ggloor/miseq_bin/blob/master/Illumina_SOP.pdf) for his original methodology.
Also see the Earth Microbiome Project 16S protocol [page](https://earthmicrobiome.org/protocols-and-standards/16s/).
This arraying protocol was originally developed by [Dr. Jordan Bisanz](https://github.com/jbisanz).


## Primers

Both forward and reverse primers contain a 12 nucleotide barcode. Therefore, paired reads from a single sample have a unique, differentiating, combination of forward and reverse barcodes.

Primer sequence without linker, pad, or barcode are as follows:
* V4 (Caporaso version, see [here](https://doi.org/10.1073/pnas.1000080107) and [here](https://doi.org/10.1038/ismej.2012.8)), 515F–806R, FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
* V3-V4 (Bakt version, see [here](https://doi.org/10.1038/ismej.2011.41)), 341F-805R, FWD:CCTACGGGNGGCWGCAG; REV:GACTACHVGGGTATCTAATCC
    * *These are used less frequently by us compared to the V4 primers and so are not typically arrayed in 96 well plates as per this protocol (step 3 onwards). Instead, they are stored in the -80 °C freezer in aliquots at concentrations of 200 and 3.2 μM for long term stock storage and ad hoc use, respectively. Tubes of 3.2 μM must be thawed, briefly centrifuged, and 10 uL of each forward and reverse added per PCR reaction. See V3-V4 specific notes in the [16S_Amplicon_PCR document](3_16S_amplicon_PCR).*


## Before you start: 
* sterilize the BSC and automated pipettor cabinet with UV/HEPA & RNAzap
* use a fresh lab coat and gloves
* use pipettes designated for only PCR and PCR-grade filter tips


## 1. Order primers

Use the [primer ordering sheet example](../Amplicon_SOPs/Golay_indexed_primers_ordering_template.xlsx). 

## 2. Prepare primer stocks

2.1 Briefly centrifuge before opening the tube of lyophilized primers. 

2.2 Resuspend at 200 μM (aka pMole/μL) in PCR-grade water. Calculate the volume to add by multiplying the number of nmoles of primer by 5. For example if there were 84.5 nmoles, add 422.5 μL water.

2.3 Give them enough time to re-suspend well (incubate at RT for 10 minutes or they could be stored at 4 °C overnight before proceeding).

2.4 Store at -80 °C. When thawing to use later, centrifuge briefly before opening the tubes.

## 3. Calculate the working concentrations and volumes

Each PCR reaction ultimately requires 10 μL of each L and R primer at 3.2 μM. There are 24 L and R nucleotide balanced barcoded Golay primers, therefore we have 576 combinations. Each primer will be used 24 times in one set of plates, and we will be prepping 5 x plate sets.

10 μL x 24 combinations per primer x 5 sets of plates = 1200 μL of each primer is required. Add 50 μL of extra for pipetting allowance. 

C<sub>1</sub> x V<sub>1</sub> = C<sub>2</sub> x V<sub>2</sub>
200 μM x V<sub>1</sub> = 3.2 μM x 1250 μL
V<sub>1</sub> = 20 μL
Therefore we add 20 μL of 200 μM primer stock to 1230 μL PCR-grade water.

## 4. Prepare 1.5 mL tubes of primers

4.1 Label the side/ front of DNase/RNase free Eppendorfs with L1-24 and R1-24.

4.2 Once the 1250 µL aliquots of 3.2 µM primers are prepared in the tubes, centrifuge briefly.

4.3 Cut the tops off the L1-24 Eppendorf tubes at the hinge with RNaseZap’ed/ UV’ed scissors.

4.4 Place L1-24 tubes in robot Eppendorf tube rack inside the biosafety cabinet.

4.5 Cover the open tubes with a plate seal or other sterile cover when transporting from BSC into the robot cabinet.

## 5. Prepare stock plates

5.1 Transfer the tube rack into the robot cabinet and place in the correct dock position.

5.2 Set out 6 x 96 well plates (Axygen) in the correct docking positions.

5.3. Dispense L primers using 'KA_L_Pri_Golay' Beckman protocol (takes approximately 2.5 hours). You can use open an open (sterile) tip box reserved for this step as this protocol only uses 24 tips.

5.4 Discard empty L1-24 Eppendorf tubes. Repeat steps 4.3-5.1 with R1-24.

5.5 Dispense R primers using "KA_R_Pri_Golay" Beckman protocol (takes approximately 30 minutes per plate). A fresh tip box is required for every plate. Order is very important! Never dispense R before L!

5.6 These protocols generate 1 “stock” plate of all 6 unique layouts, each well within the plate contains 100 µL (50 µL of each L and R unique primer combo). You can proceed directly to step 6.2, or these plates can be sealed with an aluminum plate cover and stored at -80 °C before aliquoting.

## 6. Aliquot working plates from stock

6.1 If stock plates are frozen, thaw and briefly centrifuge prior to removing aluminum seal.

6.2 Portion 20 µL from the “stock” plates into 4 duplicates so you end up with 5 x Plate 1- Plate 6 and each plate contains 20 µL (10 µL x  L, 10 µL x R) using protocol "Golay_primer_aliquots".

6.3 Label all plates, seal with aluminum plate seals, and store at -80 °C until use.
