# Burton lab process for PCR of 16S amplicons (Illumina)
###### Kait Al, PhD. 2024


## Overview
The protocol detailed here is designed to amplify prokaryotes (bacteria and archaea) using paired-end 16S community sequencing on the Illumina platform. Primers 515F–806R target the V4 region, and 341F-805R target the V3-V4 region of the 16S SSU rRNA, respectively. See the primer index preparation document [here](1_Prepare_primer_index_plates.md) for more information.
 
This protocol is based on the Earth Microbiome Project [protocols](https://earthmicrobiome.org/protocols-and-standards/16s/) and was originally developed by Dr. Greg Gloor, see [here](https://github.com/ggloor/miseq_bin/blob/master/Illumina_SOP.pdf).

## Before you start: 
* sterilize the BSC and liquid handling robot enclosure with UV/HEPA & RNAzap
* use a fresh lab coat and gloves
* use pipettes designated for only PCR and PCR-grade filter tips
* primers (stock tubes and arrayed plates) are stored in the -80 °C freezer (chest)
* Promega GoTaq hot-start colorless master mix is stored in the -20 °C chest freezer (#7) in 1 mL aliquots (labelled 'MM'). 2 x 1 mL aliquots are required for a full 96-well plate.


## Procedure

1. If using pre-arrayed V4 plates generated as per [this protocol](1_Prepare_primer_index_plates.md), thaw the plate at 22 °C and briefly centrifuge. 

Note: If using V3-V4 primers or others that are not pre-arrayed in plates, thaw frozen primers if required, briefly centrifuge, and add 10 uL of each forward and reverse primer per PCR reaction. Attempt to evenly distribute the samples among as many L and R barcodes as possible. For example, for amplification of 24 samples to be sequenced, you might use all 24 F barcodes and a single R, or all 24 R. You would not use 4 different F with 6 different R barcodes, as you are more likely to see primer-dependent biases present in downstream analyses.

2. If using DNA that has been stored frozen, thaw at 22 °C and briefly centrifuge. Similarly thaw and briefly centrifuge the master mix.

3. Remove the rubber strip caps from DNA plate inside the biosafety cabinet. 

To remove the rubber caps, use forces sterilized with RNase AWAY and gently peel the strip off, being careful not to contact the underside of the strip. If a strip is touched by anything other than the sterile forceps or work surface, it should be discarded and a new strip should be used to replace it. There are enough extra strips for the occasional replacement, but not enough for new strips to be used for every column each time the DNA must be accessed.

*If using a liquid handling robot*

* UV-sterilize a surface area for all 12 strips to be placed at once. The area for each strip must be labelled 1-12 for its corresponding column, and spaced enough so that curved strips will not touch). A large sheet of tin foil that has been labelled with a sharpie and UV sterilized in the BSC works well for this.
* Once all strip caps are removed, seal the DNA plate with a sealing mat or foil seal to transfer between the BSC and the robot enclosure. Once the plate is in position inside the enclosure, remove the seal.
    
*If pipetting manually*

 * Take off one rubber strip cap from the DNA plate at a time to minimize the potential for cross-contamination. Individual strips can be placed on a single sterile kimwipe - discard and replace the wipe between strips.
       
4. Transfer 2 μL of DNA from the DNA plate into the corresponding wells of the primer/PCR plate. The PCR negative control should not have any DNA template added. If using very low abundance samples, a PCR positive control can be included from previously (successfully) sequenced DNA samples or a [ZYMO standard](https://zymoresearch.eu/collections/zymobiomics-microbial-community-standards).

5. *If using the liquid handling robot* manually transfer 2000 μL of master mix into a new sterile reservoir inside the enclosure. Using the multichannel 20 μL head, transfer 20 μL of the master mix from the reservoir to each well of the primer/PCR plate. *If pipetting manually*, transfer 20 μL master mix from its 1.5 mL aliquot tube directly into the primer/PCR plate. Pipette up and down briefly to mix.

6. Seal the primer/PCR plate with a foil plate seal very tightly and briefly centrifuge to remove air bubbles within the wells.

7. Carry out amplification in the thermocycler with the lid temperature maintained at 104 °C and the cycle conditions as stated in Table 1 or 2 for V4 or V3-V4 primers, respectively.

Table 1. PCR cycling conditions for V4 primers

| Steps | Temperature | Time | Cycles |
|:---:|:---:|:---:|:---:|
| Warm-up | 95 °C | 4 min | 1 |
| Denaturing | 95 °C | 1 min | 25 cycles |
| Annealing | 52 °C | 1 min  |  |
| Extension | 72 °C | 1 min  |  |
| Final extension | 72 °C | 5 min | 1 |
| Hold | 4 °C | Forever |  |

Table 2. PCR cycling conditions for V3-V4 primers

| Steps | Temperature | Time | Cycles |
|:---:|:---:|:---:|:---:|
| Warm-up | 95 °C | 4 min | 1 |
| Denaturing | 95 °C | 1 min | 25 cycles |
| Annealing | 50 °C | 1 min  |  |
| Extension | 72 °C | 1 min  |  |
| Final extension | 72 °C | 5 min | 1 |
| Hold | 4 °C | Forever |  |

Note: For low abundance samples, a higher cycle number can/should be used (25-35). We normally use 25 cycles to reduce chimera formation and partial products. A test amplification should be conducted to ensure that plateau is reached with this number of cycles - 25 cycles is more than sufficient for stool. 

8. Aliquots of random samples should be run on agarose gels to ensure that the reactions proceeded as planned.  Expected band size for V4 is ~300–350 bp, and for V3-V4 is ~440-460. Low-biomass samples may yield faint or no visible bands; alternative methods such as a Bioanalyzer (at Robarts) could be used to verify presence of PCR product.

9. Bring samples to David at Robarts for quantification using PicoGreen, pooling at equimolar concentrations, PCR#2 (which adds the illumina pad) and sequencing on the miseq.

When sequencing low bacterial abundance samples like urine and ureteral stents with the described methodology, as many as 500 samples may be run at a time and still achieve read depth capable of thorough taxonomic profiling (>15 million reads per run). However, with specimens higher in diversity and bacterial abundance such as feces, higher per-sample read depth necessitates fewer samples per run. Exact read thresholds and sequencing depth will depend on the study question, environment being sampled, and sequencing technology.
