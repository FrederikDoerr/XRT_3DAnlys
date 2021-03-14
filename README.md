# XRT_3DAnlys

<!-- About The Project -->
## About The Project
This repository acts as a collection of MATLAB methodologies for the routine analysis of 3D XRT image data. This repository was created to provide a basis for others interested in getting started working with 3D image datasets using MATLAB.

<!-- Structure-->
## Structure
* _Demo: contains specific demo files with application examples (e.g. region-growing, watershed segementation) and selected demo image data (XRT_Capsule.zip). 
* _Export: default export folder
* _FileExchange: dependencies, please review indvidual file packages (Readme). Some files were modified to improve integration.
* _ImgPrc_ParameterFiles: contains parameter file 'XRT_ParameterImport.xlsx'. Please amend to add your own image datasets. Sheet '3DAnlys': general settings and filepaths, Sheet 'ImgPrc_Param': image processing parameter.
* _SubRoutines: methods/functions for routine 2D/3D image processing and analysis.
* XRT_3DAnlys_MAIN.m: main Matlab file to run XRT image analysis routine

<!-- Getting Started-->
## Getting Started
There are selected demo examples defined in file '.\_ImgPrc_ParameterFiles\XRT_ParameterImport.xlsx'. Please unpack '.\_Demo\XRT_Capsule.zip' with image data. Run XRT_3DAnlys_MAIN.m

To add your own 3D image analysis:
- Open parameter file '.\_ImgPrc_ParameterFiles\XRT_ParameterImport.xlsx' and amend/add row according to your application: define a unique identifier (ExpShorthand), ImagePixelSize in um/px is used for size calculations, add a directory with a series of images (Input_FolderPath), optional directory with ROI images (Input_FolderPath_ROI), export directory (Output_FolderPath), image format (ImgFormat), analysis routine (true/false): VolumeAnlys_On, PorosAnlys3D_On or BinSegment3D_On
- Move row to the top of the sheet '3DAnlys'
- Define a unique ImgPrc_Param_REF which links to Sheet ImgPrc_Param
- Parameters in ImgPrc_Param are used to process the image data (XRT_3DAnlys_VolumeBuild.m): (1) ImgPrc_imFilter e.g. 'localcontrast' or 'medfilt2' - see subroutine XRT_2DAnlys_ImageFilter.m, (2) ImgPrc_bw_thresMethod e.g. automatic 'ridlercalvard' or 'otsu' - included in XRT_3DAnlys_VolumeBuild.m, (3) binary noise reduction 2D/3D, (4) ROI method e.g. ShrinkWrap - see XRT_3DAnlys_Shrinkwrap.m
- Open the main file (XRT_3DAnlys_MAIN.m) and run script

The Matlab analysis will run through all datasets defined in sheet '3DAnlys' until it encounters an empty row.

<!-- References-->
## References (open access):
- Doerr, F. J. S., Florence, A. J. (2020). *A micro-XRT image analysis and machine learning methodology for the characterisation of multi-particulate capsule formulations.* International Journal of Pharmaceutics: X. [https://doi.org/10.1016/j.ijpx.2020.100041](https://doi.org/10.1016/j.ijpx.2020.100041)

Data repository: [https://doi.org/10.15129/e5d22969-77d4-46a8-83b8-818b50d8ff45](https://doi.org/10.15129/e5d22969-77d4-46a8-83b8-818b50d8ff45)  
Video Abstract: [https://strathprints.strath.ac.uk/id/eprint/71463](https://strathprints.strath.ac.uk/id/eprint/71463)  
Slide Deck: [https://doi.org/10.13140/RG.2.2.26289.20322](https://doi.org/10.13140/RG.2.2.26289.20322)  

- Doerr, F. J. S., Oswald, I. D. H., Florence, A. J. (2018). *Quantitative investigation of particle formation of a model pharmaceutical formulation using single droplet evaporation experiments and X-ray tomography.* Advanced Powder Technology. [https://doi.org/10.1016/j.apt.2018.09.027](https://doi.org/10.1016/j.apt.2018.09.027)

Pre-print: [https://strathprints.strath.ac.uk/66187](https://strathprints.strath.ac.uk/66187)  
Video Abstract: [https://strathprints.strath.ac.uk/71435/](https://strathprints.strath.ac.uk/71435/)  

- Doerr, F., Burns, L., Lee, B., Hinds, J., Davis-Harrison, R., Frank, S., & Florence, A. (2020). *Peptide isolation via spray drying : particle formation, process design and implementation for the production of spray dried glucagon.* Pharmaceutical Research, 1–19. [https://doi.org/10.1007/s11095-020-02942-5](https://doi.org/10.1007/s11095-020-02942-5)

Data repository: [https://doi.org/10.15129/dcb859db-fe0d-4a56-b001-3f3d7ac6c44a](https://doi.org/10.15129/dcb859db-fe0d-4a56-b001-3f3d7ac6c44a) 

