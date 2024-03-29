# UBLRPSUI
This repo contains MATLAB code and dataset of the article

- Md Abul Hayat, Jingxian Wu, and Yingli Cao, "[Unsupervised Bayesian learning for rice panicle segmentation with UAV images](https://plantmethods.biomedcentral.com/articles/10.1186/s13007-020-00567-8)", Plant Methods, 2020.

- Please contact Dr. Jingxian Wu (`wuj@uark.edu`) before using this dataset for any research purpose.

## Dataset 
This dataset contains 12 images taken with Unmanned Areal Vehicle (UAV) at the Super Rice achievement Transformative Base (SRTB) of the Shenyang Agricultural University (SYAU) in northeastern China during 2017 and 2018.  Filename of Image #1 starts with `I1` and so forth for the following images. The filenames of the original images are ended with `_UAV.jpg`. The rice panicle in each image was manually labeled, and the corresponding manually labeled results of these 12 images have filenames ending with `_MAN.jpg`. This dataset can also be downloaded from [here](https://wuj.hosted.uark.edu/research/datasets/panicle/UBLRPSUI.zip). Information of these 12 images in this dataset is listed as follows.

| Image |	Distance | Resolution |
| --- | --- | --- |
| I1	| 3m | 820 × 865 |
| I2 | 3m | 810 × 800 |
| I3 | 3m | 1050 × 1075 |  
| I4 | 3m | 1050 × 1050 |  
| I5 | 3m | 1080 × 1040 |  
| I6 | 3m | 1120 × 1070 | 
| I7 | 6m | 415 × 410 | 
| I8 | 6m | 440 × 415 | 
| I9 | 6m | 430 × 430 | 
| I10 | 6m | 445 x 440 |
| I11 | 6m | 430 × 430 |   
| I12 | 6m | 430 × 420 |  

## MATLAB Code
Run `main_file.m` for executing the algorithm. Please contact Dr. Jingxian Wu (`wuj@uark.edu`) if you have any questions or suggestions regarding the code or dataset.

## Cite Paper
```
@article{hayat2020unsupervised,
  title={Unsupervised Bayesian learning for rice panicle segmentation with UAV images},
  author={Hayat, Md Abul and Wu, Jingxian and Cao, Yingli},
  journal={Plant methods},
  volume={16},
  number={1},
  pages={1--13},
  year={2020},
  publisher={BioMed Central}
}
```
