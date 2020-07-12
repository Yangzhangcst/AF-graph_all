# AF-graph_all
Source code on MSRC dataset for the paper "Affinity Fusion Graph-based Framework for Natural Image Segmentation". The AF-Graph paper can be found [here](https://arxiv.org/abs/2006.13542). A list of papers and datasets about natural/color image segmentation can be found in our [repository](https://github.com/Yangzhangcst/Natural-color-image-segmentation).

AF-Graph is modified from [the offcial SAS implementation](http://www.ee.columbia.edu/ln/dvmm/SuperPixelSeg/dlform.htm).


### Requirements
The code requires the version of Matlab2018a, Ubuntu 16.04.


### Data
The original MSRC dataset can be downloaded from [here](https://www.microsoft.com/en-us/research/project/image-understanding/?from=http%3A%2F%2Fresearch.microsoft.com%2Fvision%2Fcambridge%2Frecognition%2F). We place the processed dataset in `database/MSRC/` folder.


### Demo
Run the demo `demo_AF_MSRC.m` in parallel on 4 workers (`parfor` in `KSC_graph`).


### Results of MSRC dataset
The detailed results can be found in `evaluation_MSRC.txt`.

### Other results

Coming soon!


### Citing
If you find this repository useful in your research, please consider citing:
```
@INPROCEEDINGS{AF-Graph,  
  author={Y. {Zhang} and M. {Liu} and J. {He} and F. {Pan} and Y. {Guo}},  
  booktitle={arXiv:2006.13542},   
  title={Affinity Fusion Graph-based Framework for Natural Image Segmentation},   
  year={2020}}
```
