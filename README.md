# About

A Python-based computer vision application for recognition and visualizing constellation lines in night sky images.
This project implements a famous plate-solving algorithm in a simplified way for noise-avoiding goals. Plate-solving algorithm itself based on planar triangle matching to detect similar stars and draw lines between them over input images. The output is an annotated input image where the identified constellation is highlighted.

## Dataset


The project was build on [HYG Database](https://codeberg.org/astronexus/hyg), which consists of around 120000 star over many unified star catalogues. 
For this project, the dataset was manually modified for its needs: only the brightest stars from each constellation were remained, filtered and used to build a new catalogue of unique planar triangles used for matching. The resulting catalog contains approximately 700 carefully selected stars. 

## Examples


In the example below, the program succesfully recognized and annotated a raw image from Stellarium of the Cancer Constellation:

| Before | After |
|--------|-------|
| ![Raw image](examples/image1.png) | ![Annotated image](examples/cnc.jpg) |

# How to use


First of all, make sure you have python 3.10+ installed (project tested on 3.12, but should work on 3.10+ also), then install required dependencies:

```bash
cd constellation-recognition
pip install -r requirements.txt
```

To process any night sky image, run:

```bash
python matcher.py night_sky_image.png
```

The result will be saved into `examples/results` folder with the name of the recognized constellation.

# Results


