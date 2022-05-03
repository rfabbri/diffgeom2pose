We provide both C++ and Matlab code that estimates where a picture was taken,
giving a single pair of 3D local geometry structure (differential geometry) and
its projection in a 2D image. Two particular cases are common
1. A pair of 3D-2D of SIFT features (the 3D is just the SIFT center with
   reconstructed 3D SIFT orientation from two other images)
2. A pair of 3D-2D point correspondences, where the points belong to curves
   (e.g., 2D-3D edgels, corners, or junctions)

The algorithm, P2PT, supersedes the usual P3P by using only two features (each with a
direction). This also solves a local pose problem which can be integrated along
curves to provide a more global estimate or a matching algorithm. The run time of
the C++ code is in the same order as the best P3P implementation available (tens
of microseconds).

A journal version is available at: https://ieeexplore.ieee.org/document/9057738

Camera Pose Estimation Using Curve Differential Geometry, IEEE Transactions on
Pattern Analysis and Machine Intelligence - PAMI, 2020, Ricardo Fabbri, Peter J.
Giblin and Benjamin Kimia.

This is an improvement over the research code originally written for the paper:

R. Fabbri, P. J. Giblin, B. B. Kimia, "Camera Pose Estimation Using Curve
Differential Geometry", ECCV 2012, Firenze, Italy (Lecture Notes in Computer
Science)

PDF and bibtex available at: http://multiview-3d-drawing.sourceforge.net


This work was developed at Brown University, University of Liverpool and
Rio de Janeiro State University.


## Main Function

### Matlab

```
matlab/rf_pose_from_point_tangents_root_find_function_any.m
```

This function will return all possible (Rotation, Translation) solutions for a
given pair of 3D-2D point-tangents (oriented points). These can then be tested
within standard RANSAC by you to keep only the one with the most inliers.

A simple demo with random data can be found in 
```
matlab/demo.m
```

### C++

See `cmd/diffgeom2pose.cxx`

## Disclaimer

### C++ version
This code is meant for production and is currently being incorporated to
mainstream SfM pipelines such as OpenMVG and Colmap.

### Matlab version
The matlab code is reseach code writen in a "lab" language (Matlab) and, despite
extensive experiments from our PAMI'20 paper showing it is robust and reliable,
one shouldn't expect it to be ready for production. A hardened version of this
code is the C++ implementation available in the folder `p2pt/`. The implementation
produces the results in the paper, and may well be very useful, but it may still
be *hugely* improved (as we have been doing for the C++ version). Oh, and don't
forget: please cite this paper ;^)


## Contact

Please contact Ricardo Fabbri <rfabbri@iprj.uerj.br> or any of the other authors for requests.

This work was developed at Brown University, University of Liverpool and State University
of Rio de Janeiro.
