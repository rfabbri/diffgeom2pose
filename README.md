This code estimates where a picture was taken, giving a single pair of 3D local geometry
structure (differential geometry) and its projection in the 2D image. Two particular cases are common
1. A pair of 3D-2D of SIFT features (the 3D is just the SIFT center with reconstructed 3D SIFT orientation from two other images)
2. A pair of 3D-2D point correspondences, where the points belong to curves (e.g., 2D-3D edgels, corners, or junctions)

This solves a local pose problem which can be integrated along curves to provide a more global estimate or a matching algorithm.

This is an improvement over the research code originally written for the paper:

R. Fabbri, P. J. Giblin, B. B. Kimia, "Camera Pose Estimation Using Curve
Differential Geometry", ECCV 2012, Firenze, Italy (Lecture Notes in Computer
Science)

A journal version is under review.

This work was developed at Brown University, University of Liverpool and State University
of Rio de Janeiro.

PDF and bibtex available at: http://multiview-3d-drawing.sourceforge.net


## Main Function

```
rf_pose_from_point_tangents_root_find_function_any.m
```

This function will return all possible (Rotation, Translation) solutions for a
given pair of 3D-2D point-tangents (oriented points). These can then be tested
within standard RANSAC by you to keep only the one with the most inliers.

A simple demo with random data can be found in 
```
 demo.m
```

## Disclaimer

This is reseach code writen in a "lab" language (Matlab) and, as such, one
shouldn't expect it to be ready for production. The implementation produces the
results in the paper, and may well be very useful, but it may be *hugely*
improved. We're working on better algorithms, but this one works OK. Feel free
to develop and publish better algorithms yourself. But, above all, don't blame
us if you see any imperfections in this research code. Oh, and don't forget:
please cite this paper ;^)


## Contact

Please contact Ricardo Fabbri <rfabbri@iprj.uerj.br> or any of the other authors for requests.

This work was developed at Brown University, University of Liverpool and State University
of Rio de Janeiro.
