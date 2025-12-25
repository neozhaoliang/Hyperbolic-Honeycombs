> **Requirement**: This repository is included in the official release of [FragM](https://github.com/3Dickulus/FragM). You can run the examples by installing FragM and navigating to `examples >> neozhaoliang`.

In this project,  we’ll explore visualizations of hyperbolic [Coxeter groups](https://en.wikipedia.org/wiki/Coxeter_group) across various ranks (3, 4, and 5) and levels (1, 2, and 3). The visualizations are divided into two main types:

1. Tiling display: this showcases the tiling of hyperbolic honeycombs within the space using the Poincaré ball and upper half-space models
2. Sphere packing display: This illustrates the sphere packing on the ideal boundary. The complement of this packing is known as the limit set.

The level of a Coxeter group $G$ is defined as the smallest non-negative integer $l$ such that after removing any $l$ vertices from its [Coxeter diagram](https://en.wikipedia.org/wiki/Coxeter%E2%80%93Dynkin_diagram), all the connected components of the remaining diagram are either finite or affine. For example, finite (spherical) and affine (Euclidean) Coxeter groups both have level 0.

[George Maxwell's work](https://www.sciencedirect.com/science/article/pii/0021869382903180) establishes that Coxeter groups of level 1 or 2 are hyperbolic. For level 1 groups, the limit set is the whole ideal boundary, and no sphere packing is observed. In contrast, level 2 groups exhibit a maximal sphere packing on the ideal boundary, this means the spheres fill the boundary without intersecting. For level higher than 2, the spheres still fill the boundary but they will necessarily overlapp. For further mathematical details, please refer to the paper by Chen and Labbé ([Chen and Labbé's paper](https://arxiv.org/abs/1310.8608)) on the connection between hyperbolic geometry and sphere packings.

## 3D Euclidean tilings (rank = 4, level = 0)

[Shadertoy live](https://www.shadertoy.com/view/3tccWf)

<div align="center">
<img src="https://user-images.githubusercontent.com/23307174/218312165-5377d195-4acd-4c5b-bfee-bcb15b92bc7b.png"></img>
</div>

## 2D hyperbolic tilings (rank = 3, level = 1, 2)

[Shadertoy live](https://www.shadertoy.com/view/7dcXDB)

From left to right: compact tiling, paracompact tiling (with ideal vertices on the boundary), non-compact tiling (with hyperideal vertices outside the space)

<img src="https://neozhaoliang.github.io/assets/3-3-7.png" width="30%"></img> <img src="https://neozhaoliang.github.io/assets/3-3-inf.png" width="30%"></img> <img src="https://neozhaoliang.github.io/assets/3-inf-3(inf=1.15)-Vinberg.png" width="30%"></img>

The level 2 case in the rightmost image appears less attractive. However, it can be observed that each cell, which is an unbounded triangle, intersects the ideal boundary at an arc. All these arcs pack the entire boundary circle. This phenomenon generalizes to three and four-dimensional spaces. If the group has level 2, each cell in the honeycomb will intersect the boundary at a disk/sphere, and these disks/spheres pack the entire boundary.

## 3D hyperbolic honeycombs (rank = 4, level = 1, 2)


The code used to render the following image is [here](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/blob/main/src/Hyperoblic-Honeycombs-Kn2cr11.frag). It can render any hyperbolic group of rank 4 that has all labels $m_{s,t}<\infty$ and contains a finite sub-diagram of rank 3. (Images with Poincaré disks packing the boundary are of level 2)

|     |     |
|:---:|:---:|
| ![353-1000](https://neozhaoliang.github.io/assets/353-1000.png)|![363-0100](https://neozhaoliang.github.io/assets/363-0100.png)|
|![373-0101](https://neozhaoliang.github.io/assets/373-0101.png)|![444-0011](https://neozhaoliang.github.io/assets/444-0011.png)|
|![445-0011](https://neozhaoliang.github.io/assets/445-0011.png)|![445-1100](https://neozhaoliang.github.io/assets/445-1100.png)|
|![534-1000](https://neozhaoliang.github.io/assets/534-1000.png)|![535-1100](https://neozhaoliang.github.io/assets/535-1100.png)|
|![536-1100](https://neozhaoliang.github.io/assets/536-1100.png)|![735-0011](https://neozhaoliang.github.io/assets/735-0011.png)|
|![522332-1010](https://neozhaoliang.github.io/assets/522332-1010.png)|![935-0011](https://neozhaoliang.github.io/assets/935-0011.png)|
|![445-0100](https://neozhaoliang.github.io/assets/445-0100.png)|![454-1101](https://neozhaoliang.github.io/assets/454-1101.png)|
![535-1000](https://neozhaoliang.github.io/assets/535-1000.png)|![522333](https://neozhaoliang.github.io/assets/522333.png)|



## 2D circle packings (rank = 4, level = 2)

[Shadertoy live](https://www.shadertoy.com/view/WdGBz3)

|    |    |
|:---:|:---:|
|![](https://neozhaoliang.github.io/assets/327-624.png)|![](https://neozhaoliang.github.io/assets/334-433.png)|
|![](https://neozhaoliang.github.io/assets/344-4inf2.png)|![](https://neozhaoliang.github.io/assets/433-523.png)|
|![](https://neozhaoliang.github.io/assets/633-642.png)|![](https://neozhaoliang.github.io/assets/6-3-inf.png)|


## 2D circles packings (rank = 4, level > 2)

In this case, there will be overlapping circles:

<div align="center">
<img src="https://neozhaoliang.github.io/assets/434-32inf(inf=1.3).png" width="70%"></img>
</div>


### Circle packings from platonic solids

In order (left to right, top to bottom): tetrahedron, cube, octahedron, dodecahedron, icosahedron.

[Shadertoy Live](https://www.shadertoy.com/view/7dcXWs)

|   |   |
|:---:|:---:|
|![](https://neozhaoliang.github.io/assets/cube.png)|![](https://neozhaoliang.github.io/assets/octahedron.png)|
|![](https://neozhaoliang.github.io/assets/tetrahedron.png)|![](https://neozhaoliang.github.io/assets/dodecahedron.png)|
|![](https://neozhaoliang.github.io/assets/icosahedron.png)||

### Non-reflective circle packings

These packings follow from [a preprint of Kapovich and Kontorovich](https://arxiv.org/abs/2104.13838).  Level not defined.

Extended Bianchi groups. Left: [Bi23](https://www.shadertoy.com/view/NddSWn). Right: [Bi31](https://www.shadertoy.com/view/Nd3XzN).


|   |   |
|:---:|:---:|
| Bi23 |  Bi31 |
|<img src="https://user-images.githubusercontent.com/8331208/137448747-7ddecdb0-351d-4941-8d22-fc6f9246dd8b.png"></img>|<img src="https://user-images.githubusercontent.com/8331208/137448807-379af158-2898-4c78-9d5e-fa03b8cb13ae.png"></img>|

Groups from [Mcleod's thesis](http://etheses.dur.ac.uk/7743/1/thesis31072013.pdf).  Left: [Modified f(3,6)](https://www.shadertoy.com/view/sscSDr).  Right: [f(3,14)](https://www.shadertoy.com/view/7scXWn).

|   |   |
|:---:|:---:|
| Modified f(3,6) |  f(3,14) |
|<img src="https://user-images.githubusercontent.com/8331208/137448899-feeefc6a-0206-47a3-935d-dabd30389549.png"></img>|<img src="https://user-images.githubusercontent.com/8331208/137448994-b7a0734f-3dc3-460b-be3a-98ef5bd074af.png"></img>|

## 2D slices of 3D ball packings (rank = 5, level = 2)

[Shadertoy live](https://www.shadertoy.com/view/NdK3zy)

|     |     |
|:---:|:---:|
![](https://neozhaoliang.github.io/assets/333-323-233-2.png)|![](https://neozhaoliang.github.io/assets/333-522-233-2.png)|
|![](https://neozhaoliang.github.io/assets/333-224-223-2.png)|![](https://neozhaoliang.github.io/assets/244-223-232-5.png)|


## 3D ball packings (rank = 5, level >= 2)

These are the ball packings in the next section but shown in the Poincaré unit ball model.

|     |     |    |
|:---:|:---:|:---:|
|![236-323-423-2](https://neozhaoliang.github.io/assets/236-323-423-2.png)|![244-224-243-2](https://neozhaoliang.github.io/assets/244-224-243-2.png)|![244-234-334-4](https://neozhaoliang.github.io/assets/244-232-425-2.png)|



## Fractals from 3D ball clusters (rank = 5, level = 3)

|     |     |
|:---:|:---:|
|![236-444-322-5](https://neozhaoliang.github.io/assets/236-444-322-5.png)|![244-223-22inf-inf](https://neozhaoliang.github.io/assets/244-223-22inf-inf.png)|
|![244-234-334-4](https://neozhaoliang.github.io/assets/244-234-334-4.png)|![244-442-323-3](https://neozhaoliang.github.io/assets/244-442-323-3.png)|
|![244-327-327-4](https://neozhaoliang.github.io/assets/244-327-327-4.png)|![333-224-22inf-inf](https://neozhaoliang.github.io/assets/333-224-22inf-inf.png)|
|![333-225-32inf-inf](https://neozhaoliang.github.io/assets/333-225-32inf-inf.png)|![333-227-225-inf(inf=1.3)](https://neozhaoliang.github.io/assets/333-227-225-inf(inf=1.3).png)|
|![333-227-226-7](https://neozhaoliang.github.io/assets/333-227-226-7.png)|![333-433-224-2](https://neozhaoliang.github.io/assets/333-433-224-2.png)|
|![4-4-inf-inf](https://neozhaoliang.github.io/assets/4-4-inf-inf.png)|![4-4-inf-inf](https://neozhaoliang.github.io/assets/4-4-inf-inf(2).png)|

# How to use this project in FragM

Please refer to the [official Wiki page](https://github.com/3Dickulus/FragM/wiki) of FragM for more detailed information on how to use it.

+ Download or clone this repository to your local machine.

+ Visit the [Fragmentarium release page](https://github.com/3Dickulus/FragM/releases/tag/v2.5.7-221224) and select the appropriate release for your operating system. In this tutorial, the instructions are based on a Windows environment. Therefore, download the file `Fragmentarium-2.5.7-221224-winex.7z`. Save the file and extract it to a convenient location on your disk.

    ![Screenshot releases](https://github.com/user-attachments/assets/07008d02-f83d-46b1-8ecd-7321aa861e11)

+ In the extracted folder, locate the executable file named `Fragmentarium-2.5.7.exe`. Double-click it to launch the application. Upon launching, you should see the following interface:

    ![Screenshot gui)](https://github.com/user-attachments/assets/2ccae9d4-51e8-480f-80b8-fb9a7712d48d)

    The interface is organized into four main regions:

    1. The central area displays the rendered result based on the loaded .frag file.
    2. On the left side, you will find the code editor. If you make changes to the source code, press `Ctrl + S` to save your modifications, then click `Build` to recompile and view the updated result.
    3. The right side houses the control panel, where you can adjust various parameters. These controls are defined within the .frag file using `#group` macros.
    4. The bottom section is dedicated to logging. If the code fails to compile, check the error messages here for troubleshooting information.

+ From the menu bar, select `File -> Open`. Navigate to the directory where you saved the source code of this project and choose a .frag file. For example, select `Ball-Packings-UHS.frag`. Fragmentarium will load and compile the file, displaying the rendered output on your screen:

    ![Screenshot render](https://github.com/user-attachments/assets/d06d480f-69b4-4c89-ade7-15f383a8315b)


# Authors

+ [Chen Hao](https://twitter.com/Chen_Hao)
+ [Zhao Liang](https://twitter.com/neozhaoliang)
+ [Abdelaziz Nait Merzouk](https://twitter.com/FfKnighty)

# License

The .frag code written for FragM in this repository is licensed under the [GPL License](./LICENSE). The images demonstrated by the authors in this project, including those uploaded by the authors on other platforms such as Twitter, are licensed under the [CC BY-NC-SA license](https://creativecommons.org/licenses/by-nc-sa/4.0/).
