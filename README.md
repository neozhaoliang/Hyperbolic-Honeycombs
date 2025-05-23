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

<img src="https://user-images.githubusercontent.com/23307174/137573698-507e1abc-bb23-4a9c-b5ac-0a4bb12e6527.png" width="30%"></img> <img src="https://user-images.githubusercontent.com/23307174/137573695-eb58bf45-dbac-499e-a04e-c74a967f0a42.png" width="30%"></img> <img src="https://user-images.githubusercontent.com/23307174/137573687-7cdfa82c-6129-4599-bc61-ec98f0c838d2.png" width="30%"></img>

The level 2 case in the rightmost image appears less attractive. However, it can be observed that each cell, which is an unbounded triangle, intersects the ideal boundary at an arc. All these arcs pack the entire boundary circle. This phenomenon generalizes to three and four-dimensional spaces. If the group has level 2, each cell in the honeycomb will intersect the boundary at a disk/sphere, and these disks/spheres pack the entire boundary.

## 3D hyperbolic honeycombs (rank = 4, level = 1, 2)


The code used to render the following image is [here](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/blob/main/src/Hyperoblic-Honeycombs-Kn2cr11.frag). It can render any hyperbolic group of rank 4 that has all labels $m_{s,t}<\infty$ and contains a finite sub-diagram of rank 3. (Images with Poincaré disks packing the boundary are of level 2)

|     |     |
|:---:|:---:|
| ![353-1000](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/3bd6f8ae-db6a-41a6-ab8e-5c8d8820c475)|![363-0100](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/bb08fcca-eab3-4df0-a685-6c4be9738c63)|
|![373-0101](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/3f893a0e-7f3c-4ff2-8442-c243681c837a)|![444-0011](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/b7458b11-11ee-4399-80b3-76c3f41d3915)|
|![445-0011](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/264196a0-5e77-4340-b2c0-384f7764c3a6)|![445-1100](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/dfdeab76-3859-44af-ab8c-618662730ff2)|
|![534-1000](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/442411c1-f880-4b4c-b9ed-d075b8334e15)|![535-1100](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/c47b5d88-6982-4a11-9508-0b199491f387)|
|![536-1100](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/63bf509d-79a0-453d-b5c7-04ae9c20d894)|![735-0011](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/90f1ee50-02c4-4a01-a91d-af39e4e731b5)|
|![522332-1010](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/9b57bca7-8dec-47db-b578-f500eb743a20)|![935-0011](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/3cc05c9a-1685-446c-bd0f-042b058590d7)|
|![532233](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/6ccfdf8b-b5d3-4794-92c5-6e75ae6c5464)|![365-1101](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/533b37a3-18ad-44b7-94ec-c7b4b4a45538)|
|![445-0100](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/2527448a-db7e-439f-8adf-0009e12f3e0e)|![454-1101(1)](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/a70de753-88a4-4b4b-8239-ef5b9f8ccf3b)|
![535-1000](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/ebeed494-72a4-4936-9aff-560172853136)|![522333](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/99742433-98e3-4e67-a50f-c0c63f750455)|
|![534-1100](https://github.com/user-attachments/assets/6261eeba-5f32-4cc0-9f7a-83ba9aa90843)||



## 2D circle packings (rank = 4, level = 2)

[Shadertoy live](https://www.shadertoy.com/view/WdGBz3)

|    |    |
|:---:|:---:|
|![cp1](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/714b9bf8-7653-479e-b9a7-1a93b8a10554)|![cp2](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/564f0ce1-9ce8-47cd-8b6a-63696221170c) |
|<img src="https://user-images.githubusercontent.com/23307174/218310651-b8b2de42-e72f-4695-a398-30a1ff00ecdc.png"></img> | <img src="https://user-images.githubusercontent.com/23307174/218310665-9ac60e78-9981-48e8-9097-08c421d92a67.png"></img> |

<img src="https://user-images.githubusercontent.com/23307174/218777067-774d763c-7377-421b-8941-0f6c34d6ff3c.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/218777081-05e686bb-9f25-40dc-9b4d-755e93fdc0c0.png" width="45%"></img>
<img src="https://user-images.githubusercontent.com/23307174/218783208-a0376e5f-7f2c-48d3-9242-c01b0fc85693.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/218783238-b8746c11-5939-43cc-8f7c-0688aa098da6.png" width="45%"></img>

## 2D circles packings (rank = 4, level > 2)

In this case, there will be overlapping circles:

<div align="center">
<img src="https://user-images.githubusercontent.com/23307174/219011578-75b156ad-cd2f-45b6-9231-91110a629595.png" width="70%"></img>
</div>


### Circle packings from platonic solids

In order (left to right, top to bottom): tetrahedron, cube, octahedron, dodecahedron, icosahedron.

[Shadertoy Live](https://www.shadertoy.com/view/7dcXWs)

<img src="https://user-images.githubusercontent.com/8331208/137447759-f7c71794-1a45-4c07-b96e-0a46f176c0f3.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/8331208/137447819-a03b7915-4bac-409a-abe6-c8fa349f9ecf.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/8331208/137447894-3c032241-52ec-4e57-a331-885c7bac551f.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/8331208/137447957-71a1b37a-8bae-4b95-9f85-6235ae58f514.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/8331208/137447648-a3f7b426-d534-4abf-b5f9-061233d468cb.png" width="45%"></img>

### Non-reflective circle packings

These packings follow from [a preprint of Kapovich and Kontorovich](https://arxiv.org/abs/2104.13838).  Level not defined.

Extended Bianchi groups. Left: [Bi23](https://www.shadertoy.com/view/NddSWn). Right: [Bi31](https://www.shadertoy.com/view/Nd3XzN).

<img src="https://user-images.githubusercontent.com/8331208/137448747-7ddecdb0-351d-4941-8d22-fc6f9246dd8b.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/8331208/137448807-379af158-2898-4c78-9d5e-fa03b8cb13ae.png" width="45%"></img>

Groups from [Mcleod's thesis](http://etheses.dur.ac.uk/7743/1/thesis31072013.pdf).  Left: [Modified f(3,6)](https://www.shadertoy.com/view/sscSDr).  Right: [f(3,14)](https://www.shadertoy.com/view/7scXWn).

<img src="https://user-images.githubusercontent.com/8331208/137448899-feeefc6a-0206-47a3-935d-dabd30389549.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/8331208/137448994-b7a0734f-3dc3-460b-be3a-98ef5bd074af.png" width="45%"></img>

## 2D slices of 3D ball packings (rank = 5, level = 2)

[Shadertoy live](https://www.shadertoy.com/view/NdK3zy)

<img src="https://user-images.githubusercontent.com/23307174/134768877-17c234ac-9ca4-4db9-a8e0-1f10e25151eb.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/134768882-6d7dba22-8dd6-4d36-a43d-cd7760876c1d.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/134768887-690e9973-0ecb-4533-bce5-1bbd206fea62.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/134768890-2b0ae7ab-608f-4c19-81bd-1eb7f48a5f38.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/134768892-f49fc79a-bd8d-41ab-99ee-e2ce127a2541.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/134768895-4dc758e7-155a-465e-a019-829e101a27da.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/134768897-463a0959-1868-4c53-b71a-18a3679831cb.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/134768900-52370daa-cf21-4fe5-b6c0-93798d240a10.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/134841424-9f4bc005-9fc8-41b2-9df3-20408e66af78.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/134841439-7e0bc918-cb34-4b53-912e-d62dd376fe8e.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/134841448-44da62f3-59ae-4640-8688-6a00b33ae714.jpeg" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/134841453-b158c1bb-06de-4e18-ae34-0739504b2d50.jpeg" width="45%"></img>

## 3D ball packings (rank = 5, level >= 2)

These are the ball packings in the next section but shown in the Poincaré unit ball model.


|     |     |    |
|:---:|:---:|:---:|
|![236-323-423-2](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/83d16468-fb49-4d04-b301-68160b3bf867)|![244-224-243-2](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/c8a58983-a578-4291-a17f-d5211dff14a9)|![244-232-425-2](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/85bda303-d67e-4c4f-ba18-237fe6a82faa)|



## Fractals from 3D ball clusters (rank = 5, level = 3)

|     |     |
|:---:|:---:|
|![4-4-inf-inf(2)(1)](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/af669f73-b3bc-44d7-91be-e1a990a94f17)|![4-4-inf-inf(1)](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/4c6df235-9a9e-4ef8-a6c0-ca6171969e9a)|
|![236-223-227-5](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/e544cc67-2c7e-4151-a896-a704347faee1)|![236-444-322-5](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/03fd4e61-50e0-43fd-a9a2-7a532b40c240)|
|![244-223-22inf-inf](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/f57346eb-d9f8-4f65-866d-0c54bd0bebe1)|![244-234-334-4](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/16b5c836-3ca6-44ea-8d22-2b8e127b7cf7)|
|![244-442-323-3](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/aef0a32b-9386-4903-a003-48f91ba8457c)|![244-327-327-4](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/a85eb266-f8fd-4949-9c88-07c65bc3cbe5)|
|![333-224-22inf-inf](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/f1d68ae6-0a7f-4018-a778-748a13aa760a)|![333-225-32inf-inf](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/9b62f444-974f-4fb3-990d-131f8c0c88ca)|
|![333-227-225-inf(inf=1 3)](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/5ea8b446-b2df-4a8d-a859-7d0ab05e688f)|![333-227-226-7](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/e07defd6-8c88-485f-82b1-26a9cd46c989)|
![333-433-224-2](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/93bcd53c-0154-4ab7-8159-faedf2585010)|![333-442-343-3(3D)](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/1f4e99f8-a7d8-46e9-98d4-d95ae4f8105a)|


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
