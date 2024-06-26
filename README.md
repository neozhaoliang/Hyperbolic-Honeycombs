> **Requirement**: This repository is already included in the official release of [FragM](https://github.com/3Dickulus/FragM). You can run the examples by installing FragM and navigating to `examples >> neozhaoliang`.

In this project, We will visualize hyperbolic [Coxeter groups](https://en.wikipedia.org/wiki/Coxeter_group) of varying ranks 3/4/5 and levels 1/2/3. The scenes can be categorized into two types:

1. Tiling display: Show the tiling of hyperbolic honeycombs inside the space in the Poincaré ball model and upper half space model.
2. Sphere packing display: Show the sphere packing on the ideal boundary. The complement of this sphere packing is called the limit set.

The level of a Coxeter group $G$ is defined as the smallest non-negative integer $l$, such that after removing any $l$ vertices from its [Coxeter diagram](https://en.wikipedia.org/wiki/Coxeter%E2%80%93Dynkin_diagram), the remaining diagram is either finite or affine. As a result, finite (spherical) and affine (Euclidean) Coxeter groups have level 0.

It's proved in a paper by George Maxwell that Coxeter groups of level 1/2 are both hyperbolic. For level 1, the limit set is the entire ideal boundary and there is no sphere packing. For level 2, there is a maximal sphere packing on the ideal boundary, meaning the spheres do not intersect and fill the boundary completely. For levels higher than 2, the spheres still fill the boundary, but they will overlap. For further mathematical details, please refer to the paper by Chen and Labbé ([Chen and Labbé's paper](https://arxiv.org/abs/1310.8608)) on the connection between hyperbolic geometry and ball packings.

## 3D Euclidean tilings (rank = 4, level = 0)

[Shadertoy live](https://www.shadertoy.com/view/3tccWf)

<div align="center">
<img src="https://user-images.githubusercontent.com/23307174/218312165-5377d195-4acd-4c5b-bfee-bcb15b92bc7b.png"></img>
</div>

## 2D hyperbolic tilings (rank = 3, level = 1, 2)

[Shadertoy live](https://www.shadertoy.com/view/7dcXDB)

From left to right: compact tiling, paracompact tiling (with ideal vertices on the boundary), non-compact tiling (with hyperideal vertices outside the space)

<img src="https://user-images.githubusercontent.com/23307174/137573698-507e1abc-bb23-4a9c-b5ac-0a4bb12e6527.png" width="30%"></img> <img src="https://user-images.githubusercontent.com/23307174/137573695-eb58bf45-dbac-499e-a04e-c74a967f0a42.png" width="30%"></img> <img src="https://user-images.githubusercontent.com/23307174/137573687-7cdfa82c-6129-4599-bc61-ec98f0c838d2.png" width="30%"></img>

The level 2 case, shown in the rightmost image, appears less attractive. However, it can be observed that each cell, which is an unbounded triangle, intersects the ideal boundary at an arc. All these arcs pack the entire boundary circle. This phenomena generalizes to three and four-dimensional spaces. If the group has level 2, each cell in the honeycomb will intersect the boundary at a disk or sphere, and these disks/spheres pack the entire boundary.

## 3D hyperbolic honeycombs (rank = 4, level = 1, 2)

(Images with Poincaré disks packing the boundary are of level 2)

<img src="https://user-images.githubusercontent.com/23307174/219947464-9c20a496-8625-411a-a047-4d2ed074cc26.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/131051401-d97b7836-fd49-48b2-826a-7cec3d9e5977.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/131051431-73d9a58c-4348-490c-bfb5-dcad80904b22.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/131051455-32b2f0fa-0ff2-4c8e-9abb-a17aa2087bf9.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/131051476-43bc562e-6d26-4a16-a368-87157d2af8ab.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/131051485-bdfc464b-fbc2-4844-be36-e91207c95ba7.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/131051505-4654aa98-ae38-49f1-b397-d38ea8de602f.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/131051509-2c29ca93-ff91-4a11-9aa5-5b7258b94e4e.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/131051514-56d85ee1-f52e-455b-b6a9-02d710540f44.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/131051519-2a5db2a2-af77-4b07-a4ea-b9b38c82fd6b.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/131051528-5dfab8ad-cbd2-48a4-b114-b4e1bd195dfc.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/131051669-17a7ee93-cc9b-44c0-b176-b603f9ab3404.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/191715828-a935680c-d13f-4105-9fca-6e84a136f16f.jpeg" width="45%"></img>
<img src="https://user-images.githubusercontent.com/23307174/131051679-55074bf4-9945-464b-985b-c4be053c473c.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/137326046-8a239cab-6760-4cfc-8ed7-67bf80b0ba67.png" width="45%"></img>
<img src="https://user-images.githubusercontent.com/23307174/191716449-765efaea-b4ef-4cc5-879b-7ba285dca436.jpeg" width="45%"></img>

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

Top row: level 2 groups give dense ball packings of the unit ball.

Second row: level > 2 groups have overlapping balls, they give fractal patterns if some of the balls are removed. Basically these are the fratals in the next section but shown in the Poincaré unit ball model.

<img src="https://user-images.githubusercontent.com/23307174/137572582-76bdb60c-7835-4e64-aec2-7f6dee52881d.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/137572621-5f2fd253-1d66-491c-97d2-33d738f33013.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/137572850-d7bf9bfa-c387-4f3f-9b6c-498ef425428e.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/137572851-5eba3051-d2da-47c6-8f70-a0623c93fb14.png" width="45%"></img>

## Fractals from 3D ball clusters (rank = 5, level = 3)

Note some Coxeter diagrams for the rendered images below are missing (I forgot them).

|   |    |
|:---:|:---:|
| ![default](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/81507748-b826-4295-9d20-7bfeba4f0ff8) |  Missing    |
|<img src="https://user-images.githubusercontent.com/23307174/134768380-4ac1abe1-fb39-4a16-8b4d-b4152ebb7c62.png" ></img>| <img src="https://user-images.githubusercontent.com/23307174/134768398-c22017da-7aac-497a-8b76-7cdc1732ccd7.png"></img>|
| ![pink-balls](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/16f452d1-bfc1-4a3f-870f-ea6656c1af93) | ![dense-bronze](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/6987f7ac-d344-473e-ac24-291090651f5b)   |
| <img src="https://user-images.githubusercontent.com/23307174/134768414-3ec34c1f-2728-4570-8348-6e05d73bed8c.png"></img> | <img src="https://user-images.githubusercontent.com/23307174/134768433-92a87e12-32f2-4685-89a4-b0cc20d05e61.png"></img>  |
| ![volcano-lake](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/d141d96c-0af6-4021-ae48-49306432a357) | ![Vinberg](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/bfe802ea-3124-4133-b007-27f2e8ae2386) |
| <img src="https://user-images.githubusercontent.com/23307174/220014316-467f01fc-58b3-4149-aabc-93dbf0d16c2c.png"></img> | <img src="https://user-images.githubusercontent.com/23307174/134812638-1fa41d7f-cf64-4d7d-bab6-a0345e1cc74e.png"></img> |
| Missing |  Missing |
| <img src="https://user-images.githubusercontent.com/23307174/134965720-617ca96c-8f5f-4977-a9c8-1f8282e0cea3.png"></img> | <img src="https://user-images.githubusercontent.com/23307174/134814347-eee983bf-9b46-4121-96a5-9dd1274baa0a.png"></img> |
| ![knighty-cave](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/0b64a3a3-841e-4d80-9ab6-24d4f1772a88) | Missing |
|<img src="https://user-images.githubusercontent.com/23307174/137325726-f022570c-4d3d-4339-8356-bb05117c81a1.png"></img> | <img src="https://user-images.githubusercontent.com/23307174/217011897-71ed5747-7659-4208-a1f2-fd135bae47ad.png"></img>|
| ![golden-2](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/7b0f533a-c373-4dc4-800d-b02a1cfa5dfc) | ![4-4-3-inf-inf](https://github.com/neozhaoliang/Hyperbolic-Honeycombs/assets/23307174/0847fa79-bf96-439e-a782-3a43c4de59fb) |
| <img src="https://user-images.githubusercontent.com/23307174/288715495-99a89c94-3d58-46bb-97a0-29ff1b54c8d2.png"></img> | <img src="https://user-images.githubusercontent.com/23307174/288643505-e0cc5429-b3f8-4fbe-a343-f98a32154255.png"></img>|
| <img src="https://user-images.githubusercontent.com/23307174/236654674-4d3256e9-4296-4694-a8a7-564df4fa8f51.png"></img> | <img src="https://user-images.githubusercontent.com/23307174/236654678-96277ced-03b6-459d-922a-89a3024f2995.png"></img>|


# Authors

+ [Chen Hao](https://twitter.com/Chen_Hao)
+ [Zhao Liang](https://twitter.com/neozhaoliang)
+ [Abdelaziz Nait Merzouk](https://twitter.com/FfKnighty)

# License

The .frag code written for FragM in this repository is licensed under the [GPL License](./LICENSE). The images demonstrated by the authors in this project, including those uploaded by the authors on other platforms such as Twitter, are licensed under the [CC BY-NC-SA license](https://creativecommons.org/licenses/by-nc-sa/4.0/).
