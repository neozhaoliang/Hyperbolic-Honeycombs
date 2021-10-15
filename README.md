> **Prerequisite**: You need to install [Fragmentarium](https://github.com/Syntopia/Fragmentarium) and run the .frag code in fragmentarium.

This project will contain about 6~7 programs that visualize hyperbolic Coxeter groups of rank 4/5 and level 1/2/3. Currently 4 programs are published here. The curious users may refer to [Chen and Labbé's paper](https://arxiv.org/abs/1310.8608) for the math between hyperbolic geometry and ball packings.

> **Note**: the rendering parameters in the code have been set to a moderate level, to render high quality images you may need to tweak with the following tabs on the right:
>
> + MaxRaySteps
> + FudgeFactor
> + Uncomment the `#define USE_IQ_CLOUDS` and `#define KN_VOLUMETRIC` in the header


## 2D hyperbolic tilings (rank = 3, level = 1)

[Shadertoy live](https://www.shadertoy.com/view/7dcXDB)

<img src="https://user-images.githubusercontent.com/23307174/137338871-2a885fe3-1574-4e76-a2f4-a118d06ab351.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/137339107-a63d5ec2-20cb-4798-a58d-1cb094dbe47c.png" width="45%"></img> 

## 3D hyperbolic honeycombs (rank = 4, level = 1, 2)

(Images with "holes" on the boundary are of level 2)

<img src="https://user-images.githubusercontent.com/23307174/131051401-d97b7836-fd49-48b2-826a-7cec3d9e5977.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/131051431-73d9a58c-4348-490c-bfb5-dcad80904b22.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/131051442-8b0120cc-e9c0-4e4d-b0d1-f78b9617b091.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/131051455-32b2f0fa-0ff2-4c8e-9abb-a17aa2087bf9.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/131051476-43bc562e-6d26-4a16-a368-87157d2af8ab.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/131051485-bdfc464b-fbc2-4844-be36-e91207c95ba7.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/131051505-4654aa98-ae38-49f1-b397-d38ea8de602f.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/131051509-2c29ca93-ff91-4a11-9aa5-5b7258b94e4e.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/131051514-56d85ee1-f52e-455b-b6a9-02d710540f44.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/131051519-2a5db2a2-af77-4b07-a4ea-b9b38c82fd6b.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/131051528-5dfab8ad-cbd2-48a4-b114-b4e1bd195dfc.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/131051669-17a7ee93-cc9b-44c0-b176-b603f9ab3404.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/131051683-283114fd-87c9-4e04-a18e-16582015626e.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/131051679-55074bf4-9945-464b-985b-c4be053c473c.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/137326046-8a239cab-6760-4cfc-8ed7-67bf80b0ba67.png" width="45%"></img> 

## 2D circle packings (rank = 4, level = 2)

[Shadertoy live](https://www.shadertoy.com/view/WdGBz3)

<img src="https://user-images.githubusercontent.com/23307174/136644057-cf42394a-a407-4525-8b8f-e2678926458d.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/136644066-60520470-9bc0-4e48-b7c8-dc32dad9b0c8.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/136644075-e23ede99-e935-423c-a8ef-e980606b8a1e.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/136644079-e3dd2b26-228d-41a5-93e6-f66e228a618d.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/136644083-f671fc91-84ac-4d16-a235-095397e11266.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/136644087-f3da996b-429e-4b31-9586-2deef957c25f.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/136649503-9611cb87-bc68-4552-9bf6-6f9a71531029.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/136649676-8a58050d-2312-4662-8f38-83d4469ddfe6.png" width="45%"></img> 

### Circle packings from platonic solids

In order (left to right, top to bottom): tetrahedron, cube, octahedron, dodecahedron, icosahedron.

<img src="https://user-images.githubusercontent.com/8331208/137447759-f7c71794-1a45-4c07-b96e-0a46f176c0f3.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/8331208/137447819-a03b7915-4bac-409a-abe6-c8fa349f9ecf.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/8331208/137447894-3c032241-52ec-4e57-a331-885c7bac551f.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/8331208/137447957-71a1b37a-8bae-4b95-9f85-6235ae58f514.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/8331208/137447648-a3f7b426-d534-4abf-b5f9-061233d468cb.png" width="45%"></img> 

### Non-reflective circle packings

These packings follow from [a preprint of Kapovich and Kontorovich](https://arxiv.org/abs/2104.13838).  Level not defined.

Top: Extended Bianchi groups. Left: Bi23. Right: Bi31.

Bottom: Groups from [Mcleod's thesis](http://etheses.dur.ac.uk/7743/1/thesis31072013.pdf).  Left: Modified f(3,6).  Right: f(3,14).

<img src="https://user-images.githubusercontent.com/8331208/137448747-7ddecdb0-351d-4941-8d22-fc6f9246dd8b.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/8331208/137448807-379af158-2898-4c78-9d5e-fa03b8cb13ae.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/8331208/137448899-feeefc6a-0206-47a3-935d-dabd30389549.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/8331208/137448994-b7a0734f-3dc3-460b-be3a-98ef5bd074af.png" width="45%"></img>

## 2D slices of 3D ball packings (rank = 5, level = 2)

[Shadertoy live](https://www.shadertoy.com/view/NdK3zy)

<img src="https://user-images.githubusercontent.com/23307174/134768877-17c234ac-9ca4-4db9-a8e0-1f10e25151eb.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/134768882-6d7dba22-8dd6-4d36-a43d-cd7760876c1d.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/134768887-690e9973-0ecb-4533-bce5-1bbd206fea62.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/134768890-2b0ae7ab-608f-4c19-81bd-1eb7f48a5f38.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/134768892-f49fc79a-bd8d-41ab-99ee-e2ce127a2541.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/134768895-4dc758e7-155a-465e-a019-829e101a27da.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/134768897-463a0959-1868-4c53-b71a-18a3679831cb.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/134768900-52370daa-cf21-4fe5-b6c0-93798d240a10.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/134841424-9f4bc005-9fc8-41b2-9df3-20408e66af78.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/134841439-7e0bc918-cb34-4b53-912e-d62dd376fe8e.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/134841448-44da62f3-59ae-4640-8688-6a00b33ae714.jpeg" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/134841453-b158c1bb-06de-4e18-ae34-0739504b2d50.jpeg" width="45%"></img> 

## Fractals from 3D ball clusters (rank = 5, level = 3)

<img src="https://user-images.githubusercontent.com/23307174/134768380-4ac1abe1-fb39-4a16-8b4d-b4152ebb7c62.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/134768398-c22017da-7aac-497a-8b76-7cdc1732ccd7.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/134768414-3ec34c1f-2728-4570-8348-6e05d73bed8c.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/134768433-92a87e12-32f2-4685-89a4-b0cc20d05e61.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/134807646-1fa40b62-927e-479c-b9dd-dddcc9c6ad26.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/134812638-1fa41d7f-cf64-4d7d-bab6-a0345e1cc74e.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/134965720-617ca96c-8f5f-4977-a9c8-1f8282e0cea3.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/134814347-eee983bf-9b46-4121-96a5-9dd1274baa0a.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/23307174/137325726-f022570c-4d3d-4339-8356-bb05117c81a1.png" width="45%"></img>

# Authors

+ [Chen Hao](https://twitter.com/Chen_Hao)
+ [Zhao Liang](https://twitter.com/neozhaoliang)
+ [Abdelaziz Nait Merzouk](https://twitter.com/FfKnighty)
