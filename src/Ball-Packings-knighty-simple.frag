#version 330
/*
A basic Boyd-Maxwell ball packing example, also widely known as the
"pseudo Kleinian fractal"

Original creator: https://github.com/knightyFF

The underlying Coxeter group is

	o---4---o---4---o---inf---o---inf---o

The `FoldSize` parameter controls the size of a fundamental region in the x/y/z directions.
*/

#define KN_VOLUMETRIC
#define USE_EIFFIE_SHADOW
#define MULTI_SAMPLE_AO
#include "MathUtils.frag"
uniform float time;
#include "DE-Kn2cr11.frag"

#group BallPackings
uniform int Max_Iterations;  slider[0,20,100]
uniform float InvSphereRadius;  slider[0,1,2]
uniform vec3 FoldSize;  slider[(0,0,0),(1,1,1),(2,2,2)]

float DE(vec3 p) {
    p.z += FoldSize.z;
    float DEfactor = 1.;
    vec3 ap;
    int i = 0;
    while(i++ < Max_Iterations && ap!=p) {
        ap = p;
        p = 2.0 * clamp(p, -FoldSize, FoldSize) - p;
        float r2 = dot(p,p);
        orbitTrap = min(orbitTrap, abs(vec4(p, r2)));
        float k = max(InvSphereRadius / r2, 1.0);
        p *= k;
        DEfactor *= k;
    }
    return abs(p.z) / DEfactor * 0.4;
}


#preset Default
FOV = 0.4
Eye = 1.4,0.2,-0.5
Target = 1.4,1.2,-0.5
Up = 0,0,1
FocalPlane = 1.46621
Aperture = 0.00619
InFocusAWidth = 0.35484
ApertureNbrSides = 10 NotLocked
ApertureRot = 130.802
ApStarShaped = true NotLocked
Gamma = 2.2
ToneMapping = 5
Exposure = 1
Brightness = 1
Contrast = 1
Saturation = 1
GaussianWeight = 1
AntiAliasScale = 1.5
Bloom = false
BloomIntensity = 0.16666
BloomPow = 0.5618
BloomTaps = 5
Detail = -4.00064
RefineSteps = 4
FudgeFactor = 0.8
MaxRaySteps = 500
MaxDistance = 85.53
Dither = 0.52542
NormalBackStep = 1
DetailAO = -0.85715
coneApertureAO = 0.66129
maxIterAO = 10
AO_ambient = 0.88096
AO_camlight = 0
AO_pointlight = 0
AoCorrect = 0
Specular = 0.01176
SpecularExp = 27.275
CamLight = 1,1,1,0
AmbiantLight = 1,1,1,0.70588
Glow = 1,1,1,0
GlowMax = 0
Reflection = 0.380392,0.380392,0.380392
ReflectionsNumber = 3 Locked
SpotGlow = true
SpotLight = 1,0.972549,0.807843,10
LightPos = -2.9746,3.038,0.5064
LightSize = 0.0297
LightFallOff = 0.49438
LightGlowRad = 0.3226
LightGlowExp = 0.9524
HardShadow = 1 Locked
ShadowSoft = 0
BaseColor = 0.380392,0.380392,0.380392
OrbitStrength = 0.35065
X = 0.129412,0.239216,0.6,1
Y = 1,0.247059,0.0196078,0.99416
Z = 0.384314,1,0.423529,1
R = 0.792157,0.905882,1,1
BackgroundColor = 0.25098,0.317647,0.313725
GradientBackground = 0.3
CycleColors = true
Cycles = 2.70562
EnableFloor = false
FloorNormal = 0,0,0
FloorHeight = 0
FloorColor = 1,1,1
HF_Fallof = 1.62651
HF_Const = 0.08631
HF_Intensity = 0
HF_Dir = 0,0,1
HF_Offset = -3.253
HF_Color = 0.701961,0.803922,0.956863,0.24273
HF_Scatter = 1.577
HF_Anisotropy = 0,0,0,0
HF_FogIter = 2
HF_CastShadow = true
CloudScale = 1
CloudFlatness = 0
CloudTops = 1
CloudBase = -1
CloudDensity = 1
CloudRoughness = 1
CloudContrast = 1
CloudColor = 0.65,0.68,0.7
SunLightColor = 0.7,0.5,0.3
Max_Iterations = 20
InvSphereRadius = 1
FoldSize = 0.94,0.8,0.9
#endpreset



#preset knighty-original
FOV = 0.4
Up = 0,0,1
Gamma = 2.2
Exposure = 1
Brightness = 1
Contrast = 1
Saturation = 1
GaussianWeight = 1
Bloom = false
BloomPow = 0.5618
Detail = -4.00064
RefineSteps = 4
FudgeFactor = 0.8
MaxDistance = 85.53
NormalBackStep = 1
DetailAO = -0.85715
coneApertureAO = 0.66129
maxIterAO = 10
AoCorrect = 0
Specular = 0.01176
SpecularExp = 27.275
Reflection = 0.380392,0.380392,0.380392
ReflectionsNumber = 3 Locked
SpotGlow = true
SpotLight = 1,0.972549,0.807843,10
LightSize = 0.0297
LightFallOff = 0.49438
LightGlowRad = 0.3226
LightGlowExp = 0.9524
HardShadow = 1 Locked
BaseColor = 0.380392,0.380392,0.380392
OrbitStrength = 0.35065
X = 0.129412,0.239216,0.6,1
Y = 1,0.247059,0.0196078,0.99416
Z = 0.384314,1,0.423529,1
R = 0.792157,0.905882,1,1
BackgroundColor = 0.25098,0.317647,0.313725
GradientBackground = 0.3
CycleColors = true
Cycles = 2.70562
EnableFloor = false
FloorNormal = 0,0,0
FloorHeight = 0
FloorColor = 1,1,1
HF_Fallof = 1.62651
HF_Const = 0.08631
HF_Intensity = 0
HF_Dir = 0,0,1
HF_Offset = -3.253
HF_Color = 0.701961,0.803922,0.956863,0.24273
HF_FogIter = 2
HF_CastShadow = true
CloudScale = 1
CloudContrast = 1
CloudColor = 0.65,0.68,0.7
SunLightColor = 0.7,0.5,0.3
Eye = 1.41939,0.540671,-0.7
Target = 1.4,1.2,-0.7
EquiRectangular = false
FocalPlane = 1.46621
Aperture = 0.00619
InFocusAWidth = 0.35484
DofCorrect = true
ApertureNbrSides = 10
ApertureRot = 130.802
ApStarShaped = true
ToneMapping = 5
AntiAliasScale = 1.5
BloomIntensity = 0.16666
BloomTaps = 5
BloomStrong = 1
LensFlare = false
FlareIntensity = 0.25
FlareSamples = 8
FlareDispersal = 0.25
FlareHaloWidth = 0.5
FlareDistortion = 1
DepthToAlpha = false
ShowDepth = false
DepthMagnitude = 1
MaxRaySteps = 500
Dither = 0.52542
FudgeAO = 1
AO_ambient = 0.88096
AO_camlight = 0
AO_pointlight = 0
CamLight = 1,1,1,0
AmbiantLight = 1,1,1,0.70588
LightPos = -2.9746,3.038,0.5064
ShadowSoft = 0
ShadowBlur = 0
perf = false
SSS = false
sss1 = 0.1
sss2 = 0.5
HF_Scatter = 1.577
HF_Anisotropy = 0,0,0
EnCloudsDir = false
CloudDir = 0,0,1
CloudOffset = 0,0,0
CloudFlatness = 0
CloudTops = 1
CloudBase = -1
CloudDensity = 1
CloudRoughness = 1
CloudBrightness = 1
CloudColor2 = 0.07,0.17,0.24
Cloudvar1 = 0.99
Cloudvar2 = 1
CloudIter = 5
CloudBgMix = 1
WindDir = 0,0,1
WindSpeed = 1
Max_Iterations = 20
InvSphereRadius = 1
FoldSize = 0.90276,1.01192,0.988
#endpreset


