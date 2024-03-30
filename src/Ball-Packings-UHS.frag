#version 330
#define USE_IQ_CLOUDS
#define KN_VOLUMETRIC
#define USE_EIFFIE_SHADOW
#define MULTI_SAMPLE_AO
#define providesInit
#include "MathUtils.frag"
uniform float time;
#include "DE-Kn2cr11.frag"

#group BallPacking-Settings
uniform int euclideanTriangleType; slider[0,1,2]
uniform vec3 dihedralAngles0_234; slider[(1,1,1),(3,2,7),(20,20,20)]
uniform vec3 dihedralAngles1_234; slider[(1,1,1),(3,2,7),(20,20,20)]
uniform float dihedralAngle0_1; slider[1,4,20]
uniform vec4 isRealBall; slider[(0,0,0,0),(1,1,0,0),(1,1,1,1)]
uniform int Iterations; slider[1,50,500]


#define inf           1.0
#define L2(x)         dot(x, x)

#define s2 1.41421356
#define s3 1.73205081

float section_height;

float dihedral(float x) {
    return x == inf ? 1. : cos(PI / x);
}

vec3 dihedral(vec3 v) {
    return vec3(dihedral(v.x), dihedral(v.y), dihedral(v.z));
}

struct Ball {
    bool isplane;
    vec3 n;
    float r;
    bool invert;
    bool isRealBall;
};


// coclusters are mirror balls, they corresponde to root vectors (space-like)
Ball[5] coclusters;
// clusters are real balls, they corresponde to space-like weight vectors
Ball[5] clusters;

Ball defaultBall() {
    return Ball(false, vec3(0, 0, -1), 0., false, false);
}

// Distance from a point to a ball
float sdistanceToBall(vec3 p, Ball B) {
    if (B.isplane) {
        float k = dot(vec4(p, 1), vec4(B.n, B.r));
        return k;
    }
    else
        return length(p - B.n) - B.r;
}

Ball from_plane(vec3 n, float d) {
    return Ball(true, n, d, false, false);
}

Ball from_sphere(vec3 cen, float r) {
    return Ball(false, cen, r, false, false);
}

void invertBall(inout Ball B) {
    B.invert = !B.invert;
}

bool try_reflect(inout vec3 p,
                 in Ball B,
                 inout float scale,
                 inout vec4 orb) {
    if (B.isplane) {
        float k = dot(vec4(p, 1), vec4(B.n, B.r));
        if (k >= 0.)
            return true;
        p -= 2. * k  * B.n;
        return false;
    }
    else {
        vec3 cen = B.n;
        float r = B.r;
        vec3 q = p - cen;
        float d2 = dot(q, q);
        float k = (r * r) / d2;
        if ( (k < 1.0 && B.invert) || (k > 1. && !B.invert) )
            return true;

        orb = min(orb, vec4(abs(p), d2));
        scale*=k;
        p = k * q + cen;
        return false;
    }
}

Ball solveBall(mat3 M, vec3 b) {
    vec3 p = b * inverse(M);
    return from_sphere(vec3(p.xy, 0.), p.z);
}

Ball solveBall(vec2 P, Ball B0, Ball B1) {
    if (B0.isplane) {
        float z = B0.r;
        vec3 cen = vec3(P, z);
        float R = sqrt(L2(cen - B1.n) - B1.r*B1.r);
        return from_sphere(cen, R);
    }
    else {
        float r1 = B1.r;
        float r0 = B0.r;
        float z0 = B0.n.z;
        float k0 = L2(P - B0.n.xy);
        float k1 = L2(P - B1.n.xy);
        float z = (r1*r1 - r0*r0 + z0*z0 + k0 - k1) / (2.*z0);
        float R = sqrt(k1 + z*z - r1*r1);
        return from_sphere(vec3(P, z), R);
    }
}

void init() {

    mat3 M0, M1;
    vec3 b;
    Ball B0, B1, B2, B3, B4;
    vec3 t0 = dihedral(dihedralAngles0_234);
    vec3 t1 = dihedral(dihedralAngles1_234);
    float t01 = dihedral(dihedralAngle0_1);
    // A, B, C are the vertices of the triangle formed by mirror plane v2, v3, v4 and z=0 plane
    vec2 A, B, C;

    // the 236 case
    if (euclideanTriangleType == 0) {
        A = vec2(0, 0), B = vec2(0, s3), C = vec2(1, 0);
        B2 = from_plane(vec3(1, 0, 0), 0.);
        B3 = from_plane(vec3(-s3/2., -0.5, 0), s3/2.);
        B4 = from_plane(vec3(0, 1, 0), 0.);
        M1 = mat3(vec3(1, 0, -t1.x), vec3(s3/2., 0.5, t1.y), vec3(0, 1, -t1.z));
        M0 = mat3(vec3(1, 0, -t0.x), vec3(s3/2., 0.5, t0.y), vec3(0, 1, -t0.z));
        b = vec3(0, s3/2., 0);
    }

    // the 244 case
    else if (euclideanTriangleType == 1) {
        A = vec2(0, 0), B = vec2(0, 1), C = vec2(1, 0);
        B2 = from_plane(vec3(1, 0, 0), 0.);
        B3 = from_plane(vec3(-s2/2., -s2/2., 0), s2/2.);
        B4 = from_plane(vec3(0, 1, 0), 0.);
        M1 = mat3(vec3(1, 0, -t1.x), vec3(1./s2, 1./s2, t1.y), vec3(0, 1, -t1.z));
        M0 = mat3(vec3(1, 0, -t0.x), vec3(1./s2, 1./s2, t0.y), vec3(0, 1, -t0.z));
        b = vec3(0, s2/2., 0);
    }

    // the 333 case
    else {
        A = vec2(-1, 0), B = vec2(0, s3), C = vec2(1, 0);
        B2 = from_plane(vec3(s3/2., -.5, 0), s3/2.);
        B3 = from_plane(vec3(-s3/2., -.5, 0), s3/2.);
        B4 = from_plane(vec3(0, 1, 0), 0.);
        M1 = mat3(vec3(-s3/2., 0.5, t1.x), vec3(s3/2., .5, t1.y), vec3(0, 1, -t1.z));
        M0 = mat3(vec3(-s3/2., 0.5, t0.x), vec3(s3/2., .5, t0.y), vec3(0, 1, -t0.z));
        b = vec3(s3, s3, 0)/2.;
    }

    // now we solve the virtual ball B1, this can't be a plane
    B1 = solveBall(M1, b);
    invertBall(B1);

    // now we solve the virtual ball B0, this can be either a plane or a sphere
    // this depends on if all entries in dihedralAngles0 are all 2
    if (dot(dihedralAngles0_234, vec3(1)) == 6.) {
        B0 = from_plane(vec3(0, 0, -1), B1.r*t01);
    }
    else {
        B0 = solveBall(M0, b);
        float r1 = B1.r, r0 = B0.r;
        B0.n.z = sqrt(r0*r0 + r1*r1 + 2.*r0*r1*t01 - L2(B1.n.xy - B0.n.xy));
        invertBall(B0);
    }
    coclusters = Ball[5] (B0, B1, B2, B3, B4);

    section_height = B0.isplane ? 2.*B0.r : B0.n.z;

    //now we process the real balls
    for (int k = 0; k < 5; k++)
        clusters[k] = defaultBall();

    if (isRealBall.x == 1.) {
        clusters[1] = from_plane(vec3(0, 0, -1.), B0.n.z);
        clusters[1].isRealBall = true;
    }
    if (isRealBall.y == 1.) {
        clusters[2] = solveBall(C, B0, B1);
        clusters[2].isRealBall = true;
    }
    if (isRealBall.z== 1.) {
        clusters[3] = solveBall(A, B0, B1);
        clusters[3].isRealBall = true;
    }

    if (isRealBall.w==1.) {
        clusters[4] = solveBall(B, B0, B1);
        clusters[4].isRealBall = true;
    }
}


float map(inout vec3 p, inout float scale, inout vec4 orb) {
    for (int i = 0; i < Iterations; i++) {
        bool cond = true;
        for (int k = 0; k < 5; k++) {
            cond = cond && try_reflect(p, coclusters[k], scale, orb);
        }
        if (cond)
            break;
    }

    float d = abs(p.z);
    for (int j = 1; j < 5; j++) {
        if (clusters[j].isRealBall) {
            d = min(abs(sdistanceToBall(p, clusters[j])), d);
        }
    }
    return d;
}


float DE(vec3 p) {
    float DEfactor=1.;
    float d = map(p, DEfactor, orbitTrap);

    //Call basic shape and scale its DE
    return 0.25*d/DEfactor;

}

















#preset default
AutoFocus = false
FOV = 0.685022
Eye = 5.52166,-0.807811,2.62172
Target = 1.21024,5.81932,-4.00749
UpLock = false
Up = 0,0,1
AutoFocus = false
FocalPlane = 1
Aperture = 0.003
InFocusAWidth = 1
DofCorrect = true
ApertureNbrSides = 5
ApertureRot = 0
ApStarShaped = true
Gamma = 1
ToneMapping = 5
Exposure = 1
Brightness = 1
Contrast = 1
Saturation = 1
GaussianWeight = 1
AntiAliasScale = 1.5
Bloom = true
BloomIntensity = 0.4331
BloomPow = 6.22093
BloomTaps = 23
BloomStrong = 3.57901
DepthToAlpha = true
Detail = -4
RefineSteps = 5
FudgeFactor = 0.58333
MaxRaySteps = 607
MaxDistance = 200
Dither = 0.5
NormalBackStep = 30.6125
DetailAO = -1.96825
coneApertureAO = 0.378985
maxIterAO = 19
FudgeAO = 0.349943
AO_ambient = 1
AO_camlight = 1.59217
AO_pointlight = 0.449594
AoCorrect = 0
Specular = 0.4
SpecularExp = 200
CamLight = 1,1,1,0.2
AmbiantLight = 0.972549,0.972549,0.972549,0.5
Reflection = 0.352941,0.352941,0.352941
ReflectionsNumber = 3
SpotGlow = false
SpotLight = 1,1,1,0.7
LightPos = 2.7888,-0.0764,10
LightSize = 0
LightFallOff = 0.38016
LightGlowRad = 0
LightGlowExp = 0
HardShadow = 1
ShadowSoft = 19.0118
ShadowBlur = 0
perf = false
SSS = false
sss1 = 0.1
sss2 = 0.5
BaseColor = 0.72549,0.72549,0.521569
OrbitStrength = 1
X = 0.898039,0.937255,0.976471,1
Y = 0.898039,0.937255,0.976471,0
Z = 0.898039,0.937255,0.976471,0
R = 0.898039,0.937255,0.976471,1
BackgroundColor = 0.270588,0.403922,0.6
GradientBackground = 0
CycleColors = true
Cycles = 0.1
EnableFloor = true
FloorNormal = 0,0,1.1
FloorHeight = 1.8
FloorColor = 0.313725,0.313725,0.313725
HF_Fallof = 0.187344
HF_Const = 0.05333
HF_Intensity = 0
HF_Dir = 0,0,1
HF_Offset = 0
HF_Color = 0.564706,0.752941,0.878431,1
HF_Scatter = 10
HF_Anisotropy = 0.133333,0.00784314,0
HF_FogIter = 1
HF_CastShadow = true
EnCloudsDir = true NotLocked
Clouds_Dir = 0.273574,-0.22061,-1 NotLocked
CloudScale = 2.50337 NotLocked
CloudFlatness = 0 NotLocked
CloudTops = 1 NotLocked
CloudBase = -5.8 NotLocked
CloudDensity = 0.484136 NotLocked
CloudRoughness = 1 NotLocked
CloudContrast = 1 NotLocked
CloudColor = 0.65,0.68,0.7 NotLocked
CloudColor2 = 0.07,0.17,0.24 NotLocked
SunLightColor = 0.968627,0.968627,0.968627 NotLocked
Cloudvar1 = 0.99 NotLocked
Cloudvar2 = 1 NotLocked
CloudIter = 3 NotLocked
CloudBgMix = 1 NotLocked
euclideanTriangleType = 2
dihedralAngles0_234 = 4,3,3
dihedralAngles1_234 = 2,2,4
dihedralAngle0_1 = 2
isRealBall = 0,0,0,0
Iterations = 500
#endpreset






#preset another
FOV = 0.685022
Eye = 5.4357,-0.74,2.71194
Target = 5.07660231,-0.111265905,2.04227468
Up = 0,0,1
EquiRectangular = false
FocalPlane = 1
Aperture = 0.003
InFocusAWidth = 1
DofCorrect = true
ApertureNbrSides = 5
ApertureRot = 0
ApStarShaped = true
Gamma = 1
ToneMapping = 5
Exposure = 1
Brightness = 1
Contrast = 1
Saturation = 1
GaussianWeight = 1
AntiAliasScale = 1.5
Bloom = true
BloomIntensity = 0.4331
BloomPow = 6.22093
BloomTaps = 23
BloomStrong = 3.57901
LensFlare = false
FlareIntensity = 0.25
FlareSamples = 8
FlareDispersal = 0.25
FlareHaloWidth = 0.5
FlareDistortion = 1
DepthToAlpha = true
ShowDepth = false
DepthMagnitude = 1
Detail = -4
RefineSteps = 5
FudgeFactor = 0.52778
MaxRaySteps = 738
MaxDistance = 200
Dither = 0.5
NormalBackStep = 10
DetailAO = -1.96825
coneApertureAO = 0.378985
maxIterAO = 19
FudgeAO = 0.349943
AO_ambient = 1
AO_camlight = 1.59217
AO_pointlight = 0.449594
AoCorrect = 0
Specular = 0.4
SpecularExp = 200
CamLight = 1,1,1,0.2
AmbiantLight = 0.972549,0.972549,0.972549,0.5
Reflection = 0.352941,0.352941,0.352941
ReflectionsNumber = 3
SpotGlow = false
SpotLight = 1,1,1,0.7
LightPos = 2.7888,-0.0764,10
LightSize = 0
LightFallOff = 0.38016
LightGlowRad = 0
LightGlowExp = 0
HardShadow = 1
ShadowSoft = 19.0118
ShadowBlur = 0
perf = false
SSS = false
sss1 = 0.1
sss2 = 0.5
BaseColor = 0.72549,0.72549,0.521569
OrbitStrength = 1
X = 0.898039,0.937255,0.976471,1
Y = 0.898039,0.937255,0.976471,0
Z = 0.898039,0.937255,0.976471,0
R = 0.898039,0.937255,0.976471,1
BackgroundColor = 0.270588,0.403922,0.6
GradientBackground = 0
CycleColors = true
Cycles = 0.1
EnableFloor = true
FloorNormal = 0,0,1
FloorHeight = 1.8
FloorColor = 0.313725,0.313725,0.313725
HF_Fallof = 0.187344
HF_Const = 0.05333
HF_Intensity = 0
HF_Dir = 0,0,1
HF_Offset = 0
HF_Color = 0.564706,0.752941,0.878431,1
HF_Scatter = 10
HF_Anisotropy = 0.133333,0.00784314,0
HF_FogIter = 1
HF_CastShadow = true
EnCloudsDir = true
CloudDir = 0,0,1
CloudScale = 2.50337
CloudOffset = 0,0,0
CloudFlatness = 0
CloudTops = 1
CloudBase = -5.8
CloudDensity = 0.484136
CloudRoughness = 1
CloudContrast = 1
CloudBrightness = 1
CloudColor = 0.65,0.68,0.7
CloudColor2 = 0.07,0.17,0.24
SunLightColor = 0.968627,0.968627,0.968627
Cloudvar1 = 0.99
Cloudvar2 = 1
CloudIter = 3
CloudBgMix = 1
WindDir = 0,0,1
WindSpeed = 1
euclideanTriangleType = 2
dihedralAngles0_234 = 4,4,2
dihedralAngles1_234 = 3,4,3
dihedralAngle0_1 = 3
isRealBall = 0,0,0,0
Iterations = 500
#endpreset

