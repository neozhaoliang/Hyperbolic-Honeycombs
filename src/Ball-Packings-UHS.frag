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
uniform int Iterations; slider[1,30,500]
uniform float DEcor; slider[0,1,1]


#define inf           1.0
#define L2(x)         dot(x, x)


//infinite foldings-----------------------------------------------------------//
//fold about line with normal ndir (normalized) and dist to origin d.
void fold(inout vec2 p, in vec2 ndir, in float d) {
    float t = dot(p, ndir) + d;
    t = 2.*min(t,0.);
    p -= t * ndir;
}

//fmod with step = stp
void infmod(inout float x, float stp) {
    x *= 1./stp;
    x -= floor(x);
    x *= stp;
}

//infinite fold. stp is the tile size.
//its like performing infinite folds about plane at 0 and plane at stp
void inffold(inout float x, float stp) {
    x *= 1./(2.*stp);
    x -= floor(x);
    x  = 0.5 - abs( x - 0.5 );
    x *= 2.*stp;
}

//Simplest: fold into unit square then about the diagonal
void fold244(inout vec2 p) {
    inffold( p.x , 1. );
    inffold( p.y , 1. );
#define VAL sqrt(2.)/2.
    fold( p, -vec2( VAL, VAL ), VAL);
#undef VAL
}

//fold into rectangle then a little sequence of line folds.
void fold236(inout vec2 p) {
#define S3 sqrt(3.)
    inffold( p.x , 3. );
    inffold( p.y , S3 );
    vec2 n = -0.5 * vec2(S3, 1.);
    float d = .5 * S3;
    fold(p, n, d);
    p = abs(p);
    fold(p, n, d);
#undef S3
}

//most complicated
void fold333(inout vec2 p) {
#define S3 sqrt(3.)
    //change origin
    p.x += 1.;
    //fold y into segment of height sqrt(3)
    inffold( p.y ,S3 );
    //change of coordinates to go to ( (2,0) ; (1,sqrt(3)) ) ish basis. We only nees x coordinate to be transformed
    p.x -= p.y * 1. / S3;
    //do an fmod instead for x.
    infmod(p.x, 6.);
    //undo the coordinates change
    p.x += p.y * 1./S3;
    //The folding sequence... 4 folds
    //There are other choices. I've choosen the one where I need only one direction instead of two.
    vec2 n = -0.5 * vec2(S3, 1.);
    float d = 2. * S3;
    //1st
    fold(p, n, d);
    //2nd
    d = S3;
    fold(p, n, d);
    //3rd
    p.y = abs(p.y);
    //4th same as 2nd
    fold(p, n, d);
    //restore origin
    p.x -=1.;
#undef S3
}

const float s2 = sqrt(2.);
const float s3 = sqrt(3.);

float dihedral(float x) {
    return x == inf ? inf : cos(PI / x);
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
    orb = min(orb, vec4(abs(p),dot(p,p)));
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
        else {
            scale *= k;
            p = k * q + cen;
            return false;
        }
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

    //now we process the real balls
    for (int k = 0; k < 5; k++) {
        clusters[k] = defaultBall();
    }

    clusters[1] = from_plane(vec3(0, 0, -1.), B0.n.z);
    clusters[1].isRealBall = isRealBall.x == 1.;

    clusters[2] = solveBall(C, B0, B1);
    clusters[2].isRealBall = isRealBall.y == 1.;

    clusters[3] = solveBall(A, B0, B1);
    clusters[3].isRealBall = isRealBall.z == 1.;

    clusters[4] = solveBall(B, B0, B1);
    clusters[4].isRealBall = isRealBall.w == 1.;
}

#if 1
// In this version, invisible (non active) real balls are used to reduce oerstepping.
// One may notice that when activating other balls the ovestepping is educed a lot.
// Not very good but better than nothing. :)
float dist2balls(vec3 p, float scale) {
    float d = abs(p.z); //Ball 0
    for (int j = 1; j < 5; j++) {
        float s = clusters[j].isRealBall ? 0. : DEcor*scale;
        d = min( s + abs( sdistanceToBall( p, clusters[j] ) ) , d );
    }
    return d/scale;
}
#else
// Normal version
float dist2balls( vec3 p, float scale) {
    float d = abs(p.z); //Ball 0
    for (int j = 1; j < 5; j++)
        if (clusters[j].isRealBall)
            d = min(abs(sdistanceToBall(p, clusters[j])), d);
    return d/scale;
}
#endif

void EuclideanFold(inout vec3 p) {
    if( euclideanTriangleType == 0)
        fold236(p.xy);
    else if ( euclideanTriangleType == 1)
        fold244(p.xy);
    else
        fold333(p.xy);
}

float map(inout vec3 p, inout float scale, inout vec4 orb) {
    float d = dist2balls(p, scale);
    // bool cond = true;
    // K: Very strange behaviour : doing cond && tru_reflect() cuts things out.
    // Remove all those tests and just see if the orbit stops. Now DE is "much" better.
    for (int i = 0; i < Iterations ; i++) {
        vec3 ap = p;
        EuclideanFold(p);
        for (int k = 0; k < 2; k++)
            try_reflect(p, coclusters[k], scale, orb);
        if (all(not(bvec3(p-ap))))
            break;
        d = min( d, dist2balls(p, scale));
    }
#if 0
    // Slow version. DEs to balls are each iterations
    return d;
#else
    // Fast version. distance to balls is done only once.
    return dist2balls(p, scale);
#endif
}

float DE(vec3 p) {
    float DEfactor=1.;
    //orbiTrap is not initialized by default. That's what was producing the strange haze and colors.
    orbitTrap = vec4(1);
    float d = map(p, DEfactor, orbitTrap);

    //Call basic shape and scale its DE
    return 0.3 * d ;// Fudging here because fudge factor desn't seem to affect shadows :-/
}

#preset default
AutoFocus = false
FOV = 1
Eye = 4.93725,0.000805855,2.80231
Target = 4.98886,-6.4035,-4.86124
UpLock = false
Up = 0,0,1
AutoFocus = false
FocalPlane = 2
Aperture = 0
InFocusAWidth = 1
DofCorrect = true
ApertureNbrSides = 5
ApertureRot = 0
ApStarShaped = false
Gamma = 1.2
ToneMapping = 5
Exposure = 1
Brightness = 1
Contrast = 1
Saturation = 1
GaussianWeight = 1
AntiAliasScale = 1.5
Bloom = false
BloomIntensity = 0.733096
BloomPow = 5.22093
BloomTaps = 23
BloomStrong = 6.57901
DepthToAlpha = false
Detail = -4.2
RefineSteps = 5
FudgeFactor = 0.5
MaxRaySteps = 500
MaxDistance = 200
Dither = 0.5
NormalBackStep = 1
DetailAO = -1.96825
coneApertureAO = 0.56566
maxIterAO = 19
FudgeAO = 0.349943
AO_ambient = 1
AO_camlight = 0.97638
AO_pointlight = 0.449594
AoCorrect = 0
Specular = 0.08918
SpecularExp = 30.435
CamLight = 0.364706,0.364706,0.364706,0.28169
AmbiantLight = 0.709804,0.709804,0.709804,0.95652
Reflection = 0.192157,0.192157,0.192157
ReflectionsNumber = 1
SpotGlow = false
SpotLight = 1,1,1,0.4124
LightPos = 2.7888,2.2138,-2.061
LightSize = 0
LightFallOff = 0.2314
LightGlowRad = 0
LightGlowExp = 0
HardShadow = 1
ShadowSoft = 11.1712
ShadowBlur = 0
perf = false
SSS = false
sss1 = 0.1
sss2 = 0.5
BaseColor = 0.337255,0.337255,0.337255
OrbitStrength = 0.36364
X = 0.25098,0.505882,0.756863,1
Y = 0.392157,0.392157,0.584314,1
Z = 0.603922,0.164706,0.776471,1
R = 0.262745,0.482353,1,0.29412
BackgroundColor = 0.270588,0.403922,0.6
GradientBackground = 0
CycleColors = true
Cycles = 1.66106
EnableFloor = false
FloorNormal = 0,0,1
FloorHeight = 0
FloorColor = 1,1,1
HF_Fallof = 0.187344
HF_Const = 0.05333
HF_Intensity = 0
HF_Dir = 0,0,1
HF_Offset = 0
HF_Color = 0.564706,0.752941,0.878431,1
HF_Scatter = 10
HF_Anisotropy = 0.168627,0.168627,0.168627
HF_FogIter = 2
HF_CastShadow = true
EnCloudsDir = true NotLocked
CloudDir = 0.273574,-0.720605,-0.780451 NotLocked
CloudScale = 1 NotLocked
CloudFlatness = 0 NotLocked
CloudTops = 1 NotLocked
CloudBase = -1 NotLocked
CloudDensity = 0.484136 NotLocked
CloudRoughness = 1 NotLocked
CloudContrast = 1 NotLocked
CloudColor = 0.65,0.68,0.7 NotLocked
CloudColor2 = 0.07,0.17,0.24 NotLocked
SunLightColor = 0.7,0.5,0.3 NotLocked
Cloudvar1 = 0.99 NotLocked
Cloudvar2 = 1 NotLocked
CloudIter = 3 NotLocked
CloudBgMix = 1 NotLocked
euclideanTriangleType = 1
dihedralAngles0_234 = 3,2,7
dihedralAngles1_234 = 3,2,7
dihedralAngle0_1 = 4
isRealBall = 0,0,0,0
Iterations = 50
DEcor = 1
#endpreset


#preset 333-pseudoKleinian-like
AutoFocus = false
FOV = 0.7
Eye = 5.6304,-0.403237,0.692656
Target = 9.07877,-9.34463,-2.66041
UpLock = false
Up = 0,0,1
AutoFocus = false
FocalPlane = 1
Aperture = 0.01
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
BloomIntensity = 0.733096
BloomPow = 5.22093
BloomTaps = 23
BloomStrong = 6.57901
DepthToAlpha = true
Detail = -4
RefineSteps = 4
FudgeFactor = 1
MaxRaySteps = 300
MaxDistance = 300
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
SpecularExp = 178.26
CamLight = 1,0.945098,0.898039,0.28169
AmbiantLight = 1,0.972549,0.917647,0.95652
Reflection = 0.352941,0.352941,0.352941
ReflectionsNumber = 1
SpotGlow = false
SpotLight = 1,0.901961,0.827451,0.73438
LightPos = 10,0.6502,2.6654
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
BaseColor = 0.666667,0.666667,0.498039
OrbitStrength = 0
X = 0.666667,1,0,1
Y = 1,0.533333,0,1
Z = 0.603922,0.164706,0.776471,1
R = 0.262745,0.482353,1,0.29412
BackgroundColor = 0.270588,0.403922,0.6
GradientBackground = 0
CycleColors = true
Cycles = 0.1
EnableFloor = false
FloorNormal = 0,0,1
FloorHeight = -4.5349
FloorColor = 0.533333,0.533333,0.533333
HF_Fallof = 0.187344
HF_Const = 0.05333
HF_Intensity = 0
HF_Dir = 0,0,1
HF_Offset = 0
HF_Color = 1,1,1,0.83607
HF_Scatter = 10
HF_Anisotropy = 0.133333,0.00784314,0
HF_FogIter = 2
HF_CastShadow = true
EnCloudsDir = true NotLocked
CloudDir = 0.57357,-0.720605,0.78045 NotLocked
CloudScale = 4 NotLocked
CloudFlatness = 0 NotLocked
CloudTops = 1 NotLocked
CloudBase = -1 NotLocked
CloudDensity = 0.484136 NotLocked
CloudRoughness = 1 NotLocked
CloudContrast = 1 NotLocked
CloudColor = 0.65,0.68,0.7 NotLocked
CloudColor2 = 0.07,0.17,0.24 NotLocked
SunLightColor = 0.7,0.5,0.3 NotLocked
Cloudvar1 = 0.99 NotLocked
Cloudvar2 = 1 NotLocked
CloudIter = 3 NotLocked
CloudBgMix = 1 NotLocked
euclideanTriangleType = 2
dihedralAngles0_234 = 2,2,4
dihedralAngles1_234 = 2,2,1
dihedralAngle0_1 = 1
isRealBall = 1,0,0,0
Iterations = 246
#endpreset



#preset 333-237-236-7
AutoFocus = false
FOV = 0.685022
Eye = 4.93123,0.221804,0.977257
Target = 4.52934,-9.42281,-2.16962
UpLock = false
Up = -0.0105119,-0.280342,0.86054
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
BloomIntensity = 0.733096
BloomPow = 5.22093
BloomTaps = 23
BloomStrong = 6.57901
DepthToAlpha = true
Detail = -3.5
RefineSteps = 4
FudgeFactor = 1
MaxRaySteps = 300
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
SpecularExp = 178.26
CamLight = 1,0.945098,0.898039,0.28169
AmbiantLight = 1,0.972549,0.917647,0.95652
Reflection = 0.352941,0.352941,0.352941
ReflectionsNumber = 1
SpotGlow = false
SpotLight = 1,0.901961,0.827451,0.73438
LightPos = 2.7888,-0.6502,-1.6654
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
BaseColor = 0.776471,0.776471,0.776471
OrbitStrength = 0.875
X = 0,1,0.164706,1
Y = 1,0.533333,0,1
Z = 0.603922,0.164706,0.776471,1
R = 0.262745,0.482353,1,0.29412
BackgroundColor = 0.270588,0.403922,0.6
GradientBackground = 0
CycleColors = false
Cycles = 0.1
EnableFloor = false
FloorNormal = 0,0,1
FloorHeight = 0
FloorColor = 1,1,1
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
CloudDir = 0.273574,-0.720605,-0.780451 NotLocked
CloudScale = 1 NotLocked
CloudFlatness = 0 NotLocked
CloudTops = 1 NotLocked
CloudBase = -1 NotLocked
CloudDensity = 0.484136 NotLocked
CloudRoughness = 1 NotLocked
CloudContrast = 1 NotLocked
CloudColor = 0.65,0.68,0.7 NotLocked
CloudColor2 = 0.07,0.17,0.24 NotLocked
SunLightColor = 0.7,0.5,0.3 NotLocked
Cloudvar1 = 0.99 NotLocked
Cloudvar2 = 1 NotLocked
CloudIter = 3 NotLocked
CloudBgMix = 1 NotLocked
euclideanTriangleType = 2
dihedralAngles0_234 = 2,2,7
dihedralAngles1_234 = 2,2,6
dihedralAngle0_1 = 7
isRealBall = 0,0,0,0
Iterations = 250
#endpreset


#preset volcano-lake
FOV = 0.685022
Eye = 5.52166,-0.807811,2.62172
Target = 5.10595904,-0.127263602,2.01835614
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
Bloom = false Locked
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
DepthToAlpha = false
ShowDepth = false
DepthMagnitude = 1
Detail = -4
RefineSteps = 5
FudgeFactor = 0.57407
MaxRaySteps = 918
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
Specular = 0.1
SpecularExp = 84.21
CamLight = 0.996078,0.996078,0.996078,0.5
AmbiantLight = 0.972549,0.972549,0.972549,0.3
Reflection = 0.478431,0.478431,0.478431
ReflectionsNumber = 3
SpotGlow = false
SpotLight = 1,1,1,0.2
LightPos = 10,1.1476,10
LightSize = 0
LightFallOff = 0.72072
LightGlowRad = 0
LightGlowExp = 0
HardShadow = 1
ShadowSoft = 19.0118
ShadowBlur = 0
perf = false
SSS = false
sss1 = 0.1
sss2 = 0.5
BaseColor = 0.552941,0.552941,0.380392
OrbitStrength = 0.38384
X = 0.686275,0.721569,0.74902,1
Y = 0.670588,0.701961,0.729412,0
Z = 0.662745,0.694118,0.721569,0
R = 0.635294,0.662745,0.690196,1
BackgroundColor = 0.168627,0.258824,0.384314
GradientBackground = 0
CycleColors = true
Cycles = 1.21034487
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
dihedralAngles0_234 = 4,3,3
dihedralAngles1_234 = 2,2,4
dihedralAngle0_1 = 2
isRealBall = 0,0,0,0
Iterations = 50
DEcor = 1
#endpreset



#preset vinberg-convention
AutoFocus = false
FOV = 0.685022
Eye = 4.93123,-0.4218,1.0726
Target = 11.2567,-7.30351,-1.00865
UpLock = false
Up = 0,0,1
AutoFocus = false
FocalPlane = 1
Aperture = 0.004
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
BloomIntensity = 0.733096
BloomPow = 5.22093
BloomTaps = 23
BloomStrong = 6.57901
DepthToAlpha = true
Detail = -3.8
RefineSteps = 5
FudgeFactor = 1
MaxRaySteps = 300
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
CamLight = 1,1,1,1
AmbiantLight = 0.729412,0.729412,0.729412,0.95652
Reflection = 0.352941,0.352941,0.352941
ReflectionsNumber = 2
SpotGlow = false
SpotLight = 1,1,1,0.73438
LightPos = 0.7888,-5.6502,2.6654
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
BaseColor = 0.666667,0.666667,0.498039
OrbitStrength = 0.1
X = 0.666667,0.666667,0.498039,1
Y = 1,0.533333,0,0
Z = 0.603922,0.164706,0.776471,0
R = 1,1,1,1
BackgroundColor = 0,0,0
GradientBackground = 2.83785
CycleColors = true
Cycles = 0.1
EnableFloor = true
FloorNormal = 0,0,1
FloorHeight = 4.5
FloorColor = 0.564706,0.564706,0.564706
HF_Fallof = 0.4
HF_Const = 0.1
HF_Intensity = 0
HF_Dir = 0,0,1
HF_Offset = 0
HF_Color = 0.2,0.4,0.6,0.3
HF_Scatter = 10
HF_Anisotropy = 0.168627,0.168627,0.168627
HF_FogIter = 1
HF_CastShadow = true
EnCloudsDir = true NotLocked
CloudDir = 0.273574,-0.720605,-0.780451 NotLocked
CloudScale = 10 NotLocked
CloudFlatness = 0 NotLocked
CloudTops = 1 NotLocked
CloudBase = -3.9 NotLocked
CloudDensity = 0.484136 NotLocked
CloudRoughness = 1 NotLocked
CloudContrast = 1 NotLocked
CloudColor = 0.65,0.68,0.7 NotLocked
CloudColor2 = 0.07,0.17,0.24 NotLocked
SunLightColor = 0.7,0.5,0.3 NotLocked
Cloudvar1 = 0.99 NotLocked
Cloudvar2 = 1 NotLocked
CloudIter = 3 NotLocked
CloudBgMix = 1 NotLocked
euclideanTriangleType = 2
dihedralAngles0_234 = 2,2,7
dihedralAngles1_234 = 2,2,5
dihedralAngle0_1 = 1
isRealBall = 0,0,0,0
Iterations = 250
#endpreset


#preset Knighty-cave
AutoFocus = false
FOV = 0.9169
Eye = -5.86127,-1.34218,1.6648
Target = 2.45314,-2.28201,1.17346
UpLock = false
Up = 0,0,1
AutoFocus = false
FocalPlane = 1
Aperture = 0
InFocusAWidth = 1
DofCorrect = true
ApertureNbrSides = 5
ApertureRot = 0
ApStarShaped = false
Gamma = 2.2
ToneMapping = 5
Exposure = 1
Brightness = 1
Contrast = 1
Saturation = 1
GaussianWeight = 1
AntiAliasScale = 1.5
Bloom = false
BloomIntensity = 0.733096
BloomPow = 5.22093
BloomTaps = 20
BloomStrong = 6.57901
DepthToAlpha = false
Detail = -3.5
RefineSteps = 5 Locked
FudgeFactor = 0.6
MaxRaySteps = 672
MaxDistance = 100
Dither = 0.5
NormalBackStep = 10
DetailAO = -0.80458
coneApertureAO = 0.58333
maxIterAO = 15
FudgeAO = 0.349943
AO_ambient = 1
AO_camlight = 0
AO_pointlight = 0
AoCorrect = 0
Specular = 0.2
SpecularExp = 189.44
CamLight = 0.364706,0.364706,0.364706,0
AmbiantLight = 0.490196,0.490196,0.490196,1
Reflection = 0.615686,0.615686,0.615686
ReflectionsNumber = 2 Locked
SpotGlow = true
SpotLight = 1,0.992157,0.862745,10
LightPos = -1.4868,0.3206,4.519
LightSize = 0.00285
LightFallOff = 0.2314
LightGlowRad = 0
LightGlowExp = 0
HardShadow = 1
ShadowSoft = 0
ShadowBlur = 0
perf = false
SSS = false
sss1 = 0.1
sss2 = 0.5
BaseColor = 0.0470588,0.137255,0.168627
OrbitStrength = 0.62626
X = 0.756863,0.0980392,0.0980392,0.53542
Y = 0.0745098,0.654902,0.219608,0.35978
Z = 0.152941,0.341176,0.776471,0.8527
R = 0.603922,0.611765,0.105882,0.14204
BackgroundColor = 0,0,0
GradientBackground = 0
CycleColors = true
Cycles = 2.9
EnableFloor = false
FloorNormal = 0,0,1
FloorHeight = 0
FloorColor = 1,1,1
HF_Fallof = 0.12944
HF_Const = 0.16763
HF_Intensity = 0
HF_Dir = 0,0,1
HF_Offset = 0
HF_Color = 0.584314,0.584314,0.584314,0
HF_Scatter = 1
HF_Anisotropy = 0.164706,0.12549,0.105882
HF_FogIter = 4
HF_CastShadow = true
EnCloudsDir = true NotLocked
CloudDir = 0.273574,-0.720605,-0.780451 NotLocked
CloudScale = 1 NotLocked
CloudFlatness = 0 NotLocked
CloudTops = 1 NotLocked
CloudBase = -1 NotLocked
CloudDensity = 0.484136 NotLocked
CloudRoughness = 1 NotLocked
CloudContrast = 1 NotLocked
CloudColor = 0.65,0.68,0.7 NotLocked
CloudColor2 = 0.07,0.17,0.24 NotLocked
SunLightColor = 0.7,0.5,0.3 NotLocked
Cloudvar1 = 0.99 NotLocked
Cloudvar2 = 1 NotLocked
CloudIter = 3 NotLocked
CloudBgMix = 1 NotLocked
euclideanTriangleType = 0
dihedralAngles0_234 = 2,2,3
dihedralAngles1_234 = 2,2,7
dihedralAngle0_1 = 5
isRealBall = 0,0,0,0
Iterations = 50
DEcor = 1
#endpreset

#preset Golden
FOV = 0.7
Eye = 5.44058359,-0.63471323,0.212286673
Target = 6.38604793,-0.945990288,0.30822297
Up = -0.105484887,-0.016167548,0.987109317
EquiRectangular = false
FocalPlane = 1.46621
Aperture = 0.00619
InFocusAWidth = 0.35484
DofCorrect = true
ApertureNbrSides = 10
ApertureRot = 130.802
ApStarShaped = true
Gamma = 2.2
ToneMapping = 5
Exposure = 1
Brightness = 1
Contrast = 1
Saturation = 1
GaussianWeight = 1
AntiAliasScale = 1.5
Bloom = true Locked
BloomIntensity = 0.733096
BloomPow = 5.22093
BloomTaps = 23
BloomStrong = 6.57901
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
RefineSteps = 4
FudgeFactor = 1
MaxRaySteps = 500
MaxDistance = 200
Dither = 0.5
NormalBackStep = 1
DetailAO = -1.96825
coneApertureAO = 0.378985
maxIterAO = 19
FudgeAO = 0.349943
AO_ambient = 1
AO_camlight = 1.59217
AO_pointlight = 0.449594
AoCorrect = 0
Specular = 0.01176
SpecularExp = 27.275
CamLight = 1,1,1,0
AmbiantLight = 1,1,1,0.70588
Reflection = 0.352941,0.352941,0.352941
ReflectionsNumber = 3
SpotGlow = true
SpotLight = 1,0.972549,0.807843,10
LightPos = 7.3103452,4.5517244,0.5064
LightSize = 0.0297
LightFallOff = 0.49438
LightGlowRad = 0.3226
LightGlowExp = 0.9524
HardShadow = 1 Locked
ShadowSoft = 0
ShadowBlur = 0
perf = false
SSS = false
sss1 = 0.1
sss2 = 0.5
BaseColor = 0.666667,0.666667,0.498039
OrbitStrength = 0
X = 0.666667,1,0,1
Y = 1,0.533333,0,1
Z = 0.603922,0.164706,0.776471,1
R = 0.262745,0.482353,1,0.29412
BackgroundColor = 0.270588,0.403922,0.6
GradientBackground = 0
CycleColors = true
Cycles = 0.1
EnableFloor = false
FloorNormal = 0,0,1
FloorHeight = -4.5349
FloorColor = 0.533333,0.533333,0.533333
HF_Fallof = 1.62651
HF_Const = 0.08631
HF_Intensity = 0.15
HF_Dir = 0,0,1
HF_Offset = -3.27
HF_Color = 0.701961,0.803922,0.956863,0.24273
HF_Scatter = 1.577
HF_Anisotropy = 0,0,0
HF_FogIter = 2
HF_CastShadow = true
EnCloudsDir = true
CloudDir = 0.57357,-0.720605,0.78045
CloudScale = 1
CloudOffset = 0,0,0
CloudFlatness = 0
CloudTops = 1
CloudBase = -1
CloudDensity = 1
CloudRoughness = 1
CloudContrast = 1
CloudBrightness = 1
CloudColor = 0.65,0.68,0.7
CloudColor2 = 0.07,0.17,0.24
SunLightColor = 0.7,0.5,0.3
Cloudvar1 = 0.99
Cloudvar2 = 1
CloudIter = 3
CloudBgMix = 1
WindDir = 0,0,1
WindSpeed = 1
euclideanTriangleType = 2
dihedralAngles0_234 = 2,2,5
dihedralAngles1_234 = 3,2,1
dihedralAngle0_1 = 1
isRealBall = 0,0,0,1
Iterations = 220
DEcor = 1
#endpreset


#preset Golden-2
FOV = 0.7
Eye = 5.31741896,-1.01478511,0.192899929
Target = 6.29037962,-1.22863262,0.280174241
Up = 0,0,1
EquiRectangular = false
FocalPlane = 1.46621
Aperture = 0.00619
InFocusAWidth = 0.35484
DofCorrect = true
ApertureNbrSides = 10
ApertureRot = 130.802
ApStarShaped = true
Gamma = 2.2
ToneMapping = 5
Exposure = 1
Brightness = 1
Contrast = 1
Saturation = 1
GaussianWeight = 1
AntiAliasScale = 1.5
Bloom = true Locked
BloomIntensity = 0.733096
BloomPow = 5.22093
BloomTaps = 23
BloomStrong = 6.57901
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
RefineSteps = 4
FudgeFactor = 1
MaxRaySteps = 500
MaxDistance = 200
Dither = 0.5
NormalBackStep = 1
DetailAO = -1.96825
coneApertureAO = 0.378985
maxIterAO = 19
FudgeAO = 0.349943
AO_ambient = 1
AO_camlight = 1.59217
AO_pointlight = 0.449594
AoCorrect = 0
Specular = 0.01176
SpecularExp = 27.275
CamLight = 1,1,1,0
AmbiantLight = 1,1,1,0.70588
Reflection = 0.352941,0.352941,0.352941
ReflectionsNumber = 3
SpotGlow = true
SpotLight = 1,0.972549,0.807843,10
LightPos = 7.3103452,4.5517244,0.5064
LightSize = 0.0297
LightFallOff = 0.49438
LightGlowRad = 0.3226
LightGlowExp = 0.9524
HardShadow = 1 Locked
ShadowSoft = 0
ShadowBlur = 0
perf = false
SSS = false
sss1 = 0.1
sss2 = 0.5
BaseColor = 0.666667,0.666667,0.498039
OrbitStrength = 0
X = 0.666667,1,0,1
Y = 1,0.533333,0,1
Z = 0.603922,0.164706,0.776471,1
R = 0.262745,0.482353,1,0.29412
BackgroundColor = 0.270588,0.403922,0.6
GradientBackground = 0
CycleColors = true
Cycles = 0.1
EnableFloor = false
FloorNormal = 0,0,1
FloorHeight = -4.5349
FloorColor = 0.533333,0.533333,0.533333
HF_Fallof = 1.62651
HF_Const = 0.08631
HF_Intensity = 0.15
HF_Dir = 0,0,1
HF_Offset = -3.27
HF_Color = 0.701961,0.803922,0.956863,0.24273
HF_Scatter = 1.577
HF_Anisotropy = 0,0,0
HF_FogIter = 2
HF_CastShadow = true
EnCloudsDir = true
CloudDir = 0.57357,-0.720605,0.78045
CloudScale = 1
CloudOffset = 0,0,0
CloudFlatness = 0
CloudTops = 1
CloudBase = -1
CloudDensity = 1
CloudRoughness = 1
CloudContrast = 1
CloudBrightness = 1
CloudColor = 0.65,0.68,0.7
CloudColor2 = 0.07,0.17,0.24
SunLightColor = 0.7,0.5,0.3
Cloudvar1 = 0.99
Cloudvar2 = 1
CloudIter = 3
CloudBgMix = 1
WindDir = 0,0,1
WindSpeed = 1
euclideanTriangleType = 2
dihedralAngles0_234 = 2,2,5
dihedralAngles1_234 = 3,2,1
dihedralAngle0_1 = 1
isRealBall = 0,0,0,1
Iterations = 220
DEcor = 1
#endpreset
