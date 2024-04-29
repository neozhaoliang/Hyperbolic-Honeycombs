#version 330
#define providesInit
#define providesColor
#define KN_VOLUMETRIC
#define USE_EIFFIE_SHADOW
#define MULTI_SAMPLE_AO
#include "MathUtils.frag"
#include "DE-Kn2cr11.frag"

uniform vec3 Eye;
uniform vec3 Target;

#group Object-Colors
uniform vec3 segAColor; color[0.0,0.0,0.0]
uniform vec3 segBColor; color[0.0,0.0,0.0]
uniform vec3 segCColor; color[0.0,0.0,0.0]
uniform vec3 segDColor; color[0.0,0.0,0.0]
uniform vec3 faceABCol; color[0.0,0.0,0.0]
uniform vec3 faceACCol; color[0.0,0.0,0.0]
uniform vec3 faceADCol; color[0.0,0.0,0.0]
uniform vec3 faceBCCol; color[0.0,0.0,0.0]
uniform vec3 faceBDCol; color[0.0,0.0,0.0]
uniform vec3 faceCDCol; color[0.0,0.0,0.0]
uniform vec3 verticesColor; color[0.0,0.0,0.0]
uniform vec3 facesColor; color[0.0,0.0,0.0]
uniform vec3 floorCheckerColor1; color[0.82,0.196078,0.33]
uniform vec3 floorCheckerColor2; color[0.0,0.0,0.0]
uniform vec3 floorLineColor; color[0.0,0.0,0.0]

#group Honeycomb-Settings
// draw the honeycomb in upper half space model
uniform bool doUpperHalfSpace; checkbox[false]
// this option should only be true for hyperideal cases
uniform bool drawFloorPattern; checkbox[false]
uniform int AB; slider[2,5,9]
uniform int AC; slider[2,2,9]
uniform int AD; slider[2,2,9]
uniform int BC; slider[2,3,9]
uniform int BD; slider[2,2,9]
uniform int CD; slider[2,4,9]
uniform bool enableFaceAB; checkbox[false]
uniform bool enableFaceAC; checkbox[false]
uniform bool enableFaceAD; checkbox[false]
uniform bool enableFaceBC; checkbox[false]
uniform bool enableFaceBD; checkbox[false]
uniform bool enableFaceCD; checkbox[false]
// Barycentric coordinates of the initial vertex v0
uniform vec4 activeMirrors; slider[(0,0,0,0),(1,0,0,0),(1,1,1,1)]
uniform float vertexSize; slider[0,0.05,0.5]
uniform float edgeSize; slider[0,0.01,0.1]
uniform float faceThickness; slider[0,0.001,0.2]
uniform float floorLineThickness; slider[0,0.01,3.]
uniform int Iterations; slider[0,30,1000]
uniform float CSphRad;  slider[0,0.99995,1]

// normal vectors of the four reflection mirrors
// the reflections are in the 4d Minkowski space
mat4 M;
// the initial vertex
vec4 v0;

// This matrix also holds the reflection mirrors,
// except they are planes and spheres in 3d space.
// This is use for computing the floor pattern.
// The floor is the ideal boundary of the hyperbolic model hence cannot
// be lifted to the hyperboloid, but we can still reflect and inversion them
// in 3d space.
// The first three mirrors re normals of 3d planes,
// The last one contains the center and radius of the inversion sphere.
mat4 M2;

// cosh, sinh of the vertex radius and edge radius
float cvr, svr, csr, ssr;

// Minkowski inner product with Sylvester type (3, 1)
float hdot(in vec4 p, in vec4 q) {
    return dot(p.xyz, q.xyz) - p.w * q.w;
}

// normalize a time-like vector <v, v> < 0
vec4 hnormalize(in vec4 p) {
    return p / sqrt(-hdot(p, p));
}

// reflection about a plane in 4d hyperbolic space
float try_reflect(inout vec4 p, in vec4 n) {
    float k = min(0., hdot(p, n));
    p -= 2. * k * n;
    return k;
}

// reflection about a plane in 3d space, this is the usual Euclidean reflection
bool try_reflect(inout vec3 p, vec3 n, inout int count) {
    float k = dot(p, n);
    if (k >= 0.)
        return true;
    p -= 2.0 * k * n;
    count += 1;
    return false;
}

// Sphere inversion in 3d space
bool try_reflect(inout vec3 p, vec4 sphere, inout int count) {
    vec3 cen = sphere.xyz;
    float r = sphere.w;
    vec3 q = p - cen;
    float d2 = dot(q, q);
    if (d2 == 0.0)
    	return true;
    float k = (r * r) / d2;
    if (k < 1.0)
    	return true;
    p = k * q + cen;
    count += 1;
    return false;
}

// Minimal Euclidean distance to the four mirrors in M2.
float distABCD(vec3 p) {
    float dA = abs(dot(p, M2[0].xyz));
    float dB = abs(dot(p, M2[1].xyz));
    float dC = abs(dot(p, M2[2].xyz));
    float dD = abs(length(p - M2[3].xyz) - M2[3].w);
    return min(dA, min(dB, min(dC, dD)));
}

// Inverse stereo-graphic projection, from a point on plane z=0 to
// the unit ball centered at the origin
vec3 planeToSphere(vec2 p) {
    float r2 = dot(p, p);
    return vec3(2.0 * p, r2 - 1.0) / (1.0 + r2);
}

// get the color for the checker pattern on the floor
vec3 getFloorColor(bool found, int count) {
    if (found) {
        return (count % 2 == 0) ? floorCheckerColor1 : floorCheckerColor2;
    }
    return floorLineColor;
}

// initialize the data of the honeycomb
void init() {
    float c01 = -cos(PI / float(AB));
    float c02 = -cos(PI / float(AC));
    float c03 = -cos(PI / float(AD));
    float c12 = -cos(PI / float(BC));
    float c13 = -cos(PI / float(BD));
    float c23 = -cos(PI / float(CD));

    vec4 A, B, C, D;
    // find the reflection mirrors A, B, C, D.
    // A can be always chosen as x-axis
    A = vec4(1, 0, 0, 0);
    B = vec4(c01, sqrt(1. - c01*c01), 0., 0.);
    C = vec4(c02, 0, 0, 0);
    C.y = (c12 - C.x * B.x) / B.y;
    C.z = sqrt(abs(1. - dot(C.xy, C.xy))); // avoid rounding error in paracompact case

    D = vec4(c03, 0, 0, 0);
    D.y = (c13 - D.x * B.x) / B.y;
    D.z = (c23 - dot(D.xy, C.xy) ) / C.z;
    // !important: if you want to make the fundamental chamber lie in the upper
    // sheet of the hyperboloid then you must use "-" sign for the last entry of D.
    D.w = -sqrt(abs(dot(D.xyz, D.xyz) - 1.));

    vec4 H = vec4(1, 1, 1, -1);
    mat4 Minv = inverse(mat4(H*A, H*B, H*C, H*D));
    v0 = hnormalize(activeMirrors * Minv);
    M = mat4(A, B, C, D);

    // cosh(vradius), sinh(vradius), cosh(sradius), sinh(sradius)
    cvr = cosh(vertexSize); svr = sinh(vertexSize);
    csr = cosh(edgeSize); ssr = sinh(edgeSize);

    // M2 and M differ only in their last column vector
	vec3 cen = -vec3(c03, c13, c23) *  inverse(mat3(A.xyz, B.xyz, C.xyz)) ;
    vec4 S = vec4(cen, 1);
    S /= sqrt(dot(cen, cen) - 1.);
    M2 = mat4(A, B, C, S);
}

// for a 3d point on the unit sphere in Euclidean space, reflect it until it falls into
// the 3d fundamental domain.
bool fold3d(inout vec3 p, inout int count) {
    bool inA, inB, inC, inD;
    for (int iter = 0; iter < 300; iter++) {
        inA = try_reflect(p, M2[0].xyz, count);
        inB = try_reflect(p, M2[1].xyz, count);
        inC = try_reflect(p, M2[2].xyz, count);
        inD = try_reflect(p, M2[3], count);
        if (inA && inB && inC && inD)
            return true;
    }
    return false;
}

// for a 4d point in hyperbolic space, reflect it until it falls into the
// 4d fundamental domain
bool fold4d(inout vec4 p) {
    float k;
    for(int i = 0; i < Iterations; i++) {
        k = 0.;
        p.x = abs(p.x);
        k += try_reflect(p, M[1]);
        k += try_reflect(p, M[2]);
        k += try_reflect(p, M[3]);
        // break as soon as we find it's already in the fundamental domain
        if(k == 0.) return true;
    }
    return false;
}

// knighty's conservative distance conversion from hyperboloid to 3d flat distance.
// for a 3d point p, it's lifted to a 4d point q on the hyperboloid:
//
//        2p     1+r^2
// q = ( -----,  ----- ).
//       1-r^2   1-r^2
//
// any point with distance d (ca=cosh(d),sa=sinh(d)) to q can be written as q*ca + v*sa,
// wheren v is a unit tangent vector at q. We want v look like (sp, t).
// so we have two unknowns (s, t) and two equations:

// 1. s^2 * r^2 - t^2 = 1  (tangent vector must be space-like)
// 2. 2*s/(1-r^2)*r^2 - (1+r^2)/(1-r^2)*t = 0 (definition of tangent vector at q)

// solve for s, t we have s = (1+r^2)/(r *(1-r^2)) and t = 2*r/(1-r^2).
// hence the point we choose to project to 3d is
// ([ 2*ca/(1-r^2) + (1+r^2)/(r *(1-r^2))*sa ] * p,
//  [ ca*(1+r^2)/(1-r^2) + sa*2*r/(1-r^2) ])
float knightyDD(float ca, float sa, float r) {
    float x = 1. + r * r;
    float y = 2. - x;
    return (2. * r * ca + x * sa) / (x * ca + 2. * r * sa + y) - r;
}

// if distance between p and q is a, C is a circle with radius VR centered at q,
// then the distance from p to C is a - VR, hence
// cosh(a - VR) = cosh(a)cosh(VR) - sinh(a)sinh(VR)
// sinh(a - VR) = sinh(a)cosh(VR) - cosh(a)sinh(VR)
float dVertex(vec4 p, float r) {
    float ca = -hdot(p, v0);
    float sa = 0.5 * sqrt(-hdot(p - v0, p - v0) * hdot(p + v0, p + v0));
    return knightyDD(ca * cvr - sa * svr,
                     sa * cvr - ca * svr, r);
}

// let pj = a * n + b * v0 be the projection of p onto the plane given by (v0, n).
// take Minkowski dot with p - pj by v0 and n:
// (p, n) = a + b * (v0, n)
// (p, v0) = a * (v0, n) - b
// then solve this 2x2 linear system.
float dSegment(vec4 p, vec4 n, float r) {
    float pn = hdot(p, n);
    float pv = hdot(p, v0);
    float nv = hdot(n, v0);
    float det = -1.0 - nv * nv;
    float a = (-nv * pv - pn) / det;
    float b = (pv - pn * nv) / det;
    vec4 pj = hnormalize(min(a, 0.) * n + b * v0);
    float ca = -hdot(p, pj);
    float sa = 0.5 * sqrt(-hdot(p - pj, p - pj) * hdot(p + pj, p + pj));
    return knightyDD(ca * csr - sa * ssr, sa * csr - ca * ssr, r);
}


float dFace(vec4 p, float r, vec4 n1, vec4 n2) {
    mat3 m = mat3(-1, hdot(v0, n1), hdot(v0, n2),
                  hdot(v0, n1), 1, hdot(n1, n2),
                  hdot(v0, n2), hdot(n1, n2), 1);
    vec3 b = vec3(hdot(p, v0), hdot(p, n1), hdot(p, n2));
    vec3 c = b * inverse(m);
    if (c.y > 0. || c.z > 0.)
	return 1e5;
    vec4 q = c.x * v0 + c.y * n1 + c.z * n2;
    q = hnormalize(q);
    float ca = -hdot(p, q);
    float sa = 0.5 * sqrt(-hdot(p - q, p - q) * hdot(p + q, p + q));
    float ct = cosh(faceThickness);
    float st = sinh(faceThickness);
    return knightyDD(ca * ct - sa * st, sa * ct - ca * st, r);
}

float dSegments(vec4 p, float r) {
    float dA = dSegment(p, M[0], r);
    float dB = dSegment(p, M[1], r);
    float dC = dSegment(p, M[2], r);
    float dD = dSegment(p, M[3], r);
    return min(min(dA, dB), min(dC, dD));
}

float dFaces(vec4 p, float r) {
    float df = 1e5;
    if (enableFaceAB)
        df = min(df, dFace(p, r, M[0], M[1]));
    if (enableFaceAC)
        df = min(df, dFace(p, r, M[0], M[2]));
    if (enableFaceAD)
        df = min(df, dFace(p, r, M[0], M[3]));
    if (enableFaceBC)
        df = min(df, dFace(p, r, M[1], M[2]));
    if (enableFaceBD)
        df = min(df, dFace(p, r, M[1], M[3]));
    if (enableFaceCD)
        df = min(df, dFace(p, r, M[2], M[3]));
    return df;
 }


float DE(vec3 p) {
    // save the distanc from p to the plane z=0 for drawing floor pattern.
    // due to floating error we cannot do this exactly, so move the plane up a bit
    float h = p.z - 1e-3;
    if (doUpperHalfSpace) {
        p.z += 1.0;
        p *= 2.0 / dot(p, p);
        p.z -= 1.0;
    }
    float r = length(p);
    vec4 q = vec4(2.*p, 1.+r*r) / (1.-r*r);
    bool found = fold4d(q);
    orbitTrap = q;
    float dV = dVertex(q, r);
    float dS = dSegments(q, r);
    float dF = dFaces(q, r);
    if (doUpperHalfSpace)
        return min(dF, min(h, min(dV, dS)));
    else
        return max(r - CSphRad, min(dF, min(dV, dS)));
}

vec3 baseColor(vec3 pos, vec3 normal) {
    vec3 p0 = pos;  // save pos for later use
    float h = pos.z - 1e-3;
    if (doUpperHalfSpace) {
        pos.z += 1.;
        pos *= 2. / dot(pos, pos);
        pos.z -= 1.;
    }
    float r = length(pos);
    vec4 q = vec4(2.*pos, 1.+r*r) / (1.-r*r);
    bool found = fold4d(q);
    if (found) {
        float dA = dSegment(q, M[0], r);
        float dB = dSegment(q, M[1], r);
        float dC = dSegment(q, M[2], r);
        float dD = dSegment(q, M[3], r);
        float dV = dVertex(q, r);
        vec3 fcol;
        float dF = 1e5, dAB = 1e5, dAC = 1e5, dAD = 1e5, dBC = 1e5, dBD = 1e5, dCD = 1e5;
        if (enableFaceAB)
            dAB = dFace(q, r, M[0], M[1]);
        if (enableFaceAC)
            dAC = dFace(q, r, M[0], M[2]);
        if (enableFaceAD)
            dAD = dFace(q, r, M[0], M[3]);
        if (enableFaceBC)
            dBC = dFace(q, r, M[1], M[2]);
        if (enableFaceBD)
            dBD = dFace(q, r, M[1], M[3]);
        if (enableFaceCD)
            dCD = dFace(q, r, M[2], M[3]);
        dF = min(dAB, min(dAC, min(dAD, min(dBC, min(dBD, dCD)))));
        if (dF == dAB) fcol = faceABCol;
        if (dF == dAC) fcol = faceACCol;
        if (dF == dAD) fcol = faceADCol;
        if (dF == dBC) fcol = faceBCCol;
        if (dF == dBD) fcol = faceBDCol;
        if (dF == dCD) fcol = faceCDCol;
        float d = min(dF, min(min(min(dA, dB), min(dC, dD)), dV));
        if (doUpperHalfSpace)
            d = min(d, h);

        vec3 color = segAColor;
        if (d == dB) color = segBColor;
        if (d == dC) color = segCColor;
        if (d == dD) color = segDColor;
        if (d == dV) color = verticesColor;
        if (d == dF) color = fcol;
        if (doUpperHalfSpace && d == h) {
            if (drawFloorPattern) {
                int count = 0;
                vec3 q = planeToSphere(p0.xy);
                bool found =fold3d(q, count);
                color = getFloorColor(found, count);
                float edist = distABCD(q);
                color = mix(color, floorLineColor, (1.0 - smoothstep(0., 0.005, edist))*floorLineThickness);
            }
            else
                color = BackgroundColor;
        }
        return color;
    }

    return BackgroundColor;
}
































#preset Default
FOV = 0.71702
Eye = 2.58896218,-0.040928115,0.959116097
Target = 1.68544603,0.100495596,0.554569566
Up = 0,0,1
EquiRectangular = false
FocalPlane = 1
Aperture = 0
InFocusAWidth = 0
DofCorrect = true
ApertureNbrSides = 5
ApertureRot = 0.0108
ApStarShaped = false
Gamma = 2.2
ToneMapping = 4
Exposure = 1
Brightness = 1
Contrast = 1
Saturation = 1
GaussianWeight = 1
AntiAliasScale = 2
Bloom = false
BloomIntensity = 0.25
BloomPow = 2
BloomTaps = 4
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
Detail = -3.53
RefineSteps = 4
FudgeFactor = 1
MaxRaySteps = 2000
MaxDistance = 200
Dither = 0
NormalBackStep = 3.0357
DetailAO = -2.00788
coneApertureAO = 1
maxIterAO = 5
FudgeAO = 1
AO_ambient = 0.9399
AO_camlight = 0
AO_pointlight = 0
AoCorrect = 0
Specular = 0.3743
SpecularExp = 257.06
CamLight = 0.878431,0.921569,0.956863,0
AmbiantLight = 1,1,1,0.76
Reflection = 0.517647059,0.517647059,0.517647059
ReflectionsNumber = 3 Locked
SpotGlow = false
SpotLight = 0.94902,0.87451,0.721569,6
LightPos = 0.711744,0.9964416,1.3
LightSize = 0
LightFallOff = 1.09434
LightGlowRad = 0
LightGlowExp = 1
HardShadow = 1
ShadowSoft = 0
ShadowBlur = 0.03
perf = false
SSS = false
sss1 = 0.1
sss2 = 0.5
BaseColor = 0.835294,0.835294,0.835294
OrbitStrength = 0.66295
X = 0.0627451,0.0784314,0.0784314,0.69038
Y = 0.72549,0.85098,1,0.3923
Z = 0.8,0.78,1,-0.50384
R = 0.4,0.7,1,0.46872
BackgroundColor = 0.337254902,0.447058824,0.560784314
GradientBackground = 0.6
CycleColors = false
Cycles = 0.65577
EnableFloor = false
FloorNormal = 0,0,1
FloorHeight = 0
FloorColor = 1,1,1
HF_Fallof = 1.08964
HF_Const = 0.05181347
HF_Intensity = 0.1
HF_Dir = 0,0,1
HF_Offset = -2.4384
HF_Color = 0.592157,0.737255,1,1
HF_Scatter = 0.1
HF_Anisotropy = 0,0,0
HF_FogIter = 1
HF_CastShadow = true
EnCloudsDir = false
CloudDir = 0,0,1
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
CloudIter = 5
CloudBgMix = 1
WindDir = 0,0,1
WindSpeed = 1
segAColor = 0.666667,0.333333,1
segBColor = 0.403922,0.552941,0.541176
segCColor = 0.509803922,0.678431373,0.0549019608
segDColor = 0.8,0.74902,0.847059
faceABCol = 0.639216,0.639216,0.317647
faceACCol = 0.250980392,0.392156863,0.211764706
faceADCol = 0,0.152941,0.462745
faceBCCol = 0.529411765,0.678431373,1
faceBDCol = 0.811764706,0.37254902,0.607843137
faceCDCol = 0.333333333,0.333333333,1
verticesColor = 0.152941176,0.282352941,0.756862745
facesColor = 0,0,0
floorCheckerColor1 = 0.82,0.196078,0.33
floorCheckerColor2 = 0.196078,0.35,0.92
floorLineColor = 0,0,0
doUpperHalfSpace = true
drawFloorPattern = true
AB = 6
AC = 2
AD = 2
BC = 2
BD = 3
CD = 5
enableFaceAB = false
enableFaceAC = false
enableFaceAD = false
enableFaceBC = false
enableFaceBD = true
enableFaceCD = false
activeMirrors = 0,0,1,0.10623557
vertexSize = 0.20620843
edgeSize = 0.04096
faceThickness = 0.005
floorLineThickness = 1.85
Iterations = 150
CSphRad = 1
#endpreset

