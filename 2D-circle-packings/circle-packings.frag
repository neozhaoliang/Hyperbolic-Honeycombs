#version 330
#define providesColor
#define providesInit
#define MULTI_SAMPLE_AO
#include "Complex.frag"

#include "Soft-Raytracer.frag"

#group BallPacking-Settings
uniform int  Iterations; slider[1,200,1000]
uniform float orbitDivisor; slider[0,300,2000]
// Move a vertex of Euclidean type (if there is any) to infinity,
// so the pattern tiles the entire plane
uniform bool moveVertexToInf; checkbox[false]
// The label of edge CD must be finite
uniform vec3 dihedral_A_BCD; slider[(1,1,2),(4,4,4),(10,10,10)]
// The triangle BCD must be hyperbolic
uniform vec3 TriangleBCD; slider[(1,1,1),(4,4,4),(10,10,10)]
uniform float edgeSize; slider[0.,0.005,0.05]
uniform float edgeSize2; slider[0.,0.015,0.05]

uniform bool doInvert; checkbox[true];
uniform vec2 invPoint; slider[(0.5,0.5),(0.75,0.75),(0.9,0.9)]


#define inf        1.0
#define L2(x)               dot(x, x)
#define L2XY(x, y)          L2(x - y)
#define Hyperbolic          -1.0
#define Euclidean           0.
#define Spherical           1.

// 4x4 Cartan matrix for the Coxeter group
mat4 cartan;

// geometry type of each vertex in the Coxeter diagram
// -1 for hyerbolic, 0 for Euclidean, 1 for spherial
// To be determined from the Cartan matrix
vec4 vertexType = vec4(1);

// coordinates of a Euclidean vertex
vec2 euclideanVertex;

bool hasEuclideanVertex;


// compute cos(PI / x), for x = infiniy this is inf using Vinberg's notation
float dihedral(float x) {
    return x == inf ? inf : cos(PI / x);
}


// compute the vertex type of each vertex.
// for each vertex in the Coxeter diagram, its vertex type is determined by the triangle group G formed by
// the remaining three vertices. This vertex is of hyperbolic/spherical/euclidean iff G is
// hyperbolic/spherical/euclidean, respectively.
// This can be checked from the determinant of the leading minors of the Cartan matrix.
void checkCartan(mat4 M) {
    const float e = 0.001;
    float det;
    det = determinant(mat3(M[1].yzw, M[2].yzw, M[3].yzw));
    vertexType.x = step(-e, det) + step(e, det) - 1.;

    det = determinant(mat3(M[0].xzw, M[2].xzw, M[3].xzw));
    vertexType.y = step(-e, det) + step(e, det) - 1.;

    det = determinant(mat3(M[0].xyw, M[1].xyw, M[3].xyw));
    vertexType.z = step(-e, det) + step(e, det) - 1.;

    det = determinant(mat3(M[0].xyz, M[1].xyz, M[2].xyz));
    vertexType.w = step(-e, det) + step(e, det) - 1.;
}

// For spheres cen is the center, r is the radius
// For planes cen is the normal vector, r is the offset from the origin along the normal.
// if invert is true then the inside/outside of the sphere is exchanged.
struct Ball {
    vec2 cen;
    float r;
    bool isplane;
    bool invert;
};

// coclusters are mirror balls, they corresponde to root vectors (space-like)
Ball[4] coclusters;
// clusters are real balls, they corresponde to space-like weight vectors
Ball[4] clusters;

// create a default ball (unit circle)
Ball defaultBall() {
    return Ball(vec2(0), 1., false, false);
}

Ball from_plane(vec2 normal, float offset) {
    return Ball(normal, offset, true, false);
}

Ball from_sphere(vec2 cen, float r) {
    return Ball(cen, r, false, false);
}

void invertBall(inout Ball B) {
    B.invert = !B.invert;
}


// try to reflect a point p to the positive half space bounded by a ball
// if we are already in the positive half space, do nothing and return true,
// else reflect about the ball and return false
// if B is a sphere we try to reflect p into the interior of B
bool try_reflect(inout vec2 p,
                 Ball B,
                 inout int count,
                 inout float scale,
                 inout vec4 orb) {
    orb = min(orb, vec4(abs(p), dot(p,p), 1000.));
    vec2 cen = B.cen;
    float r = B.r;
    if (B.isplane) {
        float k = dot(vec3(p, 1), vec3(cen, r));
        if (k >= 0.)
            return true;
        p -= 2. * k  * cen;
        count += 1;
        return false;
    }
    else {
        vec2 q = p - cen;
        float d2 = dot(q, q);
        float k = (r * r) / d2;
        if ( (k < 1.0 && B.invert) || (k > 1. && !B.invert) )
            return true;
        p = k * q + cen;
        scale *= k;
        count += 1;
        return false;
    }
}

vec2 getIntersection(Ball B1, Ball B2, Ball B3) {
    vec2 dir = vec2(-B3.cen.y, B3.cen.x);
    float r1 = B1.r, r2 = B2.r;
    float k = (L2(B1.cen) - L2(B2.cen) - (r1*r1 - r2*r2)) / (2. * dot(B1.cen - B2.cen, dir));
    return k*dir;
}

Ball solveBall(Ball B1, Ball B2, Ball B3) {
    vec2 dir = vec2(-B3.cen.y, B3.cen.x);
    float r1 = B1.r, r2 = B2.r;
    float k = (L2(B1.cen) - L2(B2.cen) - (r1*r1 - r2*r2)) / (2. * dot(B1.cen - B2.cen, dir));
    vec2 cen = k * dir;
    float r = sqrt(L2XY(cen, B1.cen) - r1*r1);
    return from_sphere(cen, r);
}

float sdistanceToBall(vec2 p, Ball B) {
    if (B.isplane) {
        float k = dot(vec3(p, 1), vec3(B.cen, B.r));
        return k;
    }
    else {
        float k = length(p - B.cen) - B.r;
        return B.invert ? -k : k;
    }
}


void init() {
    Ball B0, B1, B2, B3;
    float c01 = dihedral(dihedral_A_BCD.x);
    float c02 = dihedral(dihedral_A_BCD.y);
    float c03 = dihedral(dihedral_A_BCD.z);
    float c12 = dihedral(TriangleBCD.x);
    float c13 = dihedral(TriangleBCD.y);
    float c23 = dihedral(TriangleBCD.z);

    cartan = mat4(1, -c01, -c02, -c03,
                  -c01, 1, -c12, -c13,
                  -c02, -c12, 1, -c23,
                  -c03, -c13, -c23, 1);
    checkCartan(cartan);

    float s23 = sqrt(1. - c23*c23);

    // The two virtuak balls B2, B3 are lines through the origin
    B2 = from_plane(vec2(1, 0), 0.);
    B3 = from_plane(vec2(-c23, s23), 0.);

    float k1 = c12;
    float k2 = (c13 + c23*c12) / s23;
    float r = 1. / sqrt(k1*k1 + k2*k2 - 1.);

    B1 = from_sphere(vec2(k1*r, k2*r), r);

    k1 = c02;
    k2 = (c03 + c23*c02) / s23;

    float a = k1*k1 + k2*k2 - 1.;
    float b = dot(vec3(k1, k2, c01), vec3(B1.cen, B1.r));
    float c = L2(B1.cen) - B1.r*B1.r;

    r = b / a - sqrt(b*b - a*c) / a;
    B0 = from_sphere(vec2(k1*r, k2*r), r);

    invertBall(B0);
    invertBall(B1);

    coclusters = Ball[4] (B0, B1, B2, B3);

    for (int k = 0; k < 4; k++) {
        clusters[k] = defaultBall();
    }

    invertBall(clusters[0]);

    if (vertexType.y == Hyperbolic) {
        float r = sqrt(L2(B0.cen) - B0.r*B0.r);
        clusters[1] = from_sphere(vec2(0), r);
    }
    if (vertexType.z == Hyperbolic) {
        clusters[2] = solveBall(B0, B1, B3);
    }
    if (vertexType.w == Hyperbolic) {
        clusters[3] = solveBall(B0, B1, B2);
    }

    if (vertexType.y == Euclidean) {
        hasEuclideanVertex = true;
        euclideanVertex = vec2(0);
        return;
    }

    if (vertexType.z == Euclidean) {
        hasEuclideanVertex = true;
        euclideanVertex = getIntersection(B0, B1, B3);
        return;
    }

    if (vertexType.w == Euclidean) {
        hasEuclideanVertex = true;
        euclideanVertex = getIntersection(B0, B1, B2);
        return;
    }
}

vec2 applyMobius(vec2 p) {
    if (hasEuclideanVertex) {
        vec2 A = euclideanVertex;
        vec2 B = vec2(0, 0);
        vec2 C = vec2(1, 0);
        vec2 D = vec2(4, 0);
        p = cDiv(cMul(p, A) + B, cMul(C, p) + D);
    }
    return p;
}

bool outside = false;

float distanceToMirrors(vec2 p) {
    float d = abs(sdistanceToBall(p, coclusters[0]));
   if (length(p)>1.) {
        p /= dot(p,p);
        outside = true;
	}
    for (int k = 1; k < 4; k++) {
        d = min(d, abs(sdistanceToBall(p, coclusters[k])));
    }
    return d;
}


void fold(inout vec2 p,
          inout int count,
          inout float scale,
          inout int index,
          inout vec4 orb) {
    if (moveVertexToInf)
        p = applyMobius(p);

    for (int i = 0; i < Iterations; i++) {
        vec2 ap = p;
        for (int k = 0; k < 4; k++) {
            try_reflect(p, coclusters[k], count, scale, orb);
        }
        if (all(not(bvec2(p - ap))))
            break;
    }
    for (int k = 0; k < 4; k++) {
        if (vertexType[k] == Hyperbolic && sdistanceToBall(p, clusters[k]) < -0.0001) {
            index = k;
            break;
        }
    }

}



// signed distance to unit ball and plane z=-1
float sdSphere(vec3 p) { return length(p) - 1.0; }
float sdPlane(vec3 p) { return p.z + 1.0; }

// project points on the unit ball to plabe z=-1
vec2 sphereToPlane(vec3 p) {
    return 2. * p.xy / (1. - p.z);
}

// you can implement your color functions here
vec3 colormap(int index, float t) {
    float c = float(index) + 1.;
    return .5 + .45*cos(2.*PI  * pow(t, 0.4) * c + vec3(c, c+2., c+1.) / .6 + vec3(0, 1, 2));
}

vec3 getColor(vec2 p) {
    orbitTrap = vec4(1);
    int index = -1;
    int count = 0;
    float scale = 1.;
    fold(p, count, scale, index, orbitTrap);
    float dist = distanceToMirrors(p);
    float t = clamp(float(count + 1) / orbitDivisor, 0., 1.);
    vec3 col = colormap(index, t);
    float ew = outside ? edgeSize2 : edgeSize;
    if (dist < ew)
        col = vec3(0);
    return col;
}

bool isFloor = false;

vec3 baseColor(vec3 p, vec3 n) {
	isFloor = false;
	DE(p);
	vec2 z = p.xy;
	if (!isFloor)
	    z = sphereToPlane(p).xy;
	if (doInvert) {
        float k = 1. / L2(invPoint);
        vec2 invCtr = k * invPoint;
        float t = (k - 1.) / L2(z -invCtr);
        z = t*z + (1. - t)*invCtr;
    }
	return getColor(z);
}


float DE(vec3 p) {
    float d1 = sdSphere(p);
    float d2 = sdPlane(p);
    if (d1 > d2)
        isFloor = true;
    return min(d1, d2);
}



#preset default
FOV = 0.4
Eye = 11.9485,-11.1753,19.9812
Target = 1.73248,-6.10803,7.96281
Up = 0,0,1
EquiRectangular = false
FocalPlane = 0.7772
Aperture = 0
Gamma = 2.2
ToneMapping = 5
Exposure = 1
Brightness = 1
Contrast = 1
Saturation = 1
GaussianWeight = 1
AntiAliasScale = 1.5
Detail = -3
DetailAO = -0.5
FudgeFactor = 1
MaxRaySteps = 357
BoundingSphere = 100
Dither = 0.5
NormalBackStep = 1
AO = 0,0,0,0.7
Specular = 0.4
SpecularExp = 71.154
SpotLight = 1,1,1,2.4
SpotLightPos = 10,-1.2,10
SpotLightSize = 0
CamLight = 1,1,1,0.43878
CamLightMin = 0
Glow = 1,1,1,0
GlowMax = 14
Fog = 0
Shadow = 0.95455 NotLocked
Sun = 0.61976,0.75831
SunSize = 0.30279
Reflection = 0
BaseColor = 0.435294,0.435294,0.435294
OrbitStrength = 1
X = 0.5,0.6,0.6,0.7
Y = 1,0.6,0,0.4
Z = 0.8,0.78,1,0.5
R = 0.4,0.7,1,0.12
BackgroundColor = 0.6,0.6,0.45
GradientBackground = 0.3
CycleColors = false
Cycles = 1.1
EnableFloor = false
FloorNormal = 0,0,0
FloorHeight = 0
FloorColor = 1,1,1
max_reflections = 155
moveVertexToInf = true
dihedral_A_BCD = 2,2,3
TriangleBCD = 2,3,7
Iterations = 391
edgeSize = 0.0015
edgeSize2 = 0.014
doInvert = false
invPoint = 0.74737,0.87719
#endpreset


#preset 433-334
FOV = 0.4
Eye = 8.41185,-4.58378,11.0765
Target = -0.185318,-1.76319,-2.8022
Up = -0.543931,0.191448,0.375846
EquiRectangular = false
FocalPlane = 0.7772
Aperture = 0
Gamma = 2.2
ToneMapping = 5
Exposure = 1
Brightness = 1
Contrast = 1
Saturation = 1
GaussianWeight = 1
AntiAliasScale = 1.5
Detail = -3
DetailAO = -0.5
FudgeFactor = 1
MaxRaySteps = 357
BoundingSphere = 100
Dither = 0.5
NormalBackStep = 1
AO = 0,0,0,0.7
Specular = 0.4
SpecularExp = 71.154
SpotLight = 1,1,1,2.4
SpotLightPos = 10,-1.2,10
SpotLightSize = 0
CamLight = 1,1,1,0.43878
CamLightMin = 0
Glow = 1,1,1,0
GlowMax = 14
Fog = 0
Shadow = 1 NotLocked
Sun = 0.61976,0.75831
SunSize = 0.30279
Reflection = 0.15054
BaseColor = 0.435294,0.435294,0.435294
OrbitStrength = 1
X = 0.5,0.6,0.6,0.7
Y = 1,0.6,0,0.4
Z = 0.8,0.78,1,0.5
R = 0.4,0.7,1,0.12
BackgroundColor = 0.6,0.6,0.45
GradientBackground = 0.3
CycleColors = false
Cycles = 1.1
EnableFloor = false
FloorNormal = 0,0,0
FloorHeight = 0
FloorColor = 1,1,1
Iterations = 200
orbitDivisor = 126.7
moveVertexToInf = false
dihedral_A_BCD = 4,3,3
TriangleBCD = 3,3,4
edgeSize = 0.0021
edgeSize2 = 0.02347
doInvert = true
invPoint = -0.69132,0.89068
#endpreset
