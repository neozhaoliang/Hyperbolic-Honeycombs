#define USE_EIFFIE_SHADOW
#define MULTI_SAMPLE_AO
#define providesInit
#include "DE-Kn2cr11.frag"

#group BallPacking-Settings
//reflectors s2, s3 and s4 are always one of the 3 types of euclidean Coxeter groups tilings.
uniform int euclideanTriangleType; slider[0,1,2]
// K: Shouldn't min values be 2 ? it seems that a value of 1 means infinity.
//Dihedral angle between reflector s0 and reflectors s2, s3 and s4.
uniform vec3 dihedralAngles0_234; slider[(1,1,1),(3,2,7),(20,20,20)]
//Dihedral angle between reflector s2 and reflectors s2, s3 and s4
uniform vec3 dihedralAngles1_234; slider[(1,1,1),(3,2,7),(20,20,20)]
//Dihedral angle between reflectors s0 and s1
uniform float dihedralAngle0_1; slider[1,4,20]
//Which reflected balls are active. the one associated to s0 is always active (horizontal plane)
uniform vec4 isRealBall; slider[(0,0,0,0),(1,1,0,0),(1,1,1,1)]
// With infinite foldings, 20 iterations are usually sufficient.
//Number of iterations.
uniform int Iterations; slider[0,10,100]

//infinite foldings----------------------------------------------------------------------------------------------------------------//
//- 
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

//----------------------------------------------------------------------------------------------------------------------------------//

// K: Ah! OK! :)
#define inf           1.0 
#define L2(x)         dot(x, x) 

// K: Why not say sqrt(2.) an sqrt(3.) instead ?
#define s2 1.41421356 
#define s3 1.73205081 

float section_height;

float dihedral(float x) {
    return x == inf ? 1. : cos(PI / x);
}

vec3 dihedral(vec3 v) {
    return vec3(dihedral(v.x), dihedral(v.y), dihedral(v.z));
}

// K: Invert member variable removed because not really used here. 
// isRealBall could also be removed. the info could be put inside the sign of r. (?)
struct Ball {
    bool isplane;
    vec3 n;
    float r;
    bool isRealBall;
};


// coclusters are mirror balls, they corresponde to root vectors (space-like)
Ball[5] coclusters;
// clusters are real balls, they corresponde to space-like weight vectors
Ball[5] clusters;

Ball defaultBall() {
    return Ball(false, vec3(0, 0, -1), 0., false);
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
    return Ball(true, n, d, false);
}

Ball from_sphere(vec3 cen, float r) {
    return Ball(false, cen, r, false);
}

bool try_reflect(inout vec3 p,
                 in Ball B,
                 inout float scale,
                 inout vec4 orb) {
    orb = min(orb, vec4(abs(p), dot(p,p)));
    if (B.isplane) {
        float k = dot(vec4(p, 1), vec4(B.n, B.r));
        if (k >= 0.)
            return true;
	else {
            p -= 2. * k  * B.n;
            return false;
	}
    }
    else {
        vec3 cen = B.n;
        float r = B.r;
        vec3 q = p - cen;
        float d2 = dot(q, q);
        float k = (r * r) / d2;
	if ( k < 1.0 )
            return true;
	else {
            scale*=k;
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

    // now we solve the virtual ball B0, this can be either a plane or a sphere
    // this depends on if all entries in dihedralAngles0 are all 2
    //if (dot(dihedralAngles0_234, vec3(1)) == 6.) { // K :well! what happens for (1,2,3) for example ?
    if ( dihedralAngles0_234.x == 2. && dihedralAngles0_234.y == 2. && dihedralAngles0_234.z == 2. ) {
        B0 = from_plane(vec3(0, 0, -1), B1.r*t01);
    }
    else {
        B0 = solveBall(M0, b);
        float r1 = B1.r, r0 = B0.r;
        B0.n.z = sqrt(r0*r0 + r1*r1 + 2.*r0*r1*t01 - L2(B1.n.xy - B0.n.xy));
    }
    coclusters = Ball[5] (B0, B1, B2, B3, B4);

    section_height = B0.isplane ? 2.*B0.r : B0.n.z;

    //now we process the real balls
    for (int k = 0; k < 5; k++)
        clusters[k] = defaultBall();
  
    clusters[1] = from_plane(vec3(0, 0, -1.), B0.n.z);
    clusters[1].isRealBall = isRealBall.x;

    clusters[2] = solveBall(C, B0, B1);
    clusters[2].isRealBall = isRealBall.y;
    
    clusters[3] = solveBall(A, B0, B1);
    clusters[3].isRealBall = isRealBall.z;

    clusters[4] = solveBall(B, B0, B1);
    clusters[4].isRealBall = isRealBall.w;
}
 
float dist2balls( vec3 p, float scale) {
    float d = 1e5; //Ball 0
    for (int j = 1; j < 5; j++) 
        if (clusters[j].isRealBall) 
          d = min(abs(sdistanceToBall(p, clusters[j])), d);
    return d/scale;
}


void EuclideanFold(inout vec3 p) {
    if( euclideanTriangleType == 0)
        fold236(p.xy);
    else if( euclideanTriangleType == 1)
        fold244(p.xy);
    else
        fold333(p.xy);
}

float map(inout vec3 p, inout float scale, inout vec4 orb) {
    p.z += 1.;
    scale /= dot(p,p);
    p *= 2. / dot(p,p);
    p.z -= 1.;

    for (int i = 0; i < Iterations ; i++) {
	vec3 ap = p;
        EuclideanFold(p);
        for (int k = 0; k < 2; k++) 
            try_reflect(p, coclusters[k], scale, orb);
        if (all(not(bvec3(p-ap)))) break;
    }

    return dist2balls(p, scale);
}



float DE(vec3 p) {

    float DEfactor=1.;
    //orbiTrap is not initialized by default. That's what was producing the strange haze and colors.
    orbitTrap = vec4(1);

    float d = map(p, DEfactor, orbitTrap);
    return 0.25 * d;
}

















#preset default
AutoFocus = false
FOV = 0.685022
Eye = 1.89674,0.0152457,0.0385651
Target = -0.103173,-0.00352906,0.0362396
UpLock = false
Up = -0.00117871,0.00169458,0.999998
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
BloomStrong = 1
DepthToAlpha = false
Detail = -3.5
RefineSteps = 3
FudgeFactor = 0.6
MaxRaySteps = 500
MaxDistance = 10
Dither = 0.5
NormalBackStep = 6.5875
DetailAO = -1.96825
coneApertureAO = 0.56566
maxIterAO = 15
FudgeAO = 1
AO_ambient = 1
AO_camlight = 0
AO_pointlight = 0
AoCorrect = 0
Specular = 0.1
SpecularExp = 326.435
CamLight = 0.364706,0.364706,0.364706,0.45
AmbiantLight = 0.709804,0.709804,0.709804,0.57338
Reflection = 0.192157,0.192157,0.192157
ReflectionsNumber = 0 Locked
SpotGlow = false
SpotLight = 1,1,1,1
LightPos = 6,8,3
LightSize = 0
LightFallOff = 0
LightGlowRad = 0
LightGlowExp = 0
HardShadow = 1
ShadowSoft = 0
ShadowBlur = 0.2
perf = false
SSS = false
sss1 = 0.1
sss2 = 0.5
BaseColor = 0.898039,0.937255,0.976471
OrbitStrength = 0.70833
X = 0.666667,1,0.498039,1
Y = 1,0.333333,0,1
Z = 0.333333,0,1,1
R = 0.333333,0.666667,1,1
BackgroundColor = 0.270588,0.403922,0.6
GradientBackground = 0
CycleColors = true
Cycles = 4
EnableFloor = false
FloorNormal = 0,0,1
FloorHeight = 0
FloorColor = 1,1,1
HF_Fallof = 0.187344
HF_Const = 0
HF_Intensity = 0
HF_Dir = 0,0,1
HF_Offset = 0
HF_Color = 0.564706,0.752941,0.878431,0
HF_Scatter = 10
HF_Anisotropy = 0.168627,0.168627,0.168627
HF_FogIter = 1
HF_CastShadow = false
EnCloudsDir = false
Clouds_Dir = 0,0,1
CloudScale = 1
CloudFlatness = 0
CloudTops = 1
CloudBase = -1
CloudDensity = 0.484136
CloudRoughness = 1
CloudContrast = 1
CloudColor = 0.65,0.68,0.7
CloudColor2 = 0.07,0.17,0.24
SunLightColor = 0.7,0.5,0.3
Cloudvar1 = 0.99
Cloudvar2 = 0.99
CloudIter = 5
CloudBgMix = 1
euclideanTriangleType = 1
dihedralAngles0_234 = 3,2,4
dihedralAngles1_234 = 3,3,3
dihedralAngle0_1 = 5
isRealBall = 1,0,0,1
Iterations = 100
#endpreset

