#version 100

// https://iquilezles.org/www/articles/distfunctions/distfunctions.htm
// http://jamie-wong.com/2016/07/15/ray-marching-signed-distance-functions/

float sphereSDF(vec3 p, float radius)
{
    return length(p)-radius;
}


float boxSDF(vec3 p, vec3 b)
{
    vec3 q = abs(p) - b;
    return length(max(q, 0.0)) + min(max(q.x, max(q.y, q.z)), 0.0);
}


float roundBoxSDF(vec3 p, vec3 b, float r)
{
    vec3 q = abs(p) - b;
    return length(max(q, 0.0)) + min(max(q.x, max(q.y, q.z)), 0.0) - r;
}


float torusSDF(vec3 p, vec2 t)
{
    vec2 q = vec2(length(p.xz)-t.x, p.y);
    return length(q)-t.y;
}


float boundingBoxSDF(vec3 p, vec3 b, float e)
{
    p = abs(p)-b;
    vec3 q = abs(p+e)-e;
    return min(min(
    length(max(vec3(p.x, q.y, q.z), 0.0))+min(max(p.x, max(q.y, q.z)), 0.0),
    length(max(vec3(q.x, p.y, q.z), 0.0))+min(max(q.x, max(p.y, q.z)), 0.0)),
    length(max(vec3(q.x, q.y, p.z), 0.0))+min(max(q.x, max(q.y, p.z)), 0.0));
}


float cappedTorusSDF(vec3 p, vec2 sc, float ra, float rb)
{
    p.x = abs(p.x);
    float k = (sc.y*p.x>sc.x*p.y) ? dot(p.xy, sc) : length(p.xy);
    return sqrt(dot(p, p) + ra*ra - 2.0*ra*k) - rb;
}


float linkSDF(vec3 p, float le, float r1, float r2)
{
    vec3 q = vec3(p.x, max(abs(p.y)-le, 0.0), p.z);
    return length(vec2(length(q.xy)-r1, q.z)) - r2;
}


float cylinderSDF(vec3 p, vec3 c)
{
    return length(p.xz-c.xy)-c.z;
}


float coneSDF(vec3 p, vec2 c, float h)
{
    // c is the sin/cos of the angle, h is height
    // Alternatively pass q instead of (c,h),
    // which is the point at the base in 2D
    vec2 q = h*vec2(c.x/c.y, -1.0);

    vec2 w = vec2(length(p.xz), p.y);
    vec2 a = w - q*clamp(dot(w, q)/dot(q, q), 0.0, 1.0);
    vec2 b = w - q*vec2(clamp(w.x/q.x, 0.0, 1.0), 1.0);
    float k = sign(q.y);
    float d = min(dot(a, a), dot(b, b));
    float s = max(k*(w.x*q.y-w.y*q.x), k*(w.y-q.y));
    return sqrt(d)*sign(s);
}

float simpleConeSDF(vec3 p, vec2 c, float h)
{
    float q = length(p.xz);
    return max(dot(c.xy, vec2(q, p.y)), -h-p.y);
}


float infiniteConeSDF(vec3 p, vec2 c)
{
    // c is the sin/cos of the angle
    vec2 q = vec2(length(p.xz), -p.y);
    float d = length(q-c*max(dot(q, c), 0.0));
    return d * ((q.x*c.y-q.y*c.x<0.0)?-1.0:1.0);
}


float planeSDF(vec3 p, vec3 n, float h) //n is the normal vector of the plane and h is the distance from origin.
{
    // n must be normalized
    return dot(p, n) + h;
}


float hexPrismSDF(vec3 p, vec2 h)
{
    const vec3 k = vec3(-0.8660254, 0.5, 0.57735);
    p = abs(p);
    p.xy -= 2.0*min(dot(k.xy, p.xy), 0.0)*k.xy;
    vec2 d = vec2(
    length(p.xy-vec2(clamp(p.x, -k.z*h.x, k.z*h.x), h.x))*sign(p.y-h.x),
    p.z-h.y);
    return min(max(d.x, d.y), 0.0) + length(max(d, 0.0));
}


float triPrismSDF(vec3 p, vec2 h)
{
    vec3 q = abs(p);
    return max(q.z-h.y, max(q.x*0.866025+p.y*0.5, -p.y)-h.x*0.5);
}


float capsuleSDF(vec3 p, vec3 a, vec3 b, float r)
{
    vec3 pa = p - a, ba = b - a;
    float h = clamp(dot(pa, ba)/dot(ba, ba), 0.0, 1.0);
    return length(pa - ba*h) - r;
}


float verticalCapsuleSDF(vec3 p, float h, float r)
{
    p.y -= clamp(p.y, 0.0, h);
    return length(p) - r;
}


float cappedCylinderSDF(vec3 p, vec3 a, vec3 b, float r)
{
    vec3  ba = b - a;
    vec3  pa = p - a;
    float baba = dot(ba, ba);
    float paba = dot(pa, ba);
    float x = length(pa*baba-ba*paba) - r*baba;
    float y = abs(paba-baba*0.5)-baba*0.5;
    float x2 = x*x;
    float y2 = y*y*baba;
    float d = (max(x, y)<0.0)?-min(x2, y2):(((x>0.0)?x2:0.0)+((y>0.0)?y2:0.0));
    return sign(d)*sqrt(abs(d))/baba;
}


float verticalCappedCylinderSDF(vec3 p, float h, float r)
{
    vec2 d = abs(vec2(length(p.xz), p.y)) - vec2(h, r);
    return min(max(d.x, d.y), 0.0) + length(max(d, 0.0));
}


float roundedCylinderSDF(vec3 p, float ra, float rb, float h)
{
    vec2 d = vec2(length(p.xz)-2.0*ra+rb, abs(p.y) - h);
    return min(max(d.x, d.y), 0.0) + length(max(d, 0.0)) - rb;
}


float cappedConeSDF(vec3 p, vec3 a, vec3 b, float ra, float rb)
{
    float rba  = rb-ra;
    float baba = dot(b-a, b-a);
    float papa = dot(p-a, p-a);
    float paba = dot(p-a, b-a)/baba;
    float x = sqrt(papa - paba*paba*baba);
    float cax = max(0.0, x-((paba<0.5)?ra:rb));
    float cay = abs(paba-0.5)-0.5;
    float k = rba*rba + baba;
    float f = clamp((rba*(x-ra)+paba*baba)/k, 0.0, 1.0);
    float cbx = x-ra - f*rba;
    float cby = paba - f;
    float s = (cbx < 0.0 && cay < 0.0) ? -1.0 : 1.0;
    return s*sqrt(min(cax*cax + cay*cay*baba,
    cbx*cbx + cby*cby*baba));
}


float verticalCappedConeSDF(vec3 p, float h, float r1, float r2)
{
    vec2 q = vec2(length(p.xz), p.y);
    vec2 k1 = vec2(r2, h);
    vec2 k2 = vec2(r2-r1, 2.0*h);
    vec2 ca = vec2(q.x-min(q.x, (q.y<0.0)?r1:r2), abs(q.y)-h);
    vec2 cb = q - k1 + k2*clamp(dot(k1-q, k2)/dot2(k2), 0.0, 1.0);
    float s = (cb.x<0.0 && ca.y<0.0) ? -1.0 : 1.0;
    return s*sqrt(min(dot2(ca), dot2(cb)));
}


float solidAngleSDF(vec3 p, vec2 c, float ra)
{
    // c is the sin/cos of the angle
    vec2 q = vec2(length(p.xz), p.y);
    float l = length(q) - ra;
    float m = length(q - c*clamp(dot(q, c), 0.0, ra));
    return max(l, m*sign(c.y*q.x-c.x*q.y));
}


float verticalRoundConeSDF(vec3 p, float r1, float r2, float h)
{
    vec2 q = vec2(length(p.xz), p.y);

    float b = (r1-r2)/h;
    float a = sqrt(1.0-b*b);
    float k = dot(q, vec2(-b, a));

    if (k < 0.0) return length(q) - r1;
    if (k > a*h) return length(q-vec2(0.0, h)) - r2;

    return dot(q, vec2(a, b)) - r1;
}


float roundConeSDF(vec3 p, vec3 a, vec3 b, float r1, float r2)
{
    // sampling independent computations (only depend on shape)
    vec3  ba = b - a;
    float l2 = dot(ba, ba);
    float rr = r1 - r2;
    float a2 = l2 - rr*rr;
    float il2 = 1.0/l2;

    // sampling dependant computations
    vec3 pa = p - a;
    float y = dot(pa, ba);
    float z = y - l2;
    float x2 = dot2(pa*l2 - ba*y);
    float y2 = y*y*l2;
    float z2 = z*z*l2;

    // single square root!
    float k = sign(rr)*rr*rr*x2;
    if (sign(z)*a2*z2 > k) return sqrt(x2 + z2)        *il2 - r2;
    if (sign(y)*a2*y2 < k) return sqrt(x2 + y2)        *il2 - r1;
    return (sqrt(x2*a2*il2)+y*rr)*il2 - r1;
}


float ellipsoidSDF(vec3 p, vec3 r)
{
    float k0 = length(p/r);
    float k1 = length(p/(r*r));
    return k0*(k0-1.0)/k1;
}

float rhombusSDF(vec3 p, float la, float lb, float h, float ra)
{
    p = abs(p);
    vec2 b = vec2(la, lb);
    float f = clamp((ndot(b, b-2.0*p.xz))/dot(b, b), -1.0, 1.0);
    vec2 q = vec2(length(p.xz-0.5*b*vec2(1.0-f, 1.0+f))*sign(p.x*b.y+p.z*b.x-b.x*b.y)-ra, p.y-h);
    return min(max(q.x, q.y), 0.0) + length(max(q, 0.0));
}


float octahedronSDF(vec3 p, float s)
{
    p = abs(p);
    float m = p.x+p.y+p.z-s;
    vec3 q;
    if (3.0*p.x < m) q = p.xyz;
    else if (3.0*p.y < m) q = p.yzx;
    else if (3.0*p.z < m) q = p.zxy;
    else return m*0.57735027;

    float k = clamp(0.5*(q.z-q.y+s), 0.0, s);
    return length(vec3(q.x, q.y-s+k, q.z-k));
}


float simpleOctahedronSDF(vec3 p, float s)
{
    p = abs(p);
    return (p.x+p.y+p.z-s)*0.57735027;
}


float PyramidSDF(vec3 p, float h)
{
    float m2 = h*h + 0.25;

    p.xz = abs(p.xz);
    p.xz = (p.z>p.x) ? p.zx : p.xz;
    p.xz -= 0.5;

    vec3 q = vec3(p.z, h*p.y - 0.5*p.x, h*p.x + 0.5*p.y);

    float s = max(-q.x, 0.0);
    float t = clamp((q.y-0.5*p.z)/(m2+0.25), 0.0, 1.0);

    float a = m2*(q.x+s)*(q.x+s) + q.y*q.y;
    float b = m2*(q.x+0.5*t)*(q.x+0.5*t) + (q.y-m2*t)*(q.y-m2*t);

    float d2 = min(q.y, -q.x*m2-q.y*0.5) > 0.0 ? 0.0 : min(a, b);

    return sqrt((d2+q.z*q.z)/m2) * sign(max(q.z, -p.y));
}


float triangleUDF(vec3 p, vec3 a, vec3 b, vec3 c)
{
    vec3 ba = b - a; vec3 pa = p - a;
    vec3 cb = c - b; vec3 pb = p - b;
    vec3 ac = a - c; vec3 pc = p - c;
    vec3 nor = cross(ba, ac);

    return sqrt(
    (sign(dot(cross(ba, nor), pa)) +
    sign(dot(cross(cb, nor), pb)) +
    sign(dot(cross(ac, nor), pc))<2.0)
    ?
    min(min(
    dot2(ba*clamp(dot(ba, pa)/dot2(ba), 0.0, 1.0)-pa),
    dot2(cb*clamp(dot(cb, pb)/dot2(cb), 0.0, 1.0)-pb)),
    dot2(ac*clamp(dot(ac, pc)/dot2(ac), 0.0, 1.0)-pc))
    :
    dot(nor, pa)*dot(nor, pa)/dot2(nor));
}


float quadUDF(vec3 p, vec3 a, vec3 b, vec3 c, vec3 d)
{
    vec3 ba = b - a; vec3 pa = p - a;
    vec3 cb = c - b; vec3 pb = p - b;
    vec3 dc = d - c; vec3 pc = p - c;
    vec3 ad = a - d; vec3 pd = p - d;
    vec3 nor = cross(ba, ad);

    return sqrt(
    (sign(dot(cross(ba, nor), pa)) +
    sign(dot(cross(cb, nor), pb)) +
    sign(dot(cross(dc, nor), pc)) +
    sign(dot(cross(ad, nor), pd))<3.0)
    ?
    min(min(min(
    dot2(ba*clamp(dot(ba, pa)/dot2(ba), 0.0, 1.0)-pa),
    dot2(cb*clamp(dot(cb, pb)/dot2(cb), 0.0, 1.0)-pb)),
    dot2(dc*clamp(dot(dc, pc)/dot2(dc), 0.0, 1.0)-pc)),
    dot2(ad*clamp(dot(ad, pd)/dot2(ad), 0.0, 1.0)-pd))
    :
    dot(nor, pa)*dot(nor, pa)/dot2(nor));
}


//-------------------------------------------------------------------------------------------------------------
//Operation
//-------------------------------------------------------------------------------------------------------------


float opElongatefast(SDF primitive, vec3 p, vec3 h)// Function pointer not supported in GLSL. Cannot applied directly
{
    vec3 q = p - clamp(p, -h, h);
    return primitive(q);
}

float opElongate(SDF primitive, vec3 p, vec3 h)
{
    vec3 q = abs(p)-h;
    return primitive(max(q, 0.0)) + min(max(q.x, max(q.y, q.z)), 0.0);
}

float opRound(SDF primitive, float rad)
{
    return primitive(p) - rad;
}

float opOnion(float sdf, float thickness)
{
    return abs(sdf)-thickness;
}

float opExtrusion(vec3 p, sdf2d primitive, float h)
{
    float d = primitive(p.xy);
    vec2 w = vec2(d, abs(p.z) - h);
    return min(max(w.x, w.y), 0.0) + length(max(w, 0.0));
}

float opRevolution(vec3 p, sdf2d primitive, float o)
{
    vec2 q = vec2(length(p.xz) - o, p.y);
    return primitive(q);
}

float opUnion(float d1, float d2) { return min(d1, d2); }

float opSubtraction(float d1, float d2) { return max(-d1, d2); }

float opIntersection(float d1, float d2) { return max(d1, d2); }

float opSmoothUnion(float d1, float d2, float k)
{
    float h = max(k-abs(d1-d2), 0.0);
    return min(d1, d2) - h*h*0.25/k;
}

float opPolynomialSmoothUnion(float d1, float d2, float k)//k = 0.1
{
    float h = clamp(0.5 + 0.5*(d2-d1)/k, 0.0, 1.0);
    return mix(d2, d1, h) - k*h*(1.0-h);
}

float opExponentialSmoothUnion(float a, float b, float k)//k = 32
{
    float res = exp2(-k*a) + exp2(-k*b);
    return -log2(res)/k;
}

float opPowerSmoothUnion(float a, float b, float k)//k = 8
{
    a = pow(a, k); b = pow(b, k);
    return pow((a*b)/(a+b), 1.0/k);
}

float opSmoothSubtraction(float d1, float d2, float k)
{
    float h = max(k-abs(-d1-d2), 0.0);
    return max(-d1, d2) + h*h*0.25/k;
}

float opPolynomialSmoothSubtraction(float d1, float d2, float k)
{
    float h = clamp(0.5 - 0.5*(d2+d1)/k, 0.0, 1.0);
    return mix(d2, -d1, h) + k*h*(1.0-h);
}

float opSmoothIntersection(float d1, float d2, float k)
{
    float h = max(k-abs(d1-d2), 0.0);
    return max(d1, d2) + h*h*0.25/k;
}

float opPolynomialSmoothIntersection(float d1, float d2, float k)
{
    float h = clamp(0.5 - 0.5*(d2-d1)/k, 0.0, 1.0);
    return mix(d2, d1, h) + k*h*(1.0-h);
}

vec3 opTransform(vec3 p, transform T, sdf3d primitive)
{
    return primitive(invert(T)*p);
}

float opScale(vec3 p, float s, sdf3d primitive)
{
    return primitive(p/s)*s;
}

float opSymX(vec3 p, sdf3d primitive)
{
    p.x = abs(p.x);
    return primitive(p);
}

float opSymXZ(vec3 p, sdf3d primitive)
{
    p.xz = abs(p.xz);
    return primitive(p);
}

float opRep(vec3 p, vec3 c, sdf3d primitive)//Infinite Repetition
{
    vec3 q = mod(p+0.5*c, c)-0.5*c;
    return primitive(q);
}

vec3 opRepLim(in vec3 p, in float c, in vec3 l, in sdf3d primitive)
{
    vec3 q = p-c*clamp(round(p/c), -l, l);
    return primitive(q);
}

vec3 estimateNormal(vec3 p) {
    return normalize(vec3(
    sceneSDF(vec3(p.x + EPSILON, p.y, p.z)) - sceneSDF(vec3(p.x - EPSILON, p.y, p.z)),
    sceneSDF(vec3(p.x, p.y + EPSILON, p.z)) - sceneSDF(vec3(p.x, p.y - EPSILON, p.z)),
    sceneSDF(vec3(p.x, p.y, p.z  + EPSILON)) - sceneSDF(vec3(p.x, p.y, p.z - EPSILON))
    ));
}

void main() {
    gl_FragColor = vec4(1.0, 0.0, 0.0, 1.0);
}
