#include "parser.h"
#include "ppm.h"
#include <cmath>

typedef unsigned char RGB[3];

struct Ray {
    parser::Vec3f origin;
    parser::Vec3f direction;
};

parser::Vec3f multS(parser::Vec3f a,float s)
{
    parser::Vec3f result{};
    result.x = a.x*s;
    result.y = a.y*s;
    result.z = a.z*s;
    return result;
}

float dot(parser::Vec3f a,parser::Vec3f b)
{
    return a.x*b.x+a.y*b.y+a.z*b.z;
}

parser::Vec3f normalize(parser::Vec3f a)
{
    return multS(a,1.0/sqrt(dot(a,a)));
}

parser::Vec3f hadamard(parser::Vec3f a,parser::Vec3f b)
{
    parser::Vec3f result{};
    result.x = a.x*b.x;
    result.y = a.y*b.y;
    result.z = a.z*b.z;
    return result;
}

parser::Vec3f add(parser::Vec3f a, parser::Vec3f b)
{
    parser::Vec3f result{};
    result.x = a.x+b.x;
    result.y = a.y+b.y;
    result.z = a.z+b.z;
    return result;
}

parser::Vec3f sub(parser::Vec3f a, parser::Vec3f b)
{
    parser::Vec3f result{};
    result.x = a.x-b.x;
    result.y = a.y-b.y;
    result.z = a.z-b.z;
    return result;
}
parser::Vec3f CrossProduct(parser::Vec3f a, parser::Vec3f b) {
    parser::Vec3f result{};
    result.x = a.y * b.z - a.z * b.y;
    result.y = a.z * b.x - a.x * b.z;
    result.z = a.x * b.y - a.y * b.x;
    return result;
}

parser::Vec3f NormalOfFace(parser::Face face, std::vector<parser::Vec3f> vertices) {
    parser::Vec3f v0 = vertices[face.v0_id];
    parser::Vec3f v1 = vertices[face.v1_id];
    parser::Vec3f v2 = vertices[face.v2_id];
    parser::Vec3f v01 = sub(v1, v0);
    parser::Vec3f v02 = sub(v2, v0);
    return normalize(CrossProduct(v01, v02));
}

Ray GenerateRay(int i, int j, const parser::Camera& camera, int width, int height)
{
    Ray result{};
    float su,sv;
    parser::Vec3f m{},q{},s{};

    auto left = camera.near_plane.x;
    auto right = camera.near_plane.y;
    auto bottom = camera.near_plane.z;
    auto top = camera.near_plane.w;

    su = (float) (i+0.5)*(right - left)/(float) width;
    sv = (float) (j+0.5)*(top - bottom)/(float) height;
    parser::Vec3f v = camera.up;
    parser::Vec3f w = multS(camera.gaze,-1);
    parser::Vec3f u = CrossProduct(v,w);

    m = add(camera.position,multS(camera.gaze,camera.near_distance));
    q = add(m, add(multS(u,left),multS(v,top)));
    s = add(q,add(multS(u,su),multS(v,-sv)));
    result.origin = camera.position;
    result.direction = add(s,multS(camera.position,-1));
    return result;
}



float intersectSphere (parser::Vec3f sphereCenter, float sphereRadius, Ray ray) {
    float A,B,C; //constants for the quadratic equation
    float delta;
    float t,t1,t2;

    parser::Vec3f c = sphereCenter;

    C = (ray.origin.x-c.x)*(ray.origin.x-c.x)+(ray.origin.y-c.y)*(ray.origin.y-c.y)+(ray.origin.z-c.z)*(ray.origin.z-c.z)-sphereRadius*sphereRadius;

    B = 2*ray.direction.x*(ray.origin.x-c.x)+2*ray.direction.y*(ray.origin.y-c.y)+2*ray.direction.z*(ray.origin.z-c.z);

    A = ray.direction.x*ray.direction.x+ray.direction.y*ray.direction.y+ray.direction.z*ray.direction.z;

    delta = B*B-4*A*C;

    if (delta<0) return -1;
    else if (delta==0)
    {
        t = -B / (2*A);
    }
    else
    {
        delta = sqrt(delta);
        A = 2*A;
        t1 = (-B + delta) / A;
        t2 = (-B - delta) / A;

        t1<t2 ? t=t1 : t=t2;
    }

    return t;
}

//intersect mesh function by using barycentric coordinates
parser::Vec3f intersectFace (parser::Face face, Ray ray, parser::Scene scene) {
    parser::Vec3f v0 = scene.vertex_data[face.v0_id - 1];
    parser::Vec3f v1 = scene.vertex_data[face.v1_id - 1];
    parser::Vec3f v2 = scene.vertex_data[face.v2_id - 1];
    float beta, gamma;
    float t;

    //denom_det is the determinant of the A matrix in the slide
    float denom_det = (v0.x - v1.x)*((ray.direction.z * (v0.y - v2.y)) - (v0.z - v2.z * ray.direction.y));
    denom_det += (v0.y - v1.y)*(ray.direction.x * (v0.z - v2.z) - (v0.x - v2.x) * ray.direction.z);
    denom_det += (v0.z - v1.z)*(ray.direction.y * (v0.x - v2.x) - (v0.y - v2.y) * ray.direction.x);

    float beta_det = (v0.x-ray.origin.x) * ((v0.y - v2.y) * ray.direction.z - (v0.z - v2.z) * ray.direction.y);
    beta_det += (v0.y - ray.origin.y) * ((v0.z - v2.z) * ray.direction.x - (v0.x - v2.x) * ray.direction.z);
    beta_det += (v0.z - ray.origin.z) * ((v0.x - v2.x) * ray.direction.y - (v0.y - v2.y) * ray.direction.x);

    float gamma_det = (v0.x - v1.x) * ((v0.y - ray.origin.y) * ray.direction.z - (v0.z - ray.origin.z) * ray.direction.y);
    gamma_det += (v0.y - v1.y) * ((v0.z - ray.origin.z) * ray.direction.x - (v0.x - ray.origin.x) * ray.direction.z);
    gamma_det += (v0.z - v1.z) * ((v0.x - ray.origin.x) * ray.direction.y - (v0.y - ray.origin.y) * ray.direction.x);

    float t_det = (v0.x - v1.x) * ((v0.y - v2.y) * (v0.z - ray.origin.z) - (v0.z - v2.z) * (v0.y - ray.origin.y));
    t_det += (v0.y - v1.y) * ((v0.z - v2.z) * (v0.x - ray.origin.x) - (v0.x - v2.x) * (v0.z - ray.origin.z));
    t_det += (v0.z - v1.z) * ((v0.x - v2.x) * (v0.y - ray.origin.y) - (v0.y - v2.y) * (v0.x - ray.origin.x));

    beta = beta_det/denom_det;
    gamma = gamma_det/denom_det;
    t = t_det/denom_det;

    return {t, beta, gamma};
}

parser::Vec3f ComputeColor(Ray ray, parser::Scene scene)
{
    float t = 1000000;
    parser::Vec3f color{};
    parser::Sphere intersectedSphere{};
    parser::Mesh intersectedMesh{};
    for (auto sphere : scene.spheres) {
        auto sphereCenter = scene.vertex_data[sphere.center_vertex_id-1];
        float temp = intersectSphere(sphereCenter, sphere.radius, ray);
        if (temp > 0 && temp < t) {
            t = temp;
            intersectedSphere = sphere;
        }
    }
    for (auto mesh : scene.meshes) {
        for (auto face : mesh.faces) {
            auto temp = intersectFace(face, ray, scene);
            // function returns a Vec3f where temp.x is t, temp.y is beta, temp.z is gamma
            if (temp.x > 0 && temp.x < t && temp.y >= 0 && temp.z >= 0 && (temp.y + temp.z) <= 1) {
                t = temp.x;
                intersectedMesh = mesh;
            }
        }
    }
    if (t == 1000000) {
        return color;
    }
    if (intersectedSphere.radius != 0) {
        auto P = add(ray.origin,multS(ray.direction,t));
        auto L = normalize(add(scene.point_lights[0].position,multS(P,-1)));
        auto N = normalize(add(P,multS(scene.vertex_data[intersectedSphere.center_vertex_id-1],-1)));
        auto W = normalize(add(ray.origin,multS(P,-1)));
        auto H = normalize(add(L,W));
        auto lightDistance = sqrt((scene.point_lights[0].position.x-P.x)*(scene.point_lights[0].position.x-P.x)+(scene.point_lights[0].position.y-P.y)*(scene.point_lights[0].position.y-P.y)+(scene.point_lights[0].position.z-P.z)*(scene.point_lights[0].position.z-P.z));
        auto receivedIrradiance = multS(scene.point_lights[0].intensity,1/(lightDistance*lightDistance));
        auto cosAlphaSpecular = fmax(dot(N,H),0.0);
        auto cosAlphaDiffuse = fmax(dot(N,L),0.0);
        auto specular = multS(scene.materials[intersectedSphere.material_id-1].specular,pow(cosAlphaSpecular, scene.materials[intersectedSphere.material_id-1].phong_exponent));
        auto diffuse = multS(scene.materials[intersectedSphere.material_id-1].diffuse,cosAlphaDiffuse);
        color = add(color, add(diffuse, specular));
        color = hadamard(color, receivedIrradiance);
        color = add(color, hadamard(scene.materials[intersectedSphere.material_id-1].ambient, scene.ambient_light));
    }
    if (!intersectedMesh.faces.empty()) {
        auto P = add(ray.origin,multS(ray.direction,t));
        auto L = normalize(add(scene.point_lights[0].position,multS(P,-1)));
        auto N = NormalOfFace(intersectedMesh.faces[0], scene.vertex_data);
        auto W = normalize(add(ray.origin,multS(P,-1)));
        auto H = normalize(add(L,W));
        auto lightDistance = sqrt((scene.point_lights[0].position.x-P.x)*(scene.point_lights[0].position.x-P.x)+(scene.point_lights[0].position.y-P.y)*(scene.point_lights[0].position.y-P.y)+(scene.point_lights[0].position.z-P.z)*(scene.point_lights[0].position.z-P.z));
        auto receivedIrradiance = multS(scene.point_lights[0].intensity,1/(lightDistance*lightDistance));
        auto cosAlphaSpecular = fmax(dot(N,H),0.0);
        auto cosAlphaDiffuse = fmax(dot(N,L),0.0);
        auto specular = multS(scene.materials[intersectedMesh.material_id-1].specular,pow(cosAlphaSpecular, scene.materials[intersectedMesh.material_id-1].phong_exponent));
        auto diffuse = multS(scene.materials[intersectedMesh.material_id-1].diffuse,cosAlphaDiffuse);
        color = add(color, add(diffuse, specular));
        color = hadamard(color, receivedIrradiance);
        color = add(color, hadamard(scene.materials[intersectedMesh.material_id-1].ambient, scene.ambient_light));
    }
    return color;
}

int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file
    parser::Scene scene;

    scene.loadFromXml(argv[1]);
    auto spheres = scene.spheres;

    int width = scene.cameras[0].image_width;
    int height = scene.cameras[0].image_height;

    auto backroundColor = scene.background_color;

    auto* image = new unsigned char [width * height * 3]; // 3 channels

    for (int i = 0; i < width * height; ++i) //initialize image to backround color
    {
        image[i*3] = backroundColor.x;
        image[i*3+1] = backroundColor.y;
        image[i*3+2] = backroundColor.z;
    }

    for (int j = 0; j < height; j++) { //main loop
        for (int i = 0; i < width; i++) {
            Ray ray{};
            ray = GenerateRay(i, j, scene.cameras[0], width, height);
            auto color = ComputeColor(ray, scene);
            image[(j * width + i) * 3] = color.x > 255 ? 255 : round(color.x);
            image[(j * width + i) * 3 + 1] = color.y > 255 ? 255 : round(color.y);
            image[(j * width + i) * 3 + 2] = color.z > 255 ? 255 : round(color.z);
        }
    }

    write_ppm("test.ppm", image, width, height);

    delete[] image;

    return 0;
}
