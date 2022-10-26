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

parser::Vec3f add(parser::Vec3f a, parser::Vec3f b)
{
    parser::Vec3f result{};
    result.x = a.x+b.x;
    result.y = a.y+b.y;
    result.z = a.z+b.z;
    return result;
}
parser::Vec3f CrossProduct(parser::Vec3f a, parser::Vec3f b) {
    parser::Vec3f result{};
    result.x = a.y * b.z - a.z * b.y;
    result.y = a.z * b.x - a.x * b.z;
    result.z = a.x * b.y - a.y * b.x;
    return result;
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
            auto v1 = scene.vertex_data[face.v0_id-1];
            auto v2 = scene.vertex_data[face.v1_id-1];
            auto v3 = scene.vertex_data[face.v2_id-1];
            auto e1 = add(v2,multS(v1,-1));
            auto e2 = add(v3,multS(v1,-1));
            auto pvec = CrossProduct(ray.direction,e2);
            auto det = dot(e1,pvec);
            if (det > -0.000001 && det < 0.000001) continue;
            auto inv_det = 1.0/det;
            auto tvec = add(ray.origin,multS(v1,-1));
            auto u = dot(tvec,pvec)*inv_det;
            if (u < 0 || u > 1) continue;
            auto qvec = CrossProduct(tvec,e1);
            auto v = dot(ray.direction,qvec)*inv_det;
            if (v < 0 || u + v > 1) continue;
            auto temp = dot(e2,qvec)*inv_det;
            if (temp > 0 && temp < t) {
                t = temp;
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
        auto receivedIrradiance = normalize(multS(scene.point_lights[0].intensity,1/(lightDistance*lightDistance)));
        auto cosAlphaSpecular = fmax(dot(N,H),0.0);
        auto cosAlphaDiffuse = fmax(dot(N,L),0.0);
        auto specular = multS(scene.materials[intersectedSphere.material_id-1].specular,pow(cosAlphaSpecular, scene.materials[intersectedSphere.material_id-1].phong_exponent));
        auto diffuse = multS(scene.materials[intersectedSphere.material_id-1].diffuse,cosAlphaDiffuse);
        color = add(color, add(diffuse, specular));
        color.x += color.x*receivedIrradiance.x + (scene.ambient_light.x*scene.materials[intersectedSphere.material_id-1].ambient.x)/255;
        color.y += color.y*receivedIrradiance.y + (scene.ambient_light.y*scene.materials[intersectedSphere.material_id-1].ambient.y)/255;
        color.z += color.z*receivedIrradiance.z + (scene.ambient_light.z*scene.materials[intersectedSphere.material_id-1].ambient.z)/255;
    }
    if (!intersectedMesh.faces.empty()) {
        auto meshColor = scene.materials[intersectedMesh.material_id-1].diffuse;
        color = add(CrossProduct(scene.materials[intersectedMesh.material_id-1].ambient,scene.ambient_light), meshColor);
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
            image[(j * width + i) * 3] = color.x * 255 > 255 ? 255 : color.x * 255;
            image[(j * width + i) * 3 + 1] = color.y * 255 > 255 ? 255 : color.y * 255;
            image[(j * width + i) * 3 + 2] = color.z * 255 > 255 ? 255 : color.z * 255;
        }
    }

    write_ppm("test.ppm", image, width, height);

    delete[] image;

    return 0;
}
